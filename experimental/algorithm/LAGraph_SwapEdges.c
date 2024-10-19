//------------------------------------------------------------------------------
// LAGraph_SwapEdges: Randomly Swaps edges in a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2019-2022 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Timothy A. Davis, Texas A&M University

//------------------------------------------------------------------------------

#define FREE_LOOP                               \
{                                               \
    GrB_free (&pairs) ;                         \
    GrB_free (&M) ;                             \
    GrB_free (&M_fours) ;                       \
    GrB_free (&pairs_4s) ;                      \
    GrB_free (&M_t) ;                           \
    GrB_free (&interf) ;                        \
    GrB_free (&M_i) ;                           \
    GrB_free (&pairs_i) ;                       \
    GrB_free (&M_1) ;                           \
    GrB_free (&M_2) ;                           \
    GrB_free (&exists_1) ;                      \
    GrB_free (&exists_2) ;                      \
    GrB_free (E_arranged) ;                     \
    GrB_free (E_arranged + 1) ;                 \
    GrB_free (E_arranged + 2) ;                 \
    GrB_free (&pairs_new) ;                     \
    GrB_free (&M_outdeg) ;                      \
    GrB_free (&r_interf) ;                      \
    GrB_free (&r_exists) ;                      \
}

#define LG_FREE_WORK                            \
{                                               \
    /* free any workspace used here */          \
    GrB_free (&E) ;                             \
    GrB_free (&E_half) ;                        \
    GrB_free (&E_t) ;                           \
    GrB_free (&A_tril) ;                        \
    GrB_free (&random_v) ;                      \
    GrB_free (&r_permute) ;                     \
    GrB_free (&r_sorted) ;                      \
    GrB_free (&x) ;                             \
    GrB_free (&r_pairs) ;                       \
    GrB_free (&ramp_v) ;                        \
    GrB_free (&hramp_v) ;                       \
    GrB_free (&swapVals) ;                      \
    GrB_free (&first_bit) ;                     \
    GrB_free (&swap_op) ;                       \
    GrB_free (&any_bxor) ;                      \
    LAGraph_Free((void**)&row_indices, msg) ;   \
    LAGraph_Free((void**)&col_indices, msg) ;   \
    LAGraph_Free((void**)&values, msg) ;        \
    LAGraph_Free((void**)&ramp, msg) ;          \
    LAGraph_Free((void**)&half_ramp, msg) ;     \
    LAGraph_Free((void**)&edge_perm, msg) ;     \
    LAGraph_Free((void**)&swap_type, msg) ;     \
    FREE_LOOP ;                                 \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free(A_new) ;                       \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

void first_bit_equals 
    (uint8_t *z, const uint8_t *x, int64_t i, int64_t j,  const uint8_t *bit)
    {
        (*z) = (uint8_t) (((*x) & 1) == (*bit));
    }
#define FIRST_BIT_EQ                                                            \
"void first_bit_equals                                                       \n"\
"   (uint8_t *z, const uint8_t *x, int64_t i, int64_t j,  const uint8_t *bit)\n"\
"   {                                                                        \n"\
"       (*z) = (uint8_t) (((*x) & 1) == (*bit));                             \n"\
"   }"

// creates [0,3,1,2,1,3,. . .] pattern from random vector.
void swap_pattern 
    (uint8_t *z, const uint64_t *x, int64_t i, int64_t j, const uint8_t *y)
    {
        (*z) = (uint8_t) (((i % 2) * 2) | (*x & 1));
    }
#define SWAP_PAT                                                                \
"void swap_pattern"                                                             \
    "(uint8_t *z, const uint64_t *x, int64_t i, int64_t j, const uint8_t *y)"   \
    "{"                                                                         \
        "(*z) = (uint8_t) (((i % 2) * 2) | (*x & 1));"                          \
    "}"

int LAGraph_SwapEdges
(
    // output
    GrB_Matrix *A_new, //The adjacency matrix of G with edges randomly swapped
    // input: not modified
    LAGraph_Graph G,
    GrB_Index Q, // Swaps per edge
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    GrB_Matrix A = NULL; // n x n Adjacency Matrix 
    GrB_Matrix E = NULL; // e x n Incidence Matrix
    GrB_Matrix E_t = NULL; // n x e Incidence Transposed

    // cmatrix [E_selected, M_1, M_2]
    GrB_Matrix E_arranged[3] = {NULL, NULL, NULL};

    // swaps x e
    // Selected pairs for next batch of swaps
    // Each row contains 2 entries for the edges involved in a swap.
    GrB_Matrix pairs = NULL, pairs_4s = NULL, pairs_i = NULL;
    GrB_Matrix pairs_new = NULL;

    // swaps x n
    // Each row contains 4 or less entries corresponding to the verticies 
    // that are involved in the swap.
    GrB_Matrix M = NULL;
    GrB_Matrix M_fours = NULL; // M with exactly 4 entries
    GrB_Matrix M_i = NULL; // M w/o interference

    // n x swaps
    GrB_Matrix M_t = NULL;

    // swaps x swaps
    // Interference Matrix any entry is the number of verticies the swap on its 
    // row and column share. Any entry with 2 or more could cause interference.
    GrB_Matrix interf = NULL;

    GrB_Matrix M_1 = NULL; // M_1 + M_2 = M
    GrB_Matrix M_2 = NULL; // M_1 + M_2 = M

    // Has a 2 in a certain row if the planned swap already exists in the matrix
    GrB_Matrix exists_1 = NULL, exists_2 = NULL;

    GrB_Index n, e, nvals;

    // Copied vars from LAGraph_Incidence_Matrix
    GrB_Matrix E_half = NULL ;
    GrB_Matrix A_tril = NULL ;

    // random vector 
    GrB_Vector random_v = NULL, r_permute = NULL, r_sorted = NULL;

    GrB_Index *row_indices = NULL ;
    GrB_Index *col_indices = NULL ;
    void *values = NULL ;

    // [0,2,0,3,0,3,0 . . .] numSwaps?
    GrB_Vector swapVals;

    // This ramp will likely get recycled a few times. 
    GrB_Vector ramp_v = NULL;
    GrB_Vector hramp_v = NULL;
    GrB_Index *ramp = NULL ;
    GrB_Index *half_ramp = NULL ;

    // Arrays to unpack edge permutation
    GrB_Index *edge_perm = NULL ;
    uint8_t *swap_type = NULL ;
    bool iso = false;

    // Reduced Vectors
    GrB_Vector M_outdeg = NULL, r_interf = NULL, r_exists = NULL, r_pairs;
    
    // n vector of zeroes
    GrB_Vector x = NULL;

    // Number of values kept in each phase
    GrB_Index n_keep;
    GrB_Index *arr_keep = NULL;
    void *junk = NULL;

    GrB_IndexUnaryOp first_bit = NULL;
    GrB_IndexUnaryOp swap_op = NULL;

    GrB_Semiring any_bxor = NULL;

    A = G->A ;
    //--------------------------------------------------------------------------
    // Check inputs TODO
    //--------------------------------------------------------------------------
    LG_ASSERT_MSG (
        G->kind == LAGraph_ADJACENCY_UNDIRECTED,
        LAGRAPH_INVALID_GRAPH, 
        "G must be undirected"
    ) ;

    LG_ASSERT_MSG (G->nself_edges == 0, LAGRAPH_NO_SELF_EDGES_ALLOWED, 
        "G->nself_edges must be zero") ;
    /* GRB_TRY (GrB_get(A, (void*)&type, GrB_EL_TYPE_CODE)) ;
    LG_ASSERT_MSG (type == GrB_BOOL_CODE, LAGRAPH_INVALID_GRAPH, 
        "A must be type boolean") ; */
    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    LAGRAPH_TRY(LAGraph_Random_Init(msg)) ;
    
    
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_nvals(&nvals, A)) ;
    e = nvals / 2 ;
    GRB_TRY (GrB_Matrix_new(&E, GrB_UINT8, e, n)) ;
    GRB_TRY (GrB_Matrix_new (&E_half, GrB_UINT8, e, n)) ;
    GRB_TRY (GrB_Matrix_new (&E_t, GrB_UINT8, n, e)) ;
    GRB_TRY (GrB_Matrix_new (&A_tril, GrB_BOOL, n, n)) ;
     

    // GRB_TRY (GxB_IndexUnaryOp_new (
    //     &first_bit, (void *) first_bit_equals, GrB_UINT8, GrB_UINT8, GrB_BOOL,
    //     "first_bit", FIRST_BIT_EQ));
    GRB_TRY(GrB_IndexUnaryOp_new(
        &first_bit, (void *) first_bit_equals, GrB_UINT8, GrB_UINT8, GrB_BOOL
    ));
    //GRB_TRY (GxB_IndexUnaryOp_new (
    //    &swap_op, (void *) swap_pattern, GrB_UINT8, GrB_UINT64, GrB_UINT8,
    //    "swap_op", SWAP_PAT));
    GRB_TRY(GrB_IndexUnaryOp_new(
        &swap_op, (void *) swap_pattern, GrB_UINT8, GrB_UINT64, GrB_UINT8
    ));
    GRB_TRY(GrB_Semiring_new(
        &any_bxor, GxB_ANY_UINT8_MONOID ,GrB_BXOR_UINT8 
    ));
    
    GrB_Index num_swaps = 0, num_attempts = 0; 
    // Q: Should this decrease if edges are interfering alot?
    // Bound number of swaps by E-2 * (max deg)?
    // Yes might be good to decrease this if you get alot of intrf
    // or increase if intrf low
    // maybe a * #swaps that worked last iteration
    GrB_Index swaps_per_loop = e / 3 ; // Make this a cap

    

    // This is length e to let every edge have a chance to swap.
    // QUESTION: Perhaps it's better to just do random_v % e and deal with 
    // duplicate edges and self paired edges later. 
    // (inteference matrix will detect)

    GRB_TRY (GrB_Vector_new(&random_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_permute, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&r_sorted, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_UINT8, n)) ;
    GRB_TRY (GrB_Vector_new (&r_pairs, GrB_INT64, e)) ;


    // Extract adjacency matrix to make incidence matrix - E
    // get just the lower triangular entries
    //TODO: should I remove the diagonal? change the 0 if so
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;
    // Arrays to extract A into
    // QUESTION: can't I just use an unpack and use the arrays GBLAS allocated?
    LG_TRY (LAGraph_Malloc ((void**)(&row_indices), e, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&col_indices), e, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&values), e, sizeof(bool), msg)) ;
    GRB_TRY (
        GrB_Matrix_extractTuples_BOOL (row_indices, col_indices, values, &e, A_tril)
        ) ;


    // QUESTION: I don't want to rebuild ramp every time, so do I just build one 
    // that is too large and use as needed?
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&hramp_v, GrB_UINT64, e)) ;
    GRB_TRY (GrB_Vector_new(&swapVals, GrB_UINT8, e)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (ramp_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_UINT64 (ramp_v, NULL, NULL,
        GrB_ROWINDEX_INT64, ramp_v, 0, NULL)) ;
    // [0,0,1,1,2,2,...]
    GRB_TRY (GrB_Vector_apply_BinaryOp2nd_UINT64(
        hramp_v, NULL, NULL, GrB_DIV_UINT64, ramp_v, 2UL, NULL)) ;

    GrB_Index ramp_size;
    GRB_TRY (GxB_Vector_unpack_Full (
        ramp_v, (void **)&ramp, &ramp_size, &iso, NULL)) ;
    GRB_TRY (GxB_Vector_unpack_Full (
        hramp_v, (void **)&half_ramp, &ramp_size, &iso, NULL)) ;

    GrB_Scalar zero8 ;
    GrB_Scalar one8 ;
    GRB_TRY (GrB_Scalar_new (&zero8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_new (&one8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (zero8, 0)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (one8, 1)) ;

    GRB_TRY (GxB_Matrix_build_Scalar (E_half, ramp, col_indices, zero8, e)) ;
    GRB_TRY (GxB_Matrix_build_Scalar (E, ramp, row_indices, one8, e)) ;

    // I'd like to be able to pass in a NULL addition opperator and get an error 
    // if it has to be used
    GRB_TRY (GrB_eWiseAdd (E, NULL, NULL, GrB_PLUS_UINT8, E, E_half, NULL)) ;

    LG_TRY (LAGraph_Free((void**)(&row_indices), msg));
    LG_TRY (LAGraph_Free((void**)(&col_indices), msg));
    LG_TRY (LAGraph_Free((void**)(&values), msg));

    GRB_TRY (GrB_Vector_assign_UINT8 (
        x, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    GRB_TRY (GrB_Vector_assign_UINT64 (
        random_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;

    //TODO: Change seed
    LG_TRY(
        LAGraph_Random_Seed(random_v, 1548945616ul, msg));

    printf("Entering loop, Good Luck:\n") ;
    while(num_swaps < e * Q && num_attempts < e * Q * 5)
    {
        // Coming into the loop: 
        // E must be the incidence matrix of the new graph. W/o self edges nor 
        // parallel edges. Each row must have exactly one 0 and one 1. 
        // They do not have to be randomly assigned.
        // random_v has a radom dense vector.

        GRB_TRY (GrB_Matrix_new(&pairs, GrB_UINT8, swaps_per_loop, e)) ; 
        GRB_TRY (GrB_Matrix_new(&M, GrB_UINT8, swaps_per_loop, n)) ; 
        GRB_TRY (GrB_Vector_new (&M_outdeg, GrB_UINT8, swaps_per_loop)) ;
        
        
        GRB_TRY (GxB_Vector_sort (
            r_sorted, r_permute, GrB_LT_UINT64, random_v, GrB_NULL
        )) ;
        GRB_TRY(GrB_Vector_apply_IndexOp_UINT8(swapVals, NULL, NULL, swap_op, 
            r_sorted, (uint8_t) 0, NULL)) ;
        // NOTE: Small typo on 6.11.6. Refers to vb even though its not a bitmap
        // Also 
        // TODO: handle memory
        // TODO: try and make this more efficient maybe just rand % e rather 
        // than a permutation
        GrB_Index perm_size, swap_size;
        GRB_TRY (GxB_Vector_unpack_Full(
            r_permute, (void **)&edge_perm, &perm_size, &iso, NULL
        )) ;
        uint8_t *swaps = NULL;
        LG_ASSERT(!iso, GrB_NOT_IMPLEMENTED);
        GRB_TRY (GxB_Vector_unpack_Full(
            swapVals, (void **)&swaps, &swap_size, NULL, NULL
        )) ;
        LG_ASSERT(!iso, GrB_NOT_IMPLEMENTED);


        // Pair edges. Take a random permutation and pair adjacent values.
        // pairs wants:
        // cols: [0, 0, 1, 1, 2, 2, 3, 3, . . .] (ramp / 2)
        // rows: [ random permutation 1 - nvals] (LAGraph_Random_Seed + sort)
        // vals: [1, 3, 0, 2, 0, 2, 0, 3, . . .]
            // Start with [0, 2, 0, 2, 0, 2, . . .]
            // for i from 1 to nvals by 2s:
            //      vals[i] += vals[i] > vals[i-1]
        // I know it will be alot less safe but is it more efficient to pack
        // pairs as needed? Doubt it.
        GRB_TRY (GrB_Matrix_build_UINT8(
            pairs, half_ramp, edge_perm, swaps, swaps_per_loop * 2, NULL
        ));

        // QUESTION: do I repack just to throw away?
        GRB_TRY (GxB_Vector_pack_Full(
            r_permute, (void **)&edge_perm, perm_size, false, NULL
        )) ;
        GRB_TRY (GxB_Vector_pack_Full(
            swapVals, (void **)&swaps, swap_size, false, NULL
        )) ;

        // each row of M will have 4 values [0,1,2,3]
        // pairs |  E  |  M
        //   0   | 0,1 | 0,1
        //   2   | 0,1 | 2,3
        //   3   | 0,1 | 3,2
        // This gives us randomization as to which type of swap is happening.
        // M = pairs * E (any_lxor)
        GRB_TRY (GrB_mxm(M, NULL, NULL, any_bxor, pairs, E, NULL)) ;
        
        // TODO: Should there be a function that takes in A matrix and computes
        // in or out degree 
        
        // Remove rows with less than 3 entries
        GRB_TRY (GrB_mxv (M_outdeg, NULL, NULL, LAGraph_plus_one_uint8,
            M, x, NULL)) ; 
        GRB_TRY (GrB_Vector_select_UINT8(
            M_outdeg, NULL, NULL, GrB_VALUEEQ_UINT8, M_outdeg, (uint8_t) 4, NULL)) ;
        GRB_TRY (GrB_Vector_nvals(&n_keep, M_outdeg));

        GrB_Index arr_size, junk_size;
        if(n_keep == swaps_per_loop)
        {
            M_fours = M;
            pairs_4s = pairs;
        }
        else // Remove bad swaps
        {
            GRB_TRY (GxB_Vector_unpack_CSC(
                M_outdeg, &arr_keep, &junk, &arr_size, &junk_size, &iso,
                &n_keep, NULL, NULL)) ;
            LG_ASSERT_MSG(iso, GrB_NOT_IMPLEMENTED, "GraphBLAS unexpectedly gave me"
                " a matrix that is not iso");
            LG_TRY (LAGraph_Free(&junk, msg));
            
            GRB_TRY (GrB_Matrix_new(&M_fours, GrB_UINT8, n_keep, n)) ;
            GRB_TRY (GrB_Matrix_new(&pairs_4s, GrB_UINT8, n_keep, e)) ; 

            GRB_TRY (GrB_Matrix_extract(
                M_fours, NULL, NULL, M, arr_keep, n_keep, GrB_ALL, 0, NULL)) ;
            GRB_TRY (GrB_Matrix_extract(
                pairs_4s, NULL, NULL, pairs, arr_keep, n_keep, GrB_ALL, 0, NULL)) ;
            GRB_TRY (GrB_free(&M));
            GRB_TRY (GrB_free(&pairs)); 
            LG_TRY (LAGraph_Free((void **)&arr_keep, msg));
        }

        // interf <!diagonal> = M*MT(plus_one)
        GRB_TRY (GrB_Matrix_new(&M_t, GrB_UINT8, n, n_keep)) ; 
        GRB_TRY (GrB_transpose(M_t, NULL, NULL, M_fours, NULL)) ;

        GRB_TRY(GrB_Matrix_new(&interf, GrB_UINT8, n_keep, n_keep)) ;
        GRB_TRY(GrB_Vector_new(&r_interf, GrB_UINT8, n_keep)) ;

        //Scuffed way to deal with diagonal. TODO.
        // GRB_TRY(GrB_Vector_assign_UINT8(r_interf, NULL, NULL, (uint8_t) 0, 
        //     GrB_ALL, 0, NULL));
        // GRB_TRY(GrB_Matrix_diag(&interf, r_interf, 0));


        GRB_TRY (GrB_mxm(
            interf, NULL, NULL, LAGraph_plus_one_uint8, M_fours, M_t, NULL)) ;
        // QUESTION: This select is taking forever. 
        // Can I do this with a mask instead??
        GRB_TRY (GrB_Matrix_select_UINT8(
             interf, NULL, NULL, GrB_OFFDIAG, interf, 0, NULL)) ;
        /*     
        GRB_TRY (GrB_Matrix_select_UINT8(
            interf, NULL, NULL, GrB_VALUEGE_UINT8, interf, 2, NULL)) ;
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            r_interf, NULL, NULL, GxB_ANY_UINT8_MONOID, interf, NULL)); 
        */
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            r_interf, NULL, NULL, GrB_MAX_MONOID_UINT8, interf, GrB_DESC_R));

        GRB_TRY (GrB_Vector_select_UINT8(
            r_interf, NULL, NULL, GrB_VALUEGE_UINT8, r_interf, 2, NULL)) ;
        GrB_Index badcount;
        GRB_TRY (GrB_Vector_nvals(&badcount, r_interf)) ;
        if(badcount == 0)
        {
            M_i = M_fours;
            pairs_i = pairs_4s;
        }
        else //delete any interf
        {
            GRB_TRY (GrB_Vector_assign_UINT8(
                r_interf, r_interf, NULL, (uint8_t) 0, GrB_ALL, 0, GrB_DESC_RSC
            )) ;

            GRB_TRY (GxB_Vector_unpack_CSC(
                r_interf, &arr_keep, &junk, &arr_size, &junk_size, &iso,
                &n_keep, NULL, NULL)) ;
            LG_TRY (LAGraph_Free(&junk, msg));

            // Remove bad swaps
            GRB_TRY (GrB_Matrix_new(&M_i, GrB_UINT8, n_keep, n)) ;
            GRB_TRY (GrB_Matrix_new(&pairs_i, GrB_UINT8, n_keep, e)) ; 
            // QUESTION: is it easier to extract on transpose or extract and then remake the transpose?
            GRB_TRY (GrB_Matrix_extract(
                M_i, NULL, NULL, M_fours, arr_keep, n_keep, GrB_ALL, 0, NULL)) ;
            GRB_TRY (GrB_Matrix_extract(
                pairs_i, NULL, NULL, pairs_4s, arr_keep, n_keep, GrB_ALL, 0, NULL
            )) ;
            LG_TRY (LAGraph_Free((void **)&arr_keep, msg)) ;
            GRB_TRY(GrB_free(&M_fours));
        }
        //TODO decide if M_T is still useful

        // M_1 has 2 and 0 M_2 has 1 and 3. Divide by 2 to get 0 and 1 on both.
        // M1 = select M (x & 1 == 0) / 2
        // M2 = select M (x & 1 == 1) / 2
        GRB_TRY (GrB_Matrix_new(&M_1, GrB_UINT8, n_keep, n)) ;
        GRB_TRY (GrB_Matrix_new(&M_2, GrB_UINT8, n_keep, n)) ;
        GRB_TRY (GrB_select(
            M_1, NULL, NULL, first_bit, M_i, 0, NULL)) ;
        GRB_TRY (GrB_select(
            M_2, NULL, NULL, first_bit, M_i, 1, NULL)) ;
        GRB_TRY (GrB_Matrix_apply_BinaryOp2nd_UINT8(
            M_1, NULL, NULL, GrB_DIV_UINT8, M_1, 2, NULL)) ;
        GRB_TRY (GrB_Matrix_apply_BinaryOp2nd_UINT8(
            M_2, NULL, NULL, GrB_DIV_UINT8, M_2, 2, NULL)) ;

        // Check if an edge already exists in the graph.
        // If a row of M1 or M2 has more than 1 vertex in common with a 
        // row of E, the value of that entry in the exists array will be 2. 
        // Otherwise, 1 or noval.
        // exists = M_1 * E (plus_one)
        GRB_TRY (GrB_Matrix_new (&exists_1, GrB_UINT8, n_keep, e)) ; 
        GRB_TRY (GrB_Matrix_new (&exists_2, GrB_UINT8, n_keep, e)) ; 
        GRB_TRY (GrB_Vector_new (&r_exists, GrB_UINT8, n_keep)) ; 
        GRB_TRY (GrB_transpose(E_t, NULL, NULL, E, GrB_DESC_R)) ;
        GRB_TRY (GrB_mxm(
            exists_1, NULL, NULL, LAGraph_plus_one_uint8, M_1, E_t, NULL)) ;
        // exists += M_2 * E (plus_one) (accum is max if theres an intersection)
        GRB_TRY (GrB_mxm(
            exists_2, NULL, NULL, LAGraph_plus_one_uint8, M_2, E_t, NULL)) ;
        //GRB_TRY (GrB_Matrix_select_UINT8(
        //    exists, NULL, NULL, GrB_VALUEEQ_UINT8, exists, (uint8_t) 2, NULL)) ;
        // 2 if intended swap already exists in the matrix, 1 or noval otherwise
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            r_exists, NULL, NULL, GrB_MAX_MONOID_UINT8, exists_1, NULL));
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            r_exists, NULL, GrB_MAX_UINT8, GrB_MAX_MONOID_UINT8, 
            exists_2, NULL));
        GRB_TRY (GrB_Vector_select_UINT8(
            r_exists, NULL, NULL, GrB_VALUEEQ_UINT8, r_exists, (uint8_t) 2, NULL
            )) ;

        GRB_TRY (GrB_Vector_nvals(&badcount, r_exists)) ;
        if(badcount == 0)
        {
            E_arranged[1] = M_1;
            E_arranged[2] = M_2;
            pairs_new = pairs_i;
        }
        else
        {
            // Get compliment vector
            GRB_TRY (GrB_Vector_assign_UINT8(
                r_exists, r_exists, NULL, (uint8_t) 0, GrB_ALL, 0, GrB_DESC_RSC)) ;
            GRB_TRY (GxB_Vector_unpack_CSC(
                r_exists, &arr_keep, &junk, &arr_size, &junk_size, &iso, &n_keep,
                NULL, NULL)) ;
            LG_TRY (LAGraph_Free(&junk, msg)) ;
            
            GRB_TRY (GrB_Matrix_new(E_arranged + 1, GrB_UINT8, n_keep, n)) ;
            GRB_TRY (GrB_Matrix_new(E_arranged + 2, GrB_UINT8, n_keep, n)) ;
            GRB_TRY (GrB_Matrix_new(&pairs_new, GrB_UINT8, n_keep, e)) ;

            // TODO: Check if this needs to be transposed I think it does
            GRB_TRY (GrB_Matrix_extract(
                E_arranged[1], NULL, NULL, M_1, arr_keep, n_keep, 
                GrB_ALL, 0, NULL)) ;
            GRB_TRY (GrB_Matrix_extract(
                E_arranged[2], NULL, NULL, M_2, arr_keep, n_keep, 
                GrB_ALL, 0, NULL)) ;
            GRB_TRY (GrB_Matrix_extract(
                pairs_new, NULL, NULL, pairs_i, arr_keep, n_keep, 
                GrB_ALL, 0, NULL)) ;
            LG_TRY (LAGraph_Free((void **)&arr_keep, msg));
        }

        GrB_Index n_old;
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            r_pairs, NULL, NULL, GxB_ANY_UINT8_MONOID, pairs_new, GrB_DESC_RT0));

        GRB_TRY (GrB_Vector_assign_UINT8(
            r_pairs, r_pairs, NULL, (uint8_t) 0, GrB_ALL, 0, GrB_DESC_RSC)) ;
        
        GRB_TRY (GxB_Vector_unpack_CSC(
            r_pairs, &arr_keep, &junk, &arr_size, &junk_size, &iso, &n_old,
            NULL, NULL)) ;
        LG_TRY (LAGraph_Free(&junk, msg)) ;
        
        LG_ASSERT(n_old ==  e - 2 * n_keep, GrB_NOT_IMPLEMENTED) ;
        GRB_TRY (GrB_Matrix_new(E_arranged, GrB_UINT8, n_old, n)) ;
        GRB_TRY (GrB_Matrix_extract(E_arranged[0], NULL, NULL, E, 
            arr_keep, n_old, GrB_ALL, n, NULL)) ;
        LG_TRY (LAGraph_Free((void **)&arr_keep, msg));

        
        // E = Concat(E_prime, M_1, M_2)
        // where E_prime contains no edges in the indegree of pairs after it is 
        // pruned
        GRB_TRY (GxB_Matrix_concat(E, E_arranged, 3, 1, NULL)) ;
        
        // Free Matricies that have to be rebuilt
        FREE_LOOP ;

        num_attempts += swaps_per_loop;
        num_swaps += n_keep;
        swaps_per_loop = (n_keep * 3) / 2;
        swaps_per_loop = LAGRAPH_MAX(swaps_per_loop, 16) ;
        swaps_per_loop = LAGRAPH_MIN(swaps_per_loop, e / 3) ;
        // Maintain random Vector
        LG_TRY (LAGraph_Random_Next(random_v, msg)) ;
        printf("#####Made %ld swaps. Total %ld out of %ld. Attempting %ld swaps next.#####\n\n", n_keep, num_swaps, e * Q, swaps_per_loop);
    } 

    GRB_TRY(GrB_Matrix_new(A_new, GrB_UINT8, n, n)) ;
    GRB_TRY (GrB_mxm(
        *A_new, NULL, NULL, LAGraph_any_one_uint8, E, E, GrB_DESC_T0)) ;
    GRB_TRY(GrB_Matrix_select_UINT8(
        *A_new, NULL, NULL, GrB_OFFDIAG, *A_new, (uint8_t) 0, NULL)) ;

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
