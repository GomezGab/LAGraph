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


#define LG_FREE_WORK                        \
{                                           \
    /* free any workspace used here */      \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_HelloWorld // a simple algorithm, just for illustration
(
    // output
    GrB_Matrix *A_new,
    // input: not modified
    LAGraph_Graph G,
    GrB_Index Q, // Swaps per edge
    char *msg
)
{
    return (GrB_NOT_IMPLEMENTED) ;
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    GrB_Matrix A = NULL; // nxn Adjacency Matrix 
    GrB_Matrix E = NULL; // nxe Incidence Matrix
    //GrB_Matrix E_t = NULL; // Incidence Transposed
    GrB_Matrix E_selected = NULL; // E - rows to be replaced

    // Selected pairs for next batch of swaps
    // Each row? contains 2 entries for the edges involved in a swap.
    GrB_Matrix pairs = NULL;

    // Each row? contains 4 or less entries corresponding to the verticies 
    // that are involved in the swap.
    GrB_Matrix M = NULL;
    GrB_Matrix M_t = NULL;
    GrB_Matrix M_outdiag = NULL;

    // Interference Matrix any entry is the number of verticies the swap on its 
    // row and column share. Any entry with 2 or more could cause interference.
    GrB_Matrix interf = NULL;

    GrB_Matrix M_1 = NULL; // M_1 + M_2 = M
    GrB_Matrix M_2 = NULL; // M_1 + M_2 = M

    // Diagonalized version of the reduced interference matrix
    GrB_Matrix d_interf = NULL;

    // Has a 2 in a certain row if the planned swap already exists in the matrix
    GrB_Matrix exists = NULL;

    // New edges after swap, positioned to replace edges in E 
    GrB_Matrix S = NULL; 

    GrB_Index n, e, nvals;

    // Copied vars from LAGraph_Incidence_Matrix
    GrB_Matrix E_half = NULL ;
    GrB_Matrix A_tril = NULL ;

    // random vector 
    GrB_Vector random_v = NULL, r_permute = NULL, r_sorted = NULL;

    GrB_Index *row_indices = NULL ;
    GrB_Index *col_indices = NULL ;
    void *values = NULL ;

    // This ramp will likely get recycled a few times. 
    GrB_Vector ramp_v = NULL;
    GrB_Vector hramp_v = NULL;
    GrB_Index *ramp = NULL ;
    GrB_Index *half_ramp = NULL ;

    // Arrays to unpack edge permutation
    GrB_Index *edge_perm = NULL ;
    u_char *swap_type = NULL ;
    bool iso = false;

    // Reduced Vectors
    GrB_Vector M_outdeg = NULL, r_interf = NULL, r_exists = NULL;
    
    // x = zeros (swaps_per_loop,1)
    GrB_Vector x = NULL;

    // Number of values removed in each phase
    GrB_Index n_4degree, n_interf, n_exist;
    
    //--------------------------------------------------------------------------
    // Check inputs TODO
    //--------------------------------------------------------------------------
    

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------

    A = G->A ;

    char typename[LAGRAPH_MAX_NAME_LEN] ;
    GrB_Type type ;
    LG_TRY (LAGraph_Matrix_TypeName (typename, A, msg)) ;
    LG_TRY (LAGraph_TypeFromName (&type, typename, msg)) ;

    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Matrix_nvals(&nvals, A)) ;
    e = nvals / 2 ;
    GRB_TRY (GrB_Matrix_new(&E, GrB_UINT8, n, e)) ;
    GRB_TRY (GrB_Matrix_new (&E_half, GrB_UINT8, n, e)) ;
    GRB_TRY (GrB_Matrix_new (&A_tril, type, n, n)) ;
    GRB_TRY (GrB_Matrix_new(&E_selected, GrB_UINT8, n, e)) ; 

    GrB_Index num_swaps = 0, num_attempts = 0 ; 
    // QUESTION: Should this decrease if edges don't want to swap?
    GrB_Index swaps_per_loop = e / 3 ; // 1/3 reasonable?

    GRB_TRY (GrB_Matrix_new(&pairs, GrB_UINT8, e, swaps_per_loop)) ; 
    GRB_TRY (GrB_Matrix_new(&M, GrB_UINT8, n, swaps_per_loop)) ; 


    GRB_TRY (GrB_Vector_new(&random_v, GrB_UINT64, swaps_per_loop * 2)) ;
    GRB_TRY (GrB_Vector_new(&r_permute, GrB_UINT64, swaps_per_loop * 2)) ;
    GRB_TRY (GrB_Vector_new(&r_sorted, GrB_UINT64, swaps_per_loop * 2)) ;
    GRB_TRY (GrB_Vector_new(&ramp_v, GrB_UINT64, swaps_per_loop * 2)) ;
    GRB_TRY (GrB_Vector_new(&hramp_v, GrB_UINT64, swaps_per_loop * 2)) ;
    GRB_TRY (GrB_Vector_new (&M_outdeg, GrB_INT64, n)) ;
    GRB_TRY (GrB_Vector_new (&x, GrB_INT64, swaps_per_loop)) ;
    
    // Extract adjacency matrix to make incidence matrix - E
    // get just the lower triangular entries
    GRB_TRY (GrB_select (A_tril, NULL, NULL, GrB_TRIL, A, 0, NULL)) ;
    // Arrays to extract A into
    LG_TRY (LAGraph_Malloc ((void**)(&row_indices), e, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&col_indices), e, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&values), e, sizeof(type), msg)) ;
    GRB_TRY (
        GrB_Matrix_extractTuples (row_indices, col_indices, values, &e, A_tril)
        ) ;
    
    /* LG_TRY (LAGraph_Malloc ((void**)(&ramp), e, sizeof(GrB_Index), msg)) ;
    LG_TRY (LAGraph_Malloc ((void**)(&half_ramp), e, sizeof(GrB_Index), msg)) ;


    #pragma omp parallel for
    for (GrB_Index i = 0 ; i < e ; i++) 
    {
        ramp[i] = i ;
    }

    #pragma omp parallel for
    for (GrB_Index i = 0 ; i < swaps_per_loop * 2 ; i++) 
    {
        half_ramp[i] = i / 2;
    } */

    GrB_Scalar zero8 ;
    GrB_Scalar one8 ;
    GRB_TRY (GrB_Scalar_new (&zero8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_new (&one8, GrB_UINT8)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (zero8, 0)) ;
    GRB_TRY (GrB_Scalar_setElement_UINT8 (one8, 1)) ;
    GRB_TRY (GxB_Matrix_build_Scalar (E_half, col_indices, ramp, zero8, e)) ;
    GRB_TRY (GxB_Matrix_build_Scalar (E, row_indices, ramp, one8, e)) ;

    // I'd like to be able to pass in a NULL addition opperator and get an error 
    // if it has to be used
    GRB_TRY (GrB_eWiseAdd (E, NULL, NULL, GrB_PLUS_INT8, E, E_half, NULL)) ;

    // Initialize random vector
    /* GrB_Scalar zero64;
    GRB_TRY (GrB_Scalar_new (&zero64, GrB_UINT8)) ;

    //value irrelevant
    GRB_TRY (GrB_Scalar_setElement_UINT64 (zero64, 0)) ;
    // QUESTION: Should this be GxB_Vector_pack??
    GRB_TRY (
        GxB_Vector_build_Scalar(random_v, ramp, zero64, e)); */
    // TODO check this
    uint8_t zero_8 = 0;
    
    // QUESTION: diffrence? which is better?
    /* GRB_TRY(GxB_Vector_pack_Full(
        random_v, &zero_8, sizeof(uint8_t), true, NULL
    )) ; */
    GRB_TRY (GrB_assign (x, NULL, NULL, 0, GrB_ALL, swaps_per_loop, NULL)) ;
    GRB_TRY (GrB_assign (random_v, NULL, NULL, 0, GrB_ALL, e, NULL)) ;
    GRB_TRY (GrB_Vector_apply_IndexOp_INT64 (ramp_v, NULL, NULL,
        GrB_ROWINDEX_INT64, random_v, 1, NULL)) ;
    GRB_TRY (GrB_Vector_apply (half_ramp, NULL, NULL, half, ramp_v, NULL)) ;
    GRB_TRY (GxB_Vector_unpack_Full (
        ramp_v, ramp, e * sizeof(GrB_Index), NULL, NULL)) ;
    GRB_TRY (GxB_Vector_unpack_Full (
        hramp_v, half_ramp, e * sizeof(GrB_Index), NULL, NULL)) ;

    //TODO Change seed
    LG_TRY(
        LAGraph_Random_Seed(random_v, 15489451345495616ul, msg));

    while(num_swaps < e * Q && num_attempts < e * Q * 5)
    {
        // Coming into the loop: 
        // E must be the incidence matrix of the new graph. W/o self edges nor 
        // parallel edges. Each row must have exactly one 0 and one 1. 
        // They do not have to be randomly assigned.
        // random_v has a radom dense vector.
        GRB_TRY (GxB_Vector_sort (
            r_sorted, r_permute, GrB_LT_UINT64, random_v, GrB_NULL
        )) ;
        //GrB_Vector_apply(r_sorted, NULL, NULL, GrB_UnaryOp, r_sorted, GrB_DESC_R) ;
        // NOTE: Small typo on 6.11.6. Refers to vb even though its not a bitmap
        // Also 
        GRB_TRY (GxB_Vector_unpack_Full(
            r_permute, edge_perm, e * sizeof(GrB_Index), iso, GrB_NULL
        )) ;


        // Pair edges. Take a random permutation and pair adjacent values.
        // pairs wants:
        // rows: [0, 0, 1, 1, 2, 2, 3, 3, . . .] (ramp / 2)
        // cols: [ random permutation 1 - nvals] (LAGraph_Random_Seed + sort)
        // vals: [0, 3, 0, 2, 0, 2, 0, 3, . . .]
            // Start with [0, 2, 0, 2, 0, 2, . . .]
            // for i from 1 to nvals by 2s:
            //      vals[i] += vals[i] > vals[i-1]
        // I know it will be alot less safe but is it more efficient to pack
        // pairs as needed? Doubt it.
        GRB_TRY (GrB_Matrix_build_UINT8(
            pairs, edge_perm, half_ramp, vals, swaps_per_loop, GrB_NULL
        )) ;


        // each row of M will have 4 values [0,1,2,3]
        // pairs |  E  |  M
        //   0   | 0,1 | 0,1
        //   2   | 0,1 | 2,3
        //   3   | 0,1 | 3,2
        // This gives us randomization as to which type of swap is happening.
        // M = E * pairs (any_bxor)
        GRB_TRY (GrB_mxm(M, NULL, NULL, GxB_ANY_LXOR_INT8, E, pairs, NULL)) ;
        // remove any entries of M with less than 3 values
        // This is taken from LAGraph_Cached_OutDegree. 
        
        
        GRB_TRY (GrB_mxv (M_outdeg, NULL, NULL, LAGraph_plus_one_int64,
            M, x, NULL)) ; 

        // QUESTION: Is the following method correct? Do I need a new M with the 
        // correct dimension?
        GRB_TRY (GrB_Vector_select_UINT8(
            M_outdeg, NULL, NULL, GrB_EQ_INT8, M_outdeg, 4, NULL)) ;
        GRB_TRY (GrB_Vector_nvals(&n_4degree,M_outdeg)) ;
        GrB_Index *deg_keep = NULL;
        void *junk = NULL;
        GrB_Index junk_size;
        GRB_TRY (GrB_Vector_unpack_CSC(
            M_outdeg, &deg_keep, &junk, &n_4degree, &junk_size, NULL)) ;
        GRB_TRY (GrB_Matrix_extract(
            M, NULL, NULL, M, NULL, 0, deg_keep, n_4degree, NULL)) ;


        // interf <!diagonal> = M*MT(plus_one)
        // Should I keep rebuilding matricies with the correct size??
        GRB_TRY (GrB_transpose(M_t, NULL, NULL, M, NULL)) ;
        GRB_TRY (GrB_mxm(
            interf, NULL, NULL, LAGraph_plus_one_uint8, M_t, M, GrB_DESC_R)) ;
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            r_interf, NULL, NULL, GrB_MAX_MONOID_UINT8, interf, GrB_DESC_R));
        GRB_TRY (GrB_Vector_select_UINT8(
            r_interf, NULL, NULL, GrB_EQ_INT8, r_interf, 2, NULL)) ;
        // TODO Remove anything in r_interf with best method

        // set 3rd bit to 1 if it must be removed. 
        // An accumulator might work better.
        // M = d_interf * M (any_[z = 4 * (x > 2) | y]) 
        // M = select M (x & 4 == 0)

        // M1 = select M (x & 1 == 0) / 2
        GRB_TRY (GrB_Matrix_select_UINT8(
            M_1, NULL, NULL, first_bit, M_t, 0, GrB_DESC_R)) ;
        // M2 = select M (x & 1 == 1) / 2
        GRB_TRY (GrB_Matrix_select_UINT8(
            M_2, NULL, NULL, first_bit, M_t, 1, GrB_DESC_R)) ;

        // Check if an edge already exists in the graph.
        // If a row of M1 or M2 has more than 1 vertex in common with a 
        // row of E, the value of that entry in the exists array will be 2. 
        // Otherwise, 1 or noval.
        // exists = M_1 * E (plus_one)
        GRB_TRY (GrB_mxm(
            exists, NULL, NULL, LAGraph_plus_one_uint8, M_1, E, GrB_DESC_R)) ;
        // exists += M_2 * E (plus_one) (accum is one if theres an intersection)
        GRB_TRY (GrB_mxm(
            exists, NULL, GrB_MAX_MONOID_UINT8, 
            LAGraph_plus_one_uint8, M_2, E, NULL)) ;
        // 2 if intended swap already exists in the matrix, 1 or noval otherwise
        GRB_TRY (GrB_Matrix_reduce_Monoid(
            r_exists, NULL, NULL, GrB_MAX_MONOID_UINT8, exists, GrB_DESC_R));

        // Compute S. note we could do ST if its better.
        // pairs |  M  |  S
        //   0   | 0,1 | 0,1
        // 2 / 3 | 0,1 | 2,3
        //   0   | 2,3 | 2,3
        // 2 / 3 | 2,3 | 0,1
        // S = pairsT * M (any_[z = (x & 2) ^ y])
        // S = extract S (x & 2 == 0)

        // If S has a row, Remove it from E and replace with the row from S
        // QUESTION

        // update E and possibly ET
    } 

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
