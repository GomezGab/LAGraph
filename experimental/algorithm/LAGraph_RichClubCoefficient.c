//------------------------------------------------------------------------------
// LAGraph_RichClubCoefficient: rich club coefficient
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

// Get the rich club coefficient of a graph, also allows for edge randomization
// to normalize the coefficients in a graph.

// Given a Symetric Graph with no self edges, LAGraph_RichClubCoefficient will
// first randomized edges without changing the degree pattern and will then 
// calculate the rich club coefficients of the resulting graph.

// The values will be output as a dense GrB_Vector, the rich club coefficient 
// of the kth degree found at entry k.

// References:

// Julian J. McAuley, Luciano da Fontoura Costa, and Tibério S. Caetano, “The 
// rich-club phenomenon across complex network hierarchies”, Applied Physics 
// Letters Vol 91 Issue 8, August 2007. https://arxiv.org/abs/physics/0701290

// R. Milo, N. Kashtan, S. Itzkovitz, M. E. J. Newman, U. Alon, “Uniform 
// generation of random graphs with arbitrary degree sequences”, 2006. 
// https://arxiv.org/abs/cond-mat/0312028

#define LG_FREE_WORK                                    \
{                                                       \
    /* free any workspace used here */                  \
    GrB_free (&edge_degrees) ;                          \
    GrB_free (&D) ;                                     \
    GrB_free (&degrees) ;                               \
    GrB_free (&node_edges) ;                            \
    GrB_free (&edges_per_deg) ;                         \
    GrB_free (&rcCalculation) ;                         \
    LAGraph_Free((void **) &index_edge, NULL) ;         \
    LAGraph_Free((void **) &node_edges_arr, NULL); \
    LAGraph_Free((void **) &deg_arr, NULL);             \
    LAGraph_Free((void **) &edges_per_deg_arr, NULL);   \
    LAGraph_Free((void **) &cumul_array, NULL);         \
    LAGraph_Free((void **) &ones, NULL);                \
    LAGraph_Free((void **) &deg_vertex_count, NULL);    \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (rich_club_coefficents) ;      \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

//taken from LAGraph_BF_full_mxv
typedef void (*LAGraph_binary_function) (void *, const void *, const void *) ;

/* #define TWO_ONE_ADD                                                         \
    "void two_one_add(uint64_t *z, const uint64_t *x, const uint64_t *y)"   \
    "{"                                                                     \
        "(*z) = 2 * (*x) + (*y) ;"                                          \
    "}"
void two_one_add(uint64_t *z, const uint64_t *x, const uint64_t *y)
{ 
    (*z) = 2 * (*x) + (*y);
}
 */
#define ISEQ_2ISLT                                                          \
    "void iseq_2islt(uint64_t *z, const uint64_t *x, const uint64_t *y)"    \
    "{"                                                                     \
        "(*z) = (*x < *y) + (*x <= *y) ;"                                   \
    "}"
void iseq_2islt(uint64_t *z, const uint64_t *x, const uint64_t *y)
{
    (*z) = (uint64_t)((*x < *y) + (*x <= *y)) ;
}

#define RICH_CLUB_FORMULA                                                      \
    "void rich_club_formula(double *z, const uint64_t *x, const uint64_t *y)"  \
    "{"                                                                        \
        "(*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0)) ;"  \
    "}"
//TODO: Look out for ones
void rich_club_formula(double *z, const uint64_t *x, const uint64_t *y)
{
    (*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0));
} 

int LAGraph_RichClubCoefficient
(
    // output:
    //rich_club_coefficents(i): rich club coefficents of i
    GrB_Vector *rich_club_coefficents,    

    // input: 
    LAGraph_Graph G, //input graph
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    LG_CLEAR_MSG ;
    //Matrix containing every edge 
    //With an entry cooresponding to the degree of its column
    GrB_Matrix edge_degrees = NULL;

    //A matrix with diagonal entries corresponding to degrees.
    GrB_Matrix D = NULL;

    //degrees of nodes.
    GrB_Vector degrees = NULL;

    // contains the number of edges for which the ith node is
    // the smallest degree node * 2 + # edges w/ same degree as the other node
    // to account for double counting of edges w/ same degree as the other node.
    GrB_Vector node_edges = NULL;

    // the ith entry contains the number of edges whose lowest degree is i.
    GrB_Vector edges_per_deg = NULL;

    // the ith entry contains the number of verticies whose degree is i.
    GrB_Vector verts_per_deg = NULL;

    // 2 * (x < y) + (x == y)
    GrB_BinaryOp iseq_2lt = NULL;

    // [+].[iseq_2lt]
    GrB_Semiring plus_2le;

    // 2E_K / (N_k (N_k -1))
    GrB_BinaryOp rcCalculation = NULL;

    GrB_Matrix A ; // G->A, the adjacency matrix
    GrB_Index n ;
    GrB_Index vi_size ;
    GrB_Index vx_size ;
    GrB_Index edge_vec_nvals;
    GrB_Index deg_vec_size;
    bool iso = false;

    GrB_Index *index_edge = NULL;
    
    uint64_t *node_edges_arr = NULL, *deg_arr = NULL, 
        *edges_per_deg_arr = NULL, *cumul_array = NULL, *ones = NULL, 
        *deg_vertex_count = NULL;

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rich_club_coefficents != NULL, GrB_NULL_POINTER);

    //TODO: double check this
    LG_ASSERT_MSG(
        G->kind == LAGraph_ADJACENCY_UNDIRECTED, GrB_INVALID_VALUE, 
        "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->out_degree != NULL, GrB_EMPTY_OBJECT,
        "G->out_degree must be defined") ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE, 
        "G->nself_edges must be zero") ; 

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;
    GRB_TRY(GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY(GrB_Matrix_new(&D, GrB_UINT64, n, n))
    GRB_TRY(GrB_Matrix_new(&edge_degrees, GrB_UINT64,n,n)) ;


    GRB_TRY(GrB_Vector_new(&degrees, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&node_edges, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edges_per_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&verts_per_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(rich_club_coefficents, GrB_FP64, n)) ;

    GRB_TRY(GxB_BinaryOp_new(
        &iseq_2lt, (LAGraph_binary_function) (&iseq_2islt), 
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "iseq_2islt", ISEQ_2ISLT)) ;
    GRB_TRY(GrB_Semiring_new(&plus_2le, GrB_PLUS_MONOID_UINT64, iseq_2lt)) ;
    GRB_TRY(GxB_BinaryOp_new(
        &rcCalculation, (LAGraph_binary_function) (&rich_club_formula), 
        GrB_FP64, GrB_UINT64, GrB_UINT64, 
        "rich_club_formula", RICH_CLUB_FORMULA)) ;

    // degrees = G->out_degree
    GRB_TRY (GrB_assign(
        degrees, NULL, NULL, G->out_degree, GrB_ALL, n, NULL)) ;

    GRB_TRY (GrB_Matrix_diag(&D, degrees, 0)) ;

    //--------------------------------------------------------------------------
    // Calculating time
    //--------------------------------------------------------------------------

    // QUESTION: does this Mask help out GBLAS?
    // Each edge in the graph gets the value of the degree of its column node
    GRB_TRY (GrB_mxm(
        edge_degrees, A, NULL, GxB_ANY_FIRST_UINT64,D,A, GrB_DESC_S)) ;

    // QUESTION: Is it not more fficient to simply use min here and then count up?
    GRB_TRY(GrB_mxm(
        edge_degrees, NULL, NULL, plus_2le, edge_degrees, D, GrB_NULL)) ;

    // If the nodes of an edge have different degrees, the edge is counted once.
    // If they have the same degree, that edge is double counted. So, we adjust:
    GRB_TRY(GrB_Matrix_reduce_Monoid(
        node_edges, NULL, NULL, GrB_PLUS_MONOID_UINT64, edge_degrees, NULL)) ;

    // The rest of this is indexing the number of edges and number of nodes at 
    // each degree and then doing a cummulative sum to know the amount of edges 
    // and nodes at degree geq k.
    GRB_TRY(GrB_Vector_nvals (&edge_vec_nvals, node_edges));
    vi_size = (edge_vec_nvals+1)*sizeof(GrB_Index);
    vx_size = (edge_vec_nvals+1)*sizeof(GrB_UINT64);

    // Grab the index and edge count arrays from GBLAS
    // Jubled NULL so must return sorted. Needed because arrays with 
    // # of edges and # of degrees should line up.

    GRB_TRY(GxB_Vector_unpack_CSC(
        node_edges, &index_edge, (void **) &node_edges_arr,
        &vi_size,&vx_size,&iso,&edge_vec_nvals,NULL, GrB_NULL)) ;
    LG_TRY(LAGraph_Free((void **)&index_edge, NULL)) ;

    //TODO: adjust degrees over by one
    GRB_TRY(GxB_Vector_unpack_CSC(
        degrees, &index_edge, (void **) &deg_arr,
        &vi_size,&vx_size,&iso,&deg_vec_size,NULL, GrB_NULL)) ;

    LG_TRY(LAGraph_Free((void **)&index_edge, NULL)) ;

    // TODO change what this throws
    LG_ASSERT (edge_vec_nvals == deg_vec_size, GrB_NULL_POINTER) ;

    // Build with degrees as indecies and handle duplicates via adition
    GRB_TRY(GrB_Vector_build (
        edges_per_deg, deg_arr, node_edges_arr, deg_vec_size, 
        GrB_PLUS_UINT64)) ;
    // TODO: Make ones array in a better way
    LG_TRY(
        LAGraph_Malloc((void **) &ones, deg_vec_size, sizeof(u_int64_t), NULL)) ;
    for(uint64_t i = 0; i < deg_vec_size; ++i)
        ones[i] = 1;

    GRB_TRY(GrB_Vector_build (
        verts_per_deg, deg_arr, ones, deg_vec_size, GrB_PLUS_UINT64)) ;

    GrB_Index *epd_index = NULL,  *vpd_index = NULL;
    GRB_TRY(GxB_Vector_unpack_CSC(
        edges_per_deg, &epd_index, (void **)&edges_per_deg_arr,
        &vi_size, &vx_size, &iso, &edge_vec_nvals, NULL, GrB_NULL)) ;
    GRB_TRY(GxB_Vector_unpack_CSC(
        verts_per_deg, &vpd_index, (void **)&deg_vertex_count,
        &vi_size, &vx_size, &iso, &deg_vec_size, NULL, GrB_NULL)) ;

    //run a cummulative sum (backwards) on deg_vertex_count
    for(uint64_t i = deg_vec_size - 1; i > 0; --i)
    {
        deg_vertex_count[i-1]+=deg_vertex_count[i];
    }

    //run a cummulative sum (backwards) on edges_per_deg_arr
    for(uint64_t i = edge_vec_nvals - 1; i > 0; --i)
    {
        edges_per_deg_arr[i-1]+=edges_per_deg_arr[i];
    }

    // QUESTION: can I just tell GBLAS the arrays are one smaller than they are?
    // prevent division by zero by removing a possible one from the 
    // verts_per_deg
    GrB_Index to_remove = deg_vertex_count[deg_vec_size - 1] == 1? 
        vpd_index[deg_vec_size - 1]: 0;
    printf("Last vertex count is %ld", deg_vertex_count[deg_vec_size - 1]);
    //re pack but now we're cummulative
    GRB_TRY(GxB_Vector_pack_CSC(
        edges_per_deg, &epd_index, (void **)&edges_per_deg_arr,
        vi_size, vx_size, false, edge_vec_nvals, NULL, GrB_NULL));
    GRB_TRY(GxB_Vector_pack_CSC(
        verts_per_deg, &vpd_index, (void **)&deg_vertex_count,
        vi_size, vx_size, false, deg_vec_size, NULL, GrB_NULL));
    //GxB_Vector_fprint (verts_per_deg, "vpd", GxB_SHORT, stdout);
    if(to_remove > 0) 
        GRB_TRY(GrB_Vector_removeElement(verts_per_deg, to_remove));
    //GxB_Vector_fprint (verts_per_deg, "vpd", GxB_SHORT, stdout);
    GRB_TRY(GrB_eWiseMult(*rich_club_coefficents, NULL, NULL, rcCalculation, 
        edges_per_deg, verts_per_deg, NULL)) ;
    //GxB_Vector_fprint (*rich_club_coefficents, "rcc", GxB_SHORT, stdout);
    /* The following is not nessesary because we can redefine the output vector 
    * as sparse vect
    //Construct a vector thats has edges_per_deg_arr values but repeated  
    // whenever index_edge has a skip and put into cumulative edges
    // ie. [0,1,6,7,10] & [9,7,3,2,0] ->
    // [9,7,7,7,7,7,3,2,2,2,. . . ]

    uint64_t index = 0, i = 0;
    LG_TRY (LAGraph_Malloc((void **) &cumul_array, n, sizeof(uint64_t), msg)) ;

    for(; i < n; i++) // seems easily parrallelizable but idk if #pragma is enough.
    {
        if(index + 1 <edge_vec_nvals && i == index_edge[index+1])
            ++ index;
        cumul_array[i] = edges_per_deg_arr[index];
    } 
    */
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}