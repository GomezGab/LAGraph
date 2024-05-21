//------------------------------------------------------------------------------
// LAGraph_RichClubCoefficient: rich club coefficients and edge randomization
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

#define LG_FREE_WORK                        \
{                                           \
    /* free any workspace used here */      \
    GrB_free (&edge_degrees) ;              \
    GrB_free (&D) ;                         \
    GrB_free (&degrees) ;                   \
    GrB_free (&edge_gt_deg) ;               \
    GrB_free (&edge_eq_deg) ;               \
    GrB_free (&edge_adjusted_deg) ;         \
    GrB_free (&cumulative_deg) ;            \
    GrB_free (&edges_per_deg) ;             \
    GrB_free (&cumulative_edges) ;          \
    GrB_free (&two_one) ;                   \
    GrB_free (&rcCalculation) ;             \
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

#define TWO_ONE_ADD                                                         \
    "void two_one_add(uint64_t *z, const uint64_t *x, const uint64_t *y)"   \
    "{"                                                                     \
        "(*z) = 2 * (*x) + (*y) ;"                                          \
    "}"
void two_one_add(uint64_t *z, const uint64_t *x, const uint64_t *y)
{ 
    (*z) = 2 * (*x) + (*y);
}

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
//Look out for ones
void rich_club_formula(double *z, const uint64_t *x, const uint64_t *y)
{
    (*z) = ((double)(*x)) / (((double)(*y)) * (((double)(*y)) - 1.0)) ;
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
    return (GrB_NOT_IMPLEMENTED) ;
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

    //contains the number of edges for which the ith node is
    //the smallest degree node in the pair
    GrB_Vector edge_gt_deg = NULL;

    //contains the number of edges for which the ith node is
    //the same degree as the other node in the pair
    //used to correct for undercounted nodes
    GrB_Vector edge_eq_deg = NULL;

    //Combines edge_gt_deg and edge_eq_deg to 
    // account for double counting in edge_eq_deg
    GrB_Vector edge_adjusted_deg = NULL;

    //the ith entry contains the number of nodes with degree greater than i.
    GrB_Vector cumulative_deg = NULL;

    // the ith entry contains the number of edges whos lowest degree is i.
    GrB_Vector edges_per_deg = NULL;

    //the ith entry contains twice the number of edges among 
    // nodes with degree greater than i.
    GrB_Vector cumulative_edges = NULL;

    // 2x+y
    GrB_BinaryOp two_one = NULL;

    // 2 * (x > y) + (x == y)
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

    //TODO: YOU STILL HAVE TO GRBFREE THIS WITH A LOOP?
    //Should the arrays be initialized before fed to the function?
    GrB_Index *index_edge = NULL; // does this need to be null to start with? 
    
    uint64_t *edge_count_per_node = NULL, *deg_arr = NULL, 
        *edges_per_deg_arr = NULL, *cumul_array = NULL, *ones = NULL, 
        *ramp = NULL;

    //TODO: YOU STILL HAVE TO GRBFREE THIS WITH A LOOP? 
    GrB_Index **index_degree = NULL; // does this need to be null to start with?
    void **degree_array = NULL;



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
        "G->nself_edges must be zero") ; // Probably overkill

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;
    GRB_TRY(GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY(GrB_Matrix_new(&D, GrB_UINT64, n, n))
    GRB_TRY(GrB_Matrix_new(&edge_degrees, GrB_UINT64,n,n)) ;
    GRB_TRY(GrB_Vector_new(&degrees, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edge_gt_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edge_eq_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&cumulative_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edges_per_deg, GrB_UINT64, n)) ;
    
    GRB_TRY(GxB_BinaryOp_new(
        &two_one, (LAGraph_binary_function) (&two_one_add), 
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "two_one_add", TWO_ONE_ADD)) ;
    GRB_TRY(GxB_BinaryOp_new(
        &iseq_2lt, (LAGraph_binary_function) (&iseq_2islt), 
        GrB_UINT64, GrB_UINT64, GrB_UINT64, "iseq_2islt", ISEQ_2ISLT)) ;
    GRB_TRY(GrB_Semiring_new(&plus_2le, GrB_PLUS_MONOID_UINT64, iseq_2lt)) ;
    GRB_TRY(GxB_BinaryOp_new(
        &rcCalculation, (LAGraph_binary_function) (&rich_club_formula), 
        GrB_UINT64, GrB_UINT64, GrB_UINT64, 
        "rich_club_formula", RICH_CLUB_FORMULA)) ;
    // degrees = G->out_degree
    GRB_TRY (GrB_assign(
        degrees, NULL, NULL, G->out_degree, GrB_ALL, n, NULL)) ;

    GRB_TRY (GrB_Matrix_diag(&D, degrees, 0)) ;

    // Each edge in the graph gets the value of the degree of its column node
    GRB_TRY (GrB_mxm(
        edge_degrees, NULL, NULL, GxB_ANY_SECOND_UINT64,A,D, GrB_NULL)) ;
    /* 
    // QUESTION: Mentioned GxB would be slower than JIT correct?

    // Counts the edges for which its row node is the node of smallest degree.
    // Simply, this edge would be removed from the graph once this node
    // is of degree < k. 
    GRB_TRY(GrB_vxm(
        edge_gt_deg, NULL, NULL, GxB_PLUS_ISLT_UINT64, 
        degrees, edge_degrees, GrB_NULL)) ;

    // Counts the edges for which its row node and col nodes have equal degree.
    // Simply, this edge would be removed from the graph once this node
    // is of degree < k. 
    GRB_TRY(GrB_vxm(
        edge_eq_deg, NULL, NULL, GxB_PLUS_ISEQ_UINT64, 
        degrees, edge_degrees, GrB_NULL)) ;


    // QUESTION: Is JIT faster or should I use accum as seen below?

    // If the nodes of an edge have different degrees, the edge is counted once.
    // If they have the same degree, that edge is double counted. So, we adjust:
    // edge_adjusted_deg = 2 * edge_gt_deg + edge_eq_deg 

    GRB_TRY(GrB_eWiseMult(
        edge_adjusted_deg, NULL, NULL, two_one, 
        edge_gt_deg, edge_eq_deg, GrB_NULL));
     */

    // If the nodes of an edge have different degrees, the edge is counted once.
    // If they have the same degree, that edge is double counted. So, we adjust:
    // edge_adjusted_deg = 2 * (x > y) + (x == y)
    GRB_TRY(GrB_vxm(
        edge_gt_deg, NULL, NULL, plus_2le, 
        degrees, edge_degrees, GrB_NULL)) ;

    GRB_TRY(GrB_Vector_nvals (&edge_vec_nvals, edge_adjusted_deg));
    vi_size = (edge_vec_nvals+1)*sizeof(GrB_Index);
    vx_size = (edge_vec_nvals+1)*sizeof(GrB_UINT64);

    GRB_TRY(GxB_Vector_unpack_CSC(
        edge_adjusted_deg, &index_edge, (void **) &edge_count_per_node,
        &vi_size,&vx_size,&iso,&edge_vec_nvals,NULL, GrB_NULL));
    LG_TRY (LAGraph_Free((void **)&index_edge, NULL)) ;

    GRB_TRY(GxB_Vector_unpack_CSC(
        degrees, &index_edge, (void **) &deg_arr,
        &vi_size,&vx_size,&iso,&edge_vec_nvals,NULL, GrB_NULL));

    LG_TRY(LAGraph_Free((void **)&index_edge, NULL)) ;

    //Build with degrees as indecies and handle duplicates via adition
    GRB_TRY(GrB_Vector_build (
        edges_per_deg, deg_arr, edge_count_per_node, n, GrB_PLUS_INT64)) ;

    // Start ones. Do I do this by building an iso vector? or by myself?
    // QUESTION: I want to do the following. I should duplicate deg_arr memory correct?
    GRB_TRY(GrB_Vector_build (
        cumulative_deg, deg_arr, ones, n, GrB_PLUS_INT64)) ;
    GRB_TRY(GxB_Vector_unpack_CSC(
        edges_per_deg, &index_edge, &edges_per_deg_arr,
        &vi_size, &vx_size, &iso, &edge_vec_nvals, NULL, GrB_NULL));
    
    //run a cummulative sum (backwards) on edges_per_deg_arr
    // for loops
    for(uint64_t i = edge_vec_nvals; i > 1; i--)
    {
        edges_per_deg_arr[i-1]+=edges_per_deg_arr[i];
    }
    //Construct a vector thats has edges_per_deg_arr values but repeated whenever 
    //index_edge has a skip and put into cumulative edges
    // ie. [0,1,6,7,10] & [9,7,3,2,0] ->
    // [9,7,7,7,7,7,3,2,2,2, . . .]

    uint64_t index = 0, i = 0;
    LG_TRY (LAGraph_Malloc((void **) &cumul_array, n, sizeof(uint64_t), msg)) ;

    for(; i < n; i++) // seems easily parrallelizable but idk if #pragma is enough.
    {
        cumul_array[i] = edges_per_deg_arr[index];
        if(i == index_edge[index])
        {
            index++ ;
            if(index == edge_vec_nvals || edges_per_deg_arr[index] == 0)
                break ;
        }
    }
    // QUESTION: 
    GRB_TRY(GxB_Vector_pack_Full(
        cumulative_edges, cumul_array, n * sizeof(uint64_t), false, NULL
    )) ;



    //GrB select or just do an if test
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}