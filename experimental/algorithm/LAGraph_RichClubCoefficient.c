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

// The values will be output as a sparce GrB_Vector, the rich club coefficient 
// of any value not in the vector is equivilant to the closest value above it.

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
    GrB_free (&cumulative_deg) ;            \
    GrB_free (&edges_per_deg) ;             \
    GrB_free (&two_one) ;                   \
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
//Is this best way to do this?
typedef void (*LAGraph_binary_function) (void *, const void *, const void *) ;

//I am hoping this gets done in a single leaq.
void two_one_add_uint64(uint64_t *z, const uint64_t *x, const uint64_t *y)
{ 
    (*z) = 2*(*x) + (*y);
}


int LAGraph_RichClubCoefficient // a simple algorithm, just for illustration
(
    // output
    GrB_Vector *rich_club_coefficents,    //rich_club_coefficents(i): rich club coefficents of i

    // input: not modified
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

    //contains the number of edges for which the ith node is
    //the smallest degree node in the pair
    GrB_Vector edge_gt_deg = NULL;

    //contains the number of edges for which the ith node is
    //the same degree as the other node in the pair
    //used to correct for undercounted nodes
    GrB_Vector edge_eq_deg = NULL;

    //Combines edge_gt_deg and edge_eq_deg to account for double counting in edge_eq_deg
    GrB_Vector edge_adjusted_deg = NULL;

    //the ith entry contains the number of nodes with degree greter than i.
    GrB_Vector cumulative_deg = NULL;

    //the ith entry contains the number of edges among nodes with degree greter than i.
    GrB_Vector edges_per_deg = NULL;

    GrB_BinaryOp two_one = NULL;
    GrB_Matrix A ;
    GrB_Index n ;

    // Grb_Index **edge_count_per_node = NULL;


    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rich_club_coefficents != NULL, GrB_NULL_POINTER);

    //double check this
    LG_ASSERT_MSG(G->kind == LAGraph_ADJACENCY_UNDIRECTED, GrB_INVALID_VALUE, "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->out_degree != NULL, GrB_EMPTY_OBJECT,"G->out_degree must be defined") ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE, "G->nself_edges must be zero") ;

    //--------------------------------------------------------------------------
    // Initializations
    //--------------------------------------------------------------------------
    A = G->A ;
    GRB_TRY(GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY(GrB_Matrix_new(&edge_degrees, GrB_UINT64,n,n)) ;
    GRB_TRY(GrB_Vector_new(&degrees, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edge_gt_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edge_eq_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&cumulative_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edges_per_deg, GrB_UINT64, n)) ;
    
    GRB_TRY(GxB_BinaryOp_new(&two_one, (LAGraph_binary_function) (&two_one_add_uint64), GrB_UINT64, GrB_UINT64, GrB_UINT64, "two_one_add_uint64",
                "void two_one_add_uint64(uint64_t *z, const uint64_t *x, const uint64_t *y) { (*z) = 2*(*x) + (*y); } ")) ;
            
    // code from LAGraph_MaximalIndependentSet
    // degrees = G->out_degree (check if this is needed) (I think no)
    GRB_TRY (GrB_assign (degrees, NULL, NULL, G->out_degree, GrB_ALL, n, NULL)) ;

    // code from LAGraph_SquareClustering
    // out_degrees as a diagonal matrix.
    #if LAGRAPH_SUITESPARSE
        #if GxB_IMPLEMENTATION >= GxB_VERSION (7,0,0)
        // SuiteSparse 7.x and later:
        GRB_TRY (GrB_Matrix_diag(&D, degrees, 0)) ;
        #else
        // SuiteSparse 6.x and earlier, which had the incorrect signature:
        GRB_TRY (GrB_Matrix_new(&D, GrB_INT64, n, n)) ;
        GRB_TRY (GrB_Matrix_diag(D, degrees, 0)) ;
        #endif
    #else
    // standard GrB:
    GRB_TRY (GrB_Matrix_diag(&D, degrees, 0)) ;
    #endif
    //Assigns each edge in the graph the value of their column node
    //IDK what the desctriptor is
    GRB_TRY (GrB_mxm(edge_degrees, NULL, NULL, GxB_ANY_SECOND_UINT64,A,D,GrB_DESC_T1)) ;
    
    //Adds up all the edges whose rows degrees are greater than their column degrees
    GRB_TRY(GrB_vxm(edge_gt_deg, NULL, NULL, GxB_PLUS_ISGT_UINT64, degrees, edge_degrees, GrB_DESC_T1)) ;

    //Adds up all the edges whose rows degrees are greater than their column degrees
    GRB_TRY(GrB_vxm(edge_eq_deg, NULL, NULL, GxB_PLUS_ISEQ_UINT64, degrees, edge_degrees, GrB_DESC_T1)) ;

    //Do I care if this is set intersection or union?
    GRB_TRY(GrB_eWiseMult(edge_adjusted_deg, NULL, NULL, two_one, edge_eq_deg, edge_gt_deg, GrB_DESC_T1));

    //GRB_TRY(GxB_Vector_unpack()) ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
