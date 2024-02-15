//------------------------------------------------------------------------------
// LAGraph_RichClubCoefficient: a nearly empty algorithm
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

// This is a bare-bones "algorithm" that does nearly nothing all.  It simply
// illustrates how a new algorithm can be added to the experimental/algorithm
// folder.  All it does is make a copy of the G->A matrix and return it as
// the new matrix Y.  Inside, it allocates some worspace as well (the matrix W,
// which is not used).  To illustrate the use of the error msg string, it
// returns an error if the graph not directed.

// The GRB_TRY and LG_TRY macros use the LG_FREE_ALL macro to free all
// workspace and all output variables if an error occurs.  To use these macros,
// you must define the variables before using them, or before using GRB_TRY
// or LG_TRY.  The LG_TRY macro is defined in src/utility/LG_internal.h.

// To create your own algorithm, create a copy of this file, rename it
// to LAGraph_whatever.c, and use it as a template for your own algorithm.
// Then place the prototype in include/LAGraphX.h.

// See experimental/test/test_HelloWorld.c for a test for this method, and
// experimental/benchmark/helloworld_demo.c and helloworld2_demo.c for two
// methods that benchmark the performance of this algorithm.

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
    // check inputs
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

    //the ith entry contains the number of nodes with degree greter than i.
    GrB_Vector cumulative_deg = NULL;

    //the ith entry contains the number of edges among nodes with degree greter than i.
    GrB_Vector edges_per_deg = NULL;

    GrB_Matrix A ;
    GrB_Index n ;

    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    LG_ASSERT (rich_club_coefficents != NULL, GrB_NULL_POINTER);

    //double check this
    LG_ASSERT_MSG(G->kind == LAGraph_ADJACENCY_UNDIRECTED, GrB_INVALID_VALUE, "G->A must be symmetric") ;
    LG_ASSERT_MSG (G->out_degree != NULL, GrB_EMPTY_OBJECT,"G->out_degree must be defined") ;
    LG_ASSERT_MSG (G->nself_edges == 0, GrB_INVALID_VALUE, "G->nself_edges must be zero") ;
    A = G->A ;
    GRB_TRY(GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY(GrB_Matrix_new(&edge_degrees, GrB_UINT64,n,n)) ;
    GRB_TRY(GrB_Vector_new(&degrees, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edge_gt_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edge_eq_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&cumulative_deg, GrB_UINT64, n)) ;
    GRB_TRY(GrB_Vector_new(&edges_per_deg, GrB_UINT64, n)) ;

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
    GRB_TRY(GrB_vxm(edge_gt_deg, NULL, NULL, GxB_PLUS_ISGT_UINT64, D, edge_degrees, GrB_DESC_T1)) ;

    //Adds up all the edges whose rows degrees are greater than their column degrees
    GRB_TRY(GrB_vxm(edge_eq_deg, NULL, NULL, GxB_PLUS_ISEQ_UINT64, D, edge_degrees, GrB_DESC_T1)) ;


    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
