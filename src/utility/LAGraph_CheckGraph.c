//------------------------------------------------------------------------------
// LAGraph_CheckGraph: check if a graph is valid
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_CheckGraph      // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // graph to check
    char *msg
)
{

    //--------------------------------------------------------------------------
    // clear the msg and check basic components
    //--------------------------------------------------------------------------

    LG_CHECK_INIT (G, msg) ;
    GrB_Matrix A = G->A ;
    LAGraph_Kind kind = G->kind ;

    //--------------------------------------------------------------------------
    // ensure the matrix is square for directed or undirected graphs
    //--------------------------------------------------------------------------

    GrB_Index nrows, ncols ;
    if (kind == LAGRAPH_ADJACENCY_UNDIRECTED ||
        kind == LAGRAPH_ADJACENCY_DIRECTED)
    {
        GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
        GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
        LG_CHECK (nrows != ncols, -1, "adjacency matrix invalid") ;
    }

#if defined(GxB_SUITESPARSE_GRAPHBLAS)
    // only by-row format is supported
    GxB_Format_Value fmt ;
    GrB_TRY (GxB_get (A, GxB_FORMAT, &fmt)) ;
    LG_CHECK (fmt != GxB_BY_ROW, -2, "only by-row format supported") ;
#endif

    //--------------------------------------------------------------------------
    // check the cached properties
    //--------------------------------------------------------------------------

    GrB_Matrix AT = G->AT ;
    if (AT != NULL)
    {
        GrB_Index nrows2, ncols2;
        GrB_TRY (GrB_Matrix_nrows (&nrows2, AT)) ;
        GrB_TRY (GrB_Matrix_ncols (&ncols2, AT)) ;
        LG_CHECK (nrows != ncols2 || ncols != nrows2, -3,
            "G->AT matrix invalid") ;

#if defined(GxB_SUITESPARSE_GRAPHBLAS)
        // only by-row format is supported
        GxB_Format_Value fmt ;
        GrB_TRY (GxB_get (AT, GxB_FORMAT, &fmt)) ;
        LG_CHECK (fmt != GxB_BY_ROW, -4, "only by-row format supported") ;
#endif

        // ensure the types of A and AT are the same
        LG_CHECK (G->A_type != G->AT_type, -5, "A and AT types are different") ;
    }

    GrB_Vector rowdegree = G->rowdegree ;
    if (rowdegree != NULL)
    {
        GrB_Index m ;
        GrB_TRY (GrB_Vector_size (&m, rowdegree)) ;
        LG_CHECK (m != nrows, -6, "rowdegree invalid size") ;
        LG_CHECK (G->rowdegree_type != GrB_INT64, -7,
                  "rowdegree has wrong type") ;
    }

    GrB_Vector coldegree = G->coldegree ;
    if (coldegree != NULL)
    {
        GrB_Index n ;
        GrB_TRY (GrB_Vector_size (&n, coldegree)) ;
        LG_CHECK (n != ncols, -8, "coldegree invalid size") ;
        LG_CHECK (G->coldegree_type != GrB_INT64, -9,
                  "coldegree has wrong type") ;
    }

    return (0) ;
}