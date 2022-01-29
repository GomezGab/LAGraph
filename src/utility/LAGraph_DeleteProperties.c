//------------------------------------------------------------------------------
// LAGraph_DeleteProperties: deletes the cached properties of a graph
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#include "LG_internal.h"

int LAGraph_DeleteProperties    // returns 0 if successful, -1 if failure
(
    LAGraph_Graph G,        // G stays valid, only cached properties are freed
    char *msg
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    if (G == NULL)
    {
        // success: nothing to do
        return (GrB_SUCCESS) ;
    }

    //--------------------------------------------------------------------------
    // free all cached properties of the graph
    //--------------------------------------------------------------------------

    GrB_TRY (GrB_free (&(G->AT))) ;
    GrB_TRY (GrB_free (&(G->rowdegree))) ;
    GrB_TRY (GrB_free (&(G->coldegree))) ;

    //--------------------------------------------------------------------------
    // clear the scalar properties of the graph
    //--------------------------------------------------------------------------

    G->A_structure_is_symmetric = LAGRAPH_UNKNOWN ;
    G->ndiag = LAGRAPH_UNKNOWN ;
    return (GrB_SUCCESS) ;
}
