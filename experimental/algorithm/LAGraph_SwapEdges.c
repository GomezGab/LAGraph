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
    GrB_free (&W) ;                         \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (&Y) ;                         \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_HelloWorld // a simple algorithm, just for illustration
(
    // output
    
    // input: not modified
    LAGraph_Graph G,
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    // Extract adjacency matrix to make incidence matrix - E

    // While less than nvals*swapsperedge swaps have happened
        // Pair edges into M1 and M2 matricies (HOW?)
        // M <-- M1 + M2
        // Get MT
        // remove and entries of M with less than 3 values
        // Swap interference <-- MT plus_one M
        // Extract >=2 masking the diagonal
        // stop any swaps with a row in Swap interference w/ degree > 1
        // do remaining swaps. (HOW?)
        // update E
        

    LG_FREE_WORK ;
    (*Yhandle) = Y ;
    return (GrB_SUCCESS) ;
}
