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
    
    // input: not modified
    LAGraph_Graph G,
    char *msg
)
{
    return (GrB_NOT_IMPLEMENTED) ;
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // Check inputs
    //--------------------------------------------------------------------------

    // Extract adjacency matrix to make incidence matrix - E

    // While less than nvals*swapsperedge swaps have happened

        // Coming into the loop I want E to be the incidence matrix of the new 
        // graph. W/o self edges nor parallel edges. Each row must have
        // exactly one 0 and one 1. They do not have to be randomly assigned.

        // Pair edges. Take a random permutation and pair adjacent values.
        // pairs wants:
        // rows: [0, 0, 1, 1, 2, 2, 3, 3, . . .] (ramp / 2)
        // cols: [ random permutation 1 - nvals] (LAGraph_Random_Seed + sort)
        // vals: [0, 3, 0, 2, 0, 2, 0, 3, . . .]
            // Start with [0, 2, 0, 2, 0, 2, . . .]
            // for i from 1 to nvals by 2s:
            //      vals[i] += vals[i] > vals[i-1]
        // CAUTION cols is likely smaller. Manage memory correctly.

        // each row of M will have 4 values [0,1,2,3]
        // pairs |  E  |  M
        //   0   | 0,1 | 0,1
        //   2   | 0,1 | 2,3
        //   3   | 0,1 | 3,2
        // This gives us randomization as to which type of swap is happening.
        // M = pairs * E (any_bitwisexor)

        // remove and entries of M with less than 3 values

        // probably worth it to put 0s on the interf diagonal to make the 
        // following vector dense.
        // interf <diagonal> = MT * M (plus_one)
        // reduce interf to a vector on rows
        // r_interf = interf * denseVec (max_first)
        // d_interf = diagonalize r_interf

        // set 3rd bit to 1 if it must be removed. 
        // An accumulator might work better.
        // M = d_interf * M (any_[z = 4 * (x > 2) | y]) 
        // M = extract M (x & 4 == 0)

        // M1 = extract M (x & 2 == 0)
        // M2 = extract M (x & 2 != 0)
        // Check if an edge already exists in the graph.
        // If a row of M1 or M2 has more than 1 vertex in common with a 
        // row of E, the value of that entry in the exists array will be 0. 
        // Otherwise, 1 or noval.
        // exists = E * M1T (zero_one)
        // exists += E * M2T (zero_one) (accum is one if theres an intersection)

        // 0 if intended swap already exists in the matrix, 1 or noval otherwise
        // r_exists = denseVec * exists (bitwiseAnd_second)
        // QUESTION: how to remove r_exists zeros from M matrix?

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
        

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
