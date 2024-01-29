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
    GrB_free (&reorder_g) ;                 \
    GrB_free (&D) ;                         \
    GrB_free (&deg) ;                       \
    GrB_free (&edge_gp_deg) ;               \
    GrB_free (&edge_eq_deg) ;               \
    GrB_free (&cumulative_deg) ;            \
    GrB_free (&edges_per_deg) ;             \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
    GrB_free (&rich_club_coefficents) ;     \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_RichClubCoefficient // a simple algorithm, just for illustration
(
    // output
    GrB_Vector *rich_club_coefficents,    //rich_club_coefficents(i): rich club coefficents of i

    // input: not modified
    LAGraph_Graph G,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Matrix edge_degrees = NULL;
    GrB_Matrix reorder_g = NULL;
    GrB_Matrix D = NULL;

    GrB_Vector deg = NULL;
    GrB_Vector edge_gp_deg = NULL;
    GrB_Vector edge_eq_deg = NULL;
    GrB_Vector cumulative_deg = NULL;
    GrB_Vector edges_per_deg = NULL;


    LG_CLEAR_MSG ;                      // clears the msg string, if not NULL

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
