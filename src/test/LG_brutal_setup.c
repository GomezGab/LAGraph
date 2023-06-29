//------------------------------------------------------------------------------
// LG_brutal_setup.c: setup an LAGraph test with brutal memory testing
// -----------------------------------------------------------------------------

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

#include "LG_internal.h"
#include "LG_test.h"

#if LAGRAPH_SUITESPARSE
#if GxB_IMPLEMENTATION < GxB_VERSION (7,4,1)
#error "SuiteSparse:GraphBLAS v7.4.1 or later is required for LAGraph tests"
#endif
#endif

int LG_brutal_setup (char *msg)
{
    LG_brutal = -1 ;        // disable brutal testing for now
    LG_nmalloc = 0 ;        // assuming nothing is malloc'd
    int result = LAGr_Init (GrB_NONBLOCKING,
        LG_brutal_malloc, LG_brutal_calloc,
        LG_brutal_realloc, LG_brutal_free, msg) ;
    return (result) ;
}

