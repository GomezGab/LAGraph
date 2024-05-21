//------------------------------------------------------------------------------
// LAGraph_MaskColRow: remove certain rows or columns from a matrix
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

// Small utility that removes the entries in a a specified row or column of a 
// matrix 

#define LG_FREE_WORK                        \
{                                           \
    /* free any workspace used here */      \
                             \
}

#define LG_FREE_ALL                         \
{                                           \
    /* free any workspace used here */      \
    LG_FREE_WORK ;                          \
    /* free all the output variable(s) */   \
                                            \
    /* take any other corrective action */  \
}

#include "LG_internal.h"
#include "LAGraphX.h"

int LAGraph_HelloWorld // a simple algorithm, just for illustration
(
    // output
    GrB_Matrix *M,      // A - the rows and columns requested
    // input: not modified
    GrB_Vector cols,    // column "mask"
    GrB_Matrix A,       // Matrix to be pruned
    GrB_Vector rows,    // row "mask"
    bool cols_inverted, // true if you'd like to keep entries not in cols
    bool rows_inverted, // true if you'd like to keep entries not in rows
    char *msg
)
{
    //--------------------------------------------------------------------------
    // Declorations
    //--------------------------------------------------------------------------
    GrB_Vector cols_i = NULL; // inverse of cols
    GrB_Vector rows_i = NULL; // inverse of rows

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    LG_ASSERT (A != NULL && !(cols == NULL && rows == NULL), GrB_NULL_POINTER) ;
    if(cols != NULL)
    {
        if(cols_inverted)
        {
            
        }
    }
    if(rows != NULL)
    {

    }
    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}
