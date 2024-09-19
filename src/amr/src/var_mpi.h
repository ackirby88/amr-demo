/**
 * \file   var_mpi.h
 * \ingroup amr_group
 * \author akirby
 *
 * \brief   MPI related data: 
 *          mpi rank and communicator.
 * 
 * Created on August 14, 2018, 5:02 PM
 */

#ifndef VAR_MPI_H
#define VAR_MPI_H


/* header files */
#include "main.h"

#ifndef P4_TO_P8
#include <p4est_bits.h>
#else
#include <p8est_bits.h>
#endif


/**
 * mpi_t contains MPI related data.
 */
typedef struct {
    
    int rank;               /**< MPI rank */
    sc_MPI_Comm comm;       /**< MPI communicator */
    
}
mpi_t; /**< Internal data type for mpi data */


#endif /* VAR_MPI_H */

