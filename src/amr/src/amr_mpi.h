/**
 * \file   amr_mpi.h
 * \author akirby
 *
 * Created on August 14, 2018, 4:54 PM
 */

#ifndef AMR_MPI_H
#define AMR_MPI_H


/* header files */
#ifndef P4_TO_P8
#include <p4est_bits.h>
#else
#include <p8est_bits.h>
#endif


/** MPI and p4est initialization function wrapper
 * 
 * @param [in]  argc        number of command line arguments
 * @param [in]  argv        command line arguments
 * @param [out] mpi_rank    mpi rank
 * @param [out] comm        mpi communicator
 */
int mpi_init(int argc, char** argv,int *mpi_rank,sc_MPI_Comm *comm);


/** MPI finalization function wrapper 
 * 
 * @return mpi finalize success
 */
int mpi_finalize();


#endif /* AMR_MPI_H */

