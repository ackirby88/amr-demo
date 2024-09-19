/**
 * \file   physics_mpi.h
 * \author akirby
 *
 * Created on August 22, 2018, 3:09 PM
 */

#ifndef MPI_REAL_TO_GHOST_EXTERNAL_H
#define MPI_REAL_TO_GHOST_EXTERNAL_H

/* header files */
#include "physics_var_interface.h"

#ifndef DOXYGEN_IGNORE
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif



#ifdef __cplusplus
}
#endif


/** Initializes the MPI buffers.
 * 
 * @param [out] sendbuf         MPI send buffer storage
 * @param [out] request         MPI request array
 * @param [out] counts          Number of cells to send/recv array for each mpi rank
 * @param [in] ncell            Number of real cells
 * @param [in] nghost           Number of ghost cells
 * @param [in] nsend            Number of mpi sends
 * @param [in] nrecv            Number of mpi receives
 * @param [in] send             MPI send info array
 * @param [in] sizes            MPI size per cell info array
 */
void physics_mpi_initialize_buffers(double* sendbuf,MPI_Request* request,int* counts,int ncell,int nghost,int nsend,int nrecv,int* send,int* sizes);

/** Deallocates the MPI buffers.
 * 
 * @param [in,out] sendbuf
 * @param [in,out] request
 * @param [in,out] counts
 */
void physics_mpi_free_buffers(double* sendbuf,MPI_Request* request, int* counts);

/** MPI exchange of solution data.
 * 
 * @param [in] nsend            Number of mpi sends
 * @param [in] nrecv            Number of mpi receives
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] send             MPI send info array
 * @param [in] recv             MPI receive info array
 * @param [in] sizes            MPI size per cell info array
 * @param [in,out] vec          Solution data & MPI recv buffer (after dof index)
 * @param [in] mpicomm_ext      MPI communicator
 * @param [in] sendbuf          MPI send buffer storage
 * @param [in] request          MPI request array
 * @param [in] counts           Number of cells to send/recv array for each mpi rank
 */
void physics_mpi_real_to_ghost(int nsend,int nrecv,int dof,int* send,int* recv,int* sizes,double* vec,MPI_Comm mpicomm_ext,double* sendbuf,MPI_Request* request,int* counts);

#endif /* MPI_REAL_TO_GHOST_EXTERNAL_H */

