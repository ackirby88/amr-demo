/**
 * \file   physics_runge_kutta_interface.h
 * \author akirby
 *
 * Created on August 24, 2018, 2:50 PM
 */

#ifndef PHYSICS_RUNGE_KUTTA_INTERFACE_H
#define PHYSICS_RUNGE_KUTTA_INTERFACE_H


/* header files */
#include "rk_time_integrators.h"
#include "physics_time_step.h"
#include "physics_var_simulation.h"


/**  Explicit Runge-Kutta time stepping function interface.
 * 
 * @param [in]  dim             Simulation dimension
 * @param [in]  unstructured    Unstructured grid flag
 * @param [in]  nfields         Number of unknows at a solution point
 * @param [in]  nelem_subgrid   Number of 1D sub-elements in each amr quad
 * @param [in]  ncell           Number of real cells
 * @param [in]  nghost          Number of ghost cells
 * @param [in]  nface           Number of faces
 * @param [in]  dof             Number of real degrees of freedom
 * @param [in]  cell_info       Cell info data structure
 * @param [in]  face_info       Face info data structure
 * @param [in]  nsend           Number of mpi sends
 * @param [in]  nrecv           Number of mpi receives
 * @param [in]  send            MPI send info array
 * @param [in]  recv            MPI receive info array
 * @param [in]  sizes           MPI size per cell info array
 * @param [in]  geom            Cell geometries
 * @param [in,out] Q         Cell solution
 * @param [in,out] istep        Simulation step counter
 * @param [out] regrid          Regrid flag
 * @param [out] visualize       Visualization flag
 * @param [out] checkpoint      Checkpoint solution flag
 * @param [in]  mpicomm_ext     MPI comminicator for mpi calls
 * @param [in]  counts          MPI array for how many cells to send
 * @param [in]  sendbuf         MPI array for sending data
 * @param [in]  request         MPI array for send/recv requests
 * @param [out] step_count      Number of solve steps performed
 * @param [out] t1              Start wall-clock time for solve
 * @param [out] t2              Stop wall-clock time for solve
 * @param [out] compute_time    Wall-clock time for computation
 * @param [out] residual_time   Wall-clock time for residual
 * @param [out] residual_norm   Residual discrete norm (to check for nans)
 */
void physics_runge_kutta_interface( int* dim,
                                    int* unstructured,
                                    int* nfields,
                                    int* nelem_subgrid,
                                    int* ncell,
                                    int* nghost,
                                    int* nface,
                                    int* dof,
                                    int* cell_info,
                                    int* face_info,
                                    int* nsend,
                                    int* nrecv,
                                    int* send,
                                    int* recv,
                                    int* sizes,
                                    double* geom,
                                    double* Q,
                                    int* istep,
                                    int* regrid,
                                    int* visualize,
                                    int* checkpoint,
                                    MPI_Comm mpicomm_ext,
                                    int* counts,
                                    double* sendbuf,
                                    MPI_Request* request,
                                    int* step_count,
                                    double* t1,
                                    double* t2,
                                    double* compute_time,
                                    double* residual_time,
                                    double* residual_norm);

#endif /* PHYSICS_RUNGE_KUTTA_INTERFACE_H */

