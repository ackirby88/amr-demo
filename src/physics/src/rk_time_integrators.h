/**
 * \file   rk_time_integrators.h
 * \author akirby
 *
 * Created on August 23, 2018, 12:46 PM
 */

#ifndef RK_TIME_INTEGRATORS_H
#define RK_TIME_INTEGRATORS_H


/* header files */
#include "physics_mpi.h"
#include "physics_rhs.h"
#include "physics_norm_residual.h"

#ifndef DOXYGEN_IGNORE
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#endif


//#include "rhs_solve.h"
//#include "iblank_and_norm_residual_external_solver.h"

/** Explicit Runge-Kutta 1st-Order Accurate Time Stepping Scheme (Forward Euler)
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
 * @param [in,out] Q            Solution
 * @param [out] k1              Residual Storage
 * @param [in]  dt              Time step
 * @param [in]  mpicomm_ext     MPI comminicator for mpi calls
 * @param [in]  counts          MPI array for how many cells to send
 * @param [in]  sendbuf         MPI array for sending data
 * @param [in]  request         MPI array for send/recv requests
 * @param [out] resnorm   Residual discrete norm (to check for nans)
 * @param [out] compute_time    Wall-clock time for computation
 * @param [out] residual_time   Wall-clock time for residual
 */
void rk_1(  int* dim,
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
            double* k1,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time);

/** Explicit Runge-Kutta 2nd-Order Accurate Time Stepping Scheme
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
 * @param [in,out] Q            Solution
 * @param [out] Q_stage         Temporary Solution Storage
 * @param [out] k1              Residual Storage #1
 * @param [out] k2              Residual Storage #2
 * @param [in]  dt              Time step
 * @param [in]  mpicomm_ext     MPI comminicator for mpi calls
 * @param [in]  counts          MPI array for how many cells to send
 * @param [in]  sendbuf         MPI array for sending data
 * @param [in]  request         MPI array for send/recv requests
 * @param [out] resnorm   Residual discrete norm (to check for nans)
 * @param [out] compute_time    Wall-clock time for computation
 * @param [out] residual_time   Wall-clock time for residual
 */
void rk_2(  int* dim,
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
            double* Q_stage,
            double* k1,
            double* k2,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time);

/** Explicit Runge-Kutta 3rd-Order Accurate Time Stepping Scheme
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
 * @param [in,out] Q            Solution
 * @param [out] Q_stage         Temporary Solution Storage
 * @param [out] k1              Residual Storage #1
 * @param [out] k2              Residual Storage #2
 * @param [out] k3              Residual Storage #3
 * @param [in]  dt              Time step
 * @param [in]  mpicomm_ext     MPI comminicator for mpi calls
 * @param [in]  counts          MPI array for how many cells to send
 * @param [in]  sendbuf         MPI array for sending data
 * @param [in]  request         MPI array for send/recv requests
 * @param [out] resnorm   Residual discrete norm (to check for nans)
 * @param [out] compute_time    Wall-clock time for computation
 * @param [out] residual_time   Wall-clock time for residual
 */
void rk_3(  int* dim,
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
            double* Q_stage,
            double* k1,
            double* k2,
            double* k3,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time);

/** Explicit Runge-Kutta 4th-Order Accurate Time Stepping Scheme
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
 * @param [in,out] Q            Solution
 * @param [out] Q_stage         Temporary Solution Storage
 * @param [out] k1              Residual Storage #1
 * @param [out] k2              Residual Storage #2
 * @param [out] k3              Residual Storage #3
 * @param [out] k4              Residual Storage #4
 * @param [in]  dt              Time step
 * @param [in]  mpicomm_ext     MPI comminicator for mpi calls
 * @param [in]  counts          MPI array for how many cells to send
 * @param [in]  sendbuf         MPI array for sending data
 * @param [in]  request         MPI array for send/recv requests
 * @param [out] resnorm   Residual discrete norm (to check for nans)
 * @param [out] compute_time    Wall-clock time for computation
 * @param [out] residual_time   Wall-clock time for residual
 */
void rk_4(  int* dim,
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
            double* Q_stage,
            double* k1,
            double* k2,
            double* k3,
            double* k4,
            double dt,
            MPI_Comm mpicomm_ext,
            double* sendbuf,
            MPI_Request* request,
            int* counts,
            double* resnorm,
            double* compute_time,
            double* residual_time);

/** Low Storage Explicit Runge-Kutta 4th-Order Accurate Time Stepping Scheme
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
 * @param [in,out] Q            Solution
 * @param [out] Q_stage         Temporary Solution Storage
 * @param [out] k1              Residual Storage
 * @param [in]  dt              Time step
 * @param [in]  mpicomm_ext     MPI comminicator for mpi calls
 * @param [in]  counts          MPI array for how many cells to send
 * @param [in]  sendbuf         MPI array for sending data
 * @param [in]  request         MPI array for send/recv requests
 * @param [out] resnorm   Residual discrete norm (to check for nans)
 * @param [out] compute_time    Wall-clock time for computation
 * @param [out] residual_time   Wall-clock time for residual
 */
void lserk_45(  int* dim,
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
                double* Q_stage,
                double* k1,
                double dt,
                MPI_Comm mpicomm_ext,
                double* sendbuf,
                MPI_Request* request,
                int* counts,
                double* resnorm,
                double* compute_time,
                double* residual_time);


#endif /* RK_TIME_INTEGRATORS_H */

