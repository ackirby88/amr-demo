/**
 * \file   physics_solve_interface.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Solve interface calls the numerical time stepping function:
 *        e.g. explicit time-accurate, implicit time-accurate,
 *             implicit steady-state.
 * 
 * Created on August 23, 2018, 1:02 PM
 */

/* header files */
#include "physics_solve_interface.h"


void physics_solve_interface(   int* dim,
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
                                double* soln,
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
                                double* residual_norm){
    
    
    
    physics_runge_kutta_interface(  dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,dof,
                                    cell_info,face_info,nsend,nrecv,send,recv,
                                    sizes,geom,soln,istep,regrid,visualize,checkpoint,
                                    mpicomm_ext,counts,sendbuf,request,
                                    step_count,t1,t2,compute_time,residual_time,
                                    residual_norm);
    
    
}

