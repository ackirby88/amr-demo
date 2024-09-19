/**
 * \file    hpc_amr_evolve_solver.c
 * \ingroup physics_group
 * \author  akirby
 * 
 * 
 * \brief AMR-Physics interface function for evolving (stepping) solver.
 *
 * Created on August 22, 2018, 12:40 PM
 */

/* header files */
#include "physics_evolve_interface.h"


void hpc_amr_evolve_solver( int* dim,
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
                            int* evolve_solution,
                            MPI_Comm mpicomm_ext,
                            double* total_time,
                            double* compute_time,
                            double* residual_time){
    
    physics_evolve_interface(   dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,dof,
                                cell_info,face_info,nsend,nrecv,send,recv,
                                sizes,geom,soln,istep,regrid,visualize,
                                checkpoint,evolve_solution,mpicomm_ext,
                                total_time,compute_time,residual_time);
    
}
