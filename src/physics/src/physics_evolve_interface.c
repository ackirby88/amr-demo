/**
 * \file   physics_evolve_interface.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Physics evolve interface that calls cpu or gpu physics_go function.
 * 
 * Created on August 22, 2018, 12:53 PM
 */

/* header files */
#include "physics_evolve_interface.h"


void physics_evolve_interface(  int* dim,
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
    
#if defined (_GPU_)
    physics_go_gpu( dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,dof,
                    cell_info,face_info,nsend,nrecv,send,recv,sizes,geom,soln,
                    istep,regrid,visualize,checkpoint,evolve_solution,
                    mpicomm_ext,total_time,compute_time,residual_time);
#else
    physics_go( dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,dof,
                cell_info,face_info,nsend,nrecv,send,recv,sizes,geom,soln,
                istep,regrid,visualize,checkpoint,evolve_solution,
                mpicomm_ext,total_time,compute_time,residual_time);
#endif
}
