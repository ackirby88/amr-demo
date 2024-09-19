/** 
 * \file   physics_evolve_interface.h
 * \author akirby
 *
 * Created on August 22, 2018, 12:57 PM
 */

#ifndef PHYSICS_EVOLVE_INTERFACE_H
#define PHYSICS_EVOLVE_INTERFACE_H

/* header files */
#ifndef DOXYGEN_IGNORE
#include <mpi.h>
#endif

#if defined (_GPU_)
#include "physics_go_gpu.h"
#else
#include "physics_go.h"
#endif


#ifdef __cplusplus
extern "C" {
#endif

    
/** Physics evolve interface function
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
 * @param [in,out] soln         Cell solution
 * @param [in,out] istep        Simulation step counter
 * @param [out] regrid          Regrid flag
 * @param [out] visualize       Visualization flag
 * @param [out] checkpoint      Checkpoint solution flag
 * @param [out] evolve_solution Continue simulation flag
 * @param [in]  mpicomm_ext     MPI comminicator for mpi calls
 * @param [out] total_time      Total simulation time (computation+communication)
 * @param [out] compute_time    Total computation time (update+residual)
 * @param [out] residual_time   Total residual time
 */
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
                                double* residual_time);


#ifdef __cplusplus
}
#endif

#endif /* PHYSICS_EVOLVE_INTERFACE_H */

