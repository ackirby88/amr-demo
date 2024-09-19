/**
 * \file   amr_external.h
 * \author akirby
 *
 * Created on August 20, 2018, 11:19 AM
 */

#ifndef AMR_EXTERNAL_H
#define AMR_EXTERNAL_H


/* header files */
#include "var_ctx.h"
#include "var_quad.h"
#include "var_external.h"
#include "var_defines.h"
#include "amr_p4est_utilities.h"

#include <mpi.h>

#ifndef P4_TO_P8
#include <p4est_extended.h>
#else
#include <p8est_extended.h>
#endif


#ifdef __cplusplus
extern "C" {
#endif

/* external physics callback functions */
    

/** Wrapper to call external initialize solver
 * 
 * @param [in] rank                 mpi rank
 * @param [in] dim                  simulation dimension
 * @param [in] mesh_scale           coarsest element size
 * @param [in] nlevels              number of amr levels
 * @param [in] nfields              number of unknows at a solution point
 * @param [in] nfringe              number of ghost elements in 1D layer
 * @param [in] nsub_elem            number of 1D sub-elements in each amr quad
 * @param [in] extern_cell_size     number of info for each cell
 * @param [in] extern_face_size     number of info for each face
 * @param [in] extern_edge_size     number of info for each edge
 * @param [in] extern_node_size     number of info for each node
 * @param [in] extern_mpi_comm_size number of info for the mpi communication
 */
void hpc_amr_initialize_solver(int *rank,int *dim,double *mesh_scale,int *nlevels,
                               int *nfields,int *nfringe,int *nsub_elem,
                               int *extern_cell_size,int *extern_face_size,
                               int *extern_edge_size,int *extern_node_size,
                               int *extern_mpi_comm_size);

/** Wrapper to call external initialize quadrant data
 * 
 * @param [in]  dim             simulation dimension
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [out] soln            cell solution pointer
 */
void hpc_amr_initialize_quadrant(int *dim,int *unstructured,int *nfields,int *nsub_elem,
                                 int *cell_type,double *elem_h,double *elem_xyz,double *soln);

/** Wrapper to call external read regrid points for AMR 
 * 
 * @param [out] construct_grid  build the AMR grid to new set of points flag
 * @param [out] nregrid_pts     number of regrid points to be read
 * @param [out] regrid_width    buffer width for regridding around point
 * @param [out] regrid_xyzh     point coordinates to be adapted to (x,y,z,h)
 */
void hpc_amr_regrid_points(int *construct_grid,int *nregrid_pts,
                           double *regrid_width,double *regrid_xyzh);

/** Wrapper to call external tag amr features
 * 
 * @param [in]  dim             simulation dimension
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [in]  soln            cell solution pointer
 * @param [out] tag             amr tag flag (0-no tag, 1-tag)
 */
void hpc_amr_tag_feature(int *dim,int *unstructured,int *nfields,int *nsub_elem,int *cell_type,
                         double *elem_h,double *elem_xyz,double *soln,int *tag);

/** Wrapper to call external evolve solver
 * 
 * @param [in]  dim             simulation dimension
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  ncell           number of real cells
 * @param [in]  nghost          number of ghost cells
 * @param [in]  nface           number of faces
 * @param [in]  dof             number of real degrees of freedom
 * @param [in]  cell_info       cell info data structure
 * @param [in]  face_info       face info data structure
 * @param [in]  nsend           number of mpi sends
 * @param [in]  nrecv           number of mpi receives
 * @param [in]  send            mpi send info array
 * @param [in]  recv            mpi receive info array
 * @param [in]  sizes           mpi size per cell info array
 * @param [in]  geom            cell geometries
 * @param [in,out] soln         cell solution
 * @param [in,out] istep        simulation step counter
 * @param [out] regrid          regrid flag
 * @param [out] visualize       visualization flag
 * @param [out] checkpoint      checkpoint solution flag
 * @param [out] evolve_solution continue simulation flag
 * @param [in]  mpicomm_ext     mpi comminicator for mpi calls
 * @param [out] total_time      total simulation time
 * @param [out] compute_time    total computation time
 * @param [out] residual_time   total residual time
 */
void hpc_amr_evolve_solver(int* dim,int* unstructured,int* nfields,int* nsub_elem,
                           int* ncell,int* nghost,int* nface,int* dof,int* cell_info,
                           int* face_info,int* nsend,int* nrecv,int* send,int* recv,
                           int* sizes,double* geom,double* soln,
                           int* istep,int* regrid,int* visualize,int* checkpoint,
                           int* evolve_solution,MPI_Comm mpicomm_ext,double* total_time,
                           double* compute_time,double* residual_time);

/** Wrapper to call external 2D coarsen operator
 * 
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [out] uc1             fine cell #1 solution pointer
 * @param [out] uc2             fine cell #2 solution pointer
 * @param [out] uc3             fine cell #3 solution pointer
 * @param [out] uc4             fine cell #4 solution pointer
 * @param [out] soln            coarse cell solution pointer
 */
void hpc_amr_coarsen_operator_2d(int *unstructured,int *nfields,int *nsub_elem,
                                 int *cell_type,double *elem_h,double *elem_xyz,
                                 double *uc1,double *uc2,double *uc3,double *uc4,
                                 double *soln);

/** Wrapper to call external 3D coarsen operator
 * 
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [in]  uc1             fine cell #1 solution pointer
 * @param [in]  uc2             fine cell #2 solution pointer
 * @param [in]  uc3             fine cell #3 solution pointer
 * @param [in]  uc4             fine cell #4 solution pointer
 * @param [in]  uc5             fine cell #5 solution pointer
 * @param [in]  uc6             fine cell #6 solution pointer
 * @param [in]  uc7             fine cell #7 solution pointer
 * @param [in]  uc8             fine cell #8 solution pointer
 * @param [out] soln            coarse cell solution pointer 
 */
void hpc_amr_coarsen_operator_3d(int *unstructured,int *nfields,int *nsub_elem,
                                 int *cell_type,double *elem_h,double *elem_xyz,
                                 double *uc1,double *uc2,double *uc3,double *uc4,
                                 double *uc5,double *uc6,double *uc7,double *uc8,
                                 double *soln);

/** Wrapper to call external 2D refine operator
 * 
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [in]  soln            coarse cell solution pointer
 * @param [out] uc1             fine cell #1 solution pointer
 * @param [out] uc2             fine cell #2 solution pointer
 * @param [out] uc3             fine cell #3 solution pointer
 * @param [out] uc4             fine cell #4 solution pointer
 */
void hpc_amr_refine_operator_2d(int *unstructured,int *nfields,int *nsub_elem,
                                int *cell_type,double *elem_h,double *elem_xyz,
                                double *soln,
                                double *uc1,double *uc2,double *uc3,double *uc4);

/** Wrapper to call external 3D refine operator
 * 
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [in]  soln            coarse cell solution pointer
 * @param [out] uc1             fine cell #1 solution pointer
 * @param [out] uc2             fine cell #2 solution pointer
 * @param [out] uc3             fine cell #3 solution pointer
 * @param [out] uc4             fine cell #4 solution pointer
 * @param [out] uc5             fine cell #5 solution pointer
 * @param [out] uc6             fine cell #6 solution pointer
 * @param [out] uc7             fine cell #7 solution pointer
 * @param [out] uc8             fine cell #8 solution pointer
 */
void hpc_amr_refine_operator_3d(int *unstructured,int *nfields,int *nsub_elem,
                                int *cell_type,double *elem_h,double *elem_xyz,
                                double *soln,
                                double *uc1,double *uc2,double *uc3,double *uc4,
                                double *uc5,double *uc6,double *uc7,double *uc8);
    
#ifdef __cplusplus  
}
#endif



/** External wrapper function for initializing PHYSICS code module
 * 
 * @param [in]  d_ctx   simulation context data
 */
void external_initialize_solver(ctx_t *d_ctx);

/** External wrapper function for initializing quadradrant solution data
 * 
 * @param [in]  dim             simulation dimension
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quadrant
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [out] elem_soln       cell solution pointer
 */
void external_initialize_quadrant(int dim,int unstructured,int nfields,int nsub_elem,int *cell_type,double elem_h,double *elem_xyz,double *elem_soln);

/** External wrapper function of amr regrid points function
 * 
 * @param [out] construct_grid  build the AMR grid to new set of points flag
 * @param [out] nregrid_pts     number of regrid points to be read
 * @param [out] regrid_width    buffer width for regridding around point
 * @param [out] regrid_xyzh     point coordinates to be adapted to (x,y,z,h)
 */
void external_regrid_points(int *construct_grid,int *nregrid_pts,double *regrid_width,double *regrid_xyzh);

/** External wrapper function of amr tag feature function
 * 
 * @param [in]  dim             simulation dimension
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [in]  elem_soln       cell solution pointer
 * @param [out] tag             amr tag flag (0-no tag, 1-tag)
 */
void external_tag_feature(int dim,int unstructured,int nfields,int nsub_elem,int *cell_type,double elem_h,double *elem_xyz,double *elem_soln,int *tag);

/** External wrapper function of physics evolve function
 * 
 * @param [in]  d_ctx           simulation context data
 * @param [in]  p4est           p4est tree data structure
 * @param [in,out] istep        simulation step counter
 * @param [out] regrid          regrid flag
 * @param [out] visualize       visualization flag
 * @param [out] checkpoint      checkpoint solution flag
 * @param [out] evolve_solution continue simulation flag
 */
void external_evolve(ctx_t *d_ctx,p4est_t **p4est,int *istep,int *regrid,int *visualize,int *checkpoint,int *evolve_solution);

/** Initializes external solution storage and sets up communication patterns
 * 
 * @param [in] p4est            p4est tree data structure
 * @param [in] ghost            p4est ghost info 
 * @param [in] ghost_data       p4est ghost data
 */
void external_setup(p4est_t *p4est,p4est_ghost_t *ghost,quad_data_t *ghost_data);

/** Copies p4est quad data to external solution
 * 
 * @param [in] p4est            p4est tree data structure
 */
void external_p4est_to_soln(p4est_t *p4est);

/** Copies external solution to p4est quad data
 * 
 * @param [in,out] p4est        p4est tree data structure
 */
void external_soln_to_p4est(p4est_t *p4est);

/** Allocates external solution data storage
 * 
 * @param [in] p4est            p4est tree data structure
 * @param [in] ghost            p4est ghost info 
 * @param [in] ghost_data       p4est ghost data
 */
void external_allocate_soln(p4est_t *p4est,p4est_ghost_t *ghost,quad_data_t *ghost_data);

/** Deallocations external solution data storage
 * 
 * @param [in,out] p4est        p4est tree data structure
 */
void external_deallocate_soln(p4est_t *p4est);

/** Fills in mpi communication information
 * 
 * @param [in] p4est            p4est tree data structure
 * @param [in] ghost            p4est ghost info 
 */
void external_setup_ghost_exchange_data(p4est_t *p4est, p4est_ghost_t *ghost);

/** Callback function that fills in cell information
 * 
 * @param [in,out] info         p4est_iter_volume_info_t (quad vol info)
 * @param [in,out] user_data    user data passed with volume callback
 */
void external_cellinfo_callback(p4est_iter_volume_info_t *info, void *user_data);

/** Callback function that fills in face information
 * 
 * @param [in,out] info         p4est_iter_face_info_t (quad face info)
 * @param [in,out] user_data    user data passed with face callback
 */
void external_faceinfo_callback(p4est_iter_face_info_t *info, void *user_data);

/** Wrapper function for external coarsen operator
 * 
 * @param [in]  dim             simulation dimension
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [in]  uc              fine cell solution pointers
 * @param [out] soln            coarse cell solution pointer
 */
void external_coarsen_operator(int dim,int unstructured,int nfields,int nsub_elem,int *cell_type,double elem_h,double *elem_xyz,double **uc,double *soln);

/** Wrapper function for external refine operator
 * 
 * @param [in]  dim             simulation dimension
 * @param [in]  unstructured    unstructured grid flag
 * @param [in]  nfields         number of unknows at a solution point
 * @param [in]  nsub_elem       number of 1D sub-elements in each amr quad
 * @param [in]  cell_type       cell type
 * @param [in]  elem_h          cell length
 * @param [in]  elem_xyz        cell xyz geometry (z-ordered corners)
 * @param [in]  soln            coarse cell solution pointer
 * @param [out] uc              fine cell solution pointers
 */
void external_refine_operator(int dim,int unstructured,int nfields,int nsub_elem,int *cell_type,double elem_h,double *elem_xyz,double *soln,double **uc);


#endif /* AMR_EXTERNAL_H */

