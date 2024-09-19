/**
 * \file   physics_var_interface.h
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief   AMR-PHYSICS interface related data: 
 *          ghost fringe size, number of sub-elements in each amr quadrant, 
 *          number of fields (unknowns) per point, and external 
 *          communuication constants.
 *
 * Created on August 15, 2018, 2:32 PM
 */

#ifndef VAR_INTERFACE_H
#define VAR_INTERFACE_H

/** 
 * physics_interface_t contains amr-physics interface data 
 */
typedef struct {
    
    int dim;                    /**< Simulation spatial dimension */
    int nfringe;                /**< Number of ghost cell layers */
    int nfields;                /**< Number of fields (unknowns) at each grid point */
    int nelem_subgrid;          /**< Number of 1D subgrid elements in each AMR quadrant */
    
    int extern_cell_size;       /**< Number of entries per cell in cell_info */
    int extern_face_size;       /**< Number of entries per face in face_info */
    int extern_edge_size;       /**< Number of entries per face in egde_info */
    int extern_node_size;       /**< Number of entries per face in node_info */
    int extern_mpi_comm_size;   /**< Number of entries in the extern comm mpi */
    
    int nlevels;                /**< Number of AMR levels*/
    double coarse_grid_h;       /**< Length of coarsest grid element */
    
}
physics_interface_t; /**< data type for the amr-physics interface data */

physics_interface_t d_physics_interface; /**< Primary physics interface storage */

#endif /* VAR_INTERFACE_H */

