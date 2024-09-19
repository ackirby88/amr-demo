/**
 * \file   hpc_amr_initialize_solver.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief AMR-Physics interface function for initializing 
 *        physics solver constants.
 * 
 * Created on August 21, 2018, 4:05 PM
 */

/* header files */
#include "physics_initialize.h"


void hpc_amr_initialize_solver( int *rank,int *dim,double *mesh_scale,int *nlevels,
                                int *nfields,int *nfringe,int *nelem_subgrid,
                                int *extern_cell_size,int *extern_face_size,
                                int *extern_edge_size,int *extern_node_size,
                                int *extern_mpi_comm_size){
    
    physics_initialize( *rank,*dim,*mesh_scale,*nlevels,
                        *nfields,*nfringe,*nelem_subgrid,
                        *extern_cell_size,*extern_face_size,
                        *extern_edge_size,*extern_node_size,
                        *extern_mpi_comm_size);
    
}
