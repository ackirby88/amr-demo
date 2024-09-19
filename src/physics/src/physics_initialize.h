/**
 * \file   physics_initialize.h
 * \author akirby
 *
 * Created on August 21, 2018, 2:55 PM
 */

#ifndef PHYSICS_INITIALIZE_H
#define PHYSICS_INITIALIZE_H

/* header files */
#ifndef DOXYGEN_IGNORE
#include <unistd.h>
#endif

#include "physics_utilities.h"
#include "physics_var_physics.h"
#include "physics_var_simulation.h"
#include "physics_var_interface.h"


/** Initialize physics constants from AMR code module.
 * 
 * @param [in] rank                 MPI rank
 * @param [in] dim                  Simulation spatial dimension
 * @param [in] mesh_scale           Coarse grid element length
 * @param [in] nlevels              Number of AMR levels
 * @param [in] nfields              Number of fields (unknowns) per grid point
 * @param [in] nfringe              Number of ghost layers
 * @param [in] nelem_subgrid        Number of 1D subgrid elements in AMR quad
 * @param [in] extern_cell_size     Number of fields in cell_info data structure
 * @param [in] extern_face_size     Number of fields in face_info data structure
 * @param [in] extern_edge_size     Number of fields in edge_info data structure
 * @param [in] extern_node_size     Number of fields in node_info data structure
 * @param [in] extern_mpi_comm_size Number of fields in mpi_info data structure
 */
void physics_initialize(int rank,int dim,double mesh_scale,int nlevels,
                        int nfields,int nfringe,int nelem_subgrid,
                        int extern_cell_size,int extern_face_size,
                        int extern_edge_size,int extern_node_size,
                        int extern_mpi_comm_size);

/** Read user-specified input file for physics simulation data.
 * 
 * @param [in] filename     Input file name
 * @param [in] noinput      Input file found flag
 */
void physics_initialize_inputs(char *filename, int noinput);

/** Write inputs read from user-specified input file for physics simulation data.
 * 
 * @param [in] filename     Input file name
 * @param [in] rank         MPI rank
 * @param [in] noinput      Input file found flag
 */
void physics_initialize_save_inputs(char *filename, int rank, int noinput);

#endif /* PHYSICS_INITIALIZE_H */

