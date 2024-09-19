/**
 * \file   var_defines.h
 * \ingroup amr_group
 * \author akirby
 * 
 * \brief   Global variable data related to simulation static constants:
 *          ghost fringe size, number of sub-elements in each amr quadrant, 
 *          number of fields (unknowns) per point, and external 
 *          communuication constants.
 * 
 * Created on August 15, 2018, 12:11 PM
 */

#ifndef VAR_DEFINES_H
#define VAR_DEFINES_H

/* header files */
#include "main.h"

//DO NOT CHANGE UNLESS YOU IMPLEMENT NEW DATA
#define BUFF_SIZE 256           /**< Number of string characters allowed      */
#define EXTERN_COMM_MPI_SIZE  2 /**< Number of entries in the extern comm mpi */
#define EXTERN_COMM_CELL_SIZE 3 /**< Number of entries per cell in cell_info  */
#define EXTERN_COMM_FACE_SIZE 11/**< Number of entries per face in face_info  */
#define EXTERN_COMM_EDGE_SIZE 0 /**< Number of entries per edge in edge_info  */
#define EXTERN_COMM_NODE_SIZE 0 /**< Number of entries per node in node_info  */

// Modify to your problem
#define NFRINGE 2   /**< Number of ghost layers required for communication    */

//DO NOT CHANGE
#define NXPATCH NFRINGE /**< Number of patch elements equal to fringe depth   */

// Modify to your problem
#ifndef P4_TO_P8
#define NFIELDS 4   /**< Number of field variables (unknowns) per grid point in 2D */
#else
#define NFIELDS 5   /**< Number of field variables (unknowns) per grid point in 3D */
#endif


#endif /* VAR_DEFINES_H */

