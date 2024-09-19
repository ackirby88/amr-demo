/**
 * \file   var_grid.h
 * \ingroup amr_group
 * \author akirby
 *
 * \brief   Grid/mesh related data related to AMR setup and regridding: 
 *          grid dimensions & coordinates, amr levels, 
 *          box refinement coordinates, external regridding point lists.
 * 
 * Created on August 14, 2018, 5:20 PM
 */

#ifndef VAR_GRID_H
#define VAR_GRID_H


/* header files */
#include "var_defines.h"


/**
 * grid_t contains amr grid related data.
 */
typedef struct {
    
    double scale;               /**< Scaling factor of the coarsest element */
    double min_dx;              /**< Size of the finest element */
    double xlo[3];              /**< Lower coordinates of the grid domain */
    double xhi[3];              /**< Upper coordinates of the grid domain */
    int dim;                    /**< Spatial dimension of the grid */
    int max_level;              /**< Maximum level of grid refinement */
    int min_level;              /**< Minumum level of grid refinement */
    int nelem_subgrid;          /**< Number of elements in each subgrid */
    int regrid_interval;        /**< Interval at which the simulation regrids */
    int periodic[3];            /**< Periodic boundary condition indicator */
    int unstructured_flag;      /**< Unstructured grid indicator */
    char unst_grid_file[BUFF_SIZE]; /**< Unstructured grid file name */
    
    int refined_box;            /**< Initialized a refined box region */
    double box_xlo[3];          /**< Lower coordinates of refined box */
    double box_xhi[3];          /**< Upper coordinates of refined box */
    
    int construct_grid;         /**< Initialize grid by building up levels */
    int nregrid_pts;            /**< Number of regrid points */
    double regrid_width;        /**< Regrid width around points */
    double *regrid_xyzh;         /**< Coordinates of regrid points */
    
}
grid_t; /**< Internal data type for amr grid data */


#endif /* VAR_GRID_H */

