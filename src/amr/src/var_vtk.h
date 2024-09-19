/**
 * \file   var_vtk.h
 * \ingroup amr_group
 * \author akirby
 *
 * \brief Visualization related data: visualization internal, visualization data types.
 * 
 * Created on August 14, 2018, 5:20 PM
 */

#ifndef VAR_VTK_H
#define VAR_VTK_H


/**
 * vtk_t contains visualization related data
 */
typedef struct {
    
    int visualization_interval; /**< Interval to visualize solution */
    
}
vtk_t; /**< Internal data type for visualization data */


/**
 * vtk_plot_fluid_t contains conservative fluid variable visualizatio related data
 */
typedef struct {
    
    sc_array_t *rho;        /**< Density */
    sc_array_t *rhou;       /**< Momentum (x-direction) */
    sc_array_t *rhov;       /**< Momentum (y-direction) */
    sc_array_t *rhow;       /**< Momentum (z-direction) */
    sc_array_t *rhoe;       /**< Energy Density */
    
}
vtk_plot_fluid_t; /**< Visualization Data Type: conservative fluid variables */


#endif /* VAR_VTK_H */

