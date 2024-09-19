/**
 * \file   var_ctx.h
 * \ingroup amr_group
 * \author akirby
 * 
 * \brief   Simulation context data related to other data types:
 *          mpi, visualization, grid, external solver, solution, 
 *          and initialization. 
 *
 * Created on August 15, 2018, 1:05 PM
 */

#ifndef VAR_CTX_H
#define VAR_CTX_H

/* header files */
#include "var_mpi.h"
#include "var_vtk.h"
#include "var_grid.h"
#include "var_solution.h"
#include "var_external.h"
#include "var_initialize.h"

/**
 * ctx_t contains simulation context data: 
 *    * [mpi_t] d_mpi: mpi context data
 *    * [vtk_t] d_vtk: visualization context data
 *    * [grid_t] d_grid: grid context data
 *    * [external_t] d_external: external communication & solver context data
 *    * [solution_t] d_solution: solution checkpoint context data
 *    * [initialize_t] d_initialize: user initialization context data
 * 
 */
typedef struct {
    
    mpi_t d_mpi;                        /**< MPI data */
    vtk_t d_vtk;                        /**< Visualization data */
    grid_t d_grid;                      /**< Grid data */
    external_t d_external;              /**< External Communication data */
    solution_t d_solution;              /**< Solution data */
    initialize_t d_initialize;          /**< User Input data */
    
    int istep;                          /**< Simulation time step */
    double total_time;                  /**< Simulation wall-clock time */
    double compute_time;                /**< Compute wall-clock time (no communication)*/
    double residual_time;               /**< Residual wall-clock time (no communication,no update)*/
        
}
ctx_t;  /**< Internal data type for simulation context related data*/

#endif /* VAR_CTX_H */

