/**
 * \file   amr_vtk.h
 * \author akirby
 *
 * Created on August 20, 2018, 5:24 PM
 */

#ifndef AMR_VTK_H
#define AMR_VTK_H

/* header files */
#include "var_defines.h"
#include "var_quad.h"
#include "var_ctx.h"

#ifndef P4_TO_P8
#include <p4est_vtk.h>
#include <p4est_extended.h>
#else
#include <p8est_vtk.h>
#include <p8est_extended.h>
#endif


/** Visualize function for just the AMR grid
 * 
 * @param [in] p4est        p4est forest structure
 */
void vtk_visualize_grid(p4est_t *p4est);

/** Visualize function for the AMR grid and solution
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] timestep     output time step number
 */
void vtk_write_solution(p4est_t *p4est,int timestep);

/** Interpolation function of the quadrant data to the corners (external)
 * 
 * @param [in] info             p4est volume information related to quadrant
 * @param [in,out] user_data    interpolated data
 */
void interpolate_solution(p4est_iter_volume_info_t *info, void *user_data);


#endif /* AMR_VTK_H */

