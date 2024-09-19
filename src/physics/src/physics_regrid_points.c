/**
 * \file   physics_regrid_points.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief External solver function wrapper 
 *        for user-specified point list regridding 
 * 
 * User is responsible for supplying list of points in memory as follows: 
 *      
 *      regrid_xyzh(x_pt1,  // x-coordinate of point #1
 *                  y_pt1,  // y-coordinate of point #1
 *                  z_pt1,  // z-coordinate of point #1
 *                  h_pt1,  //        width of point #1
 *                  x_pt2,  // x-coordinate of point #2
 *                  y_pt2,  // y-coordinate of point #2
 *                  z_pt2,  // z-coordinate of point #2
 *                  h_pt2,  //        width of point #2
 *                  .....)
 * 
 * Created on September 10, 2018, 8:32 AM
 */

/* header files */
#include "physics_regrid_points.h"


void physics_regrid_points(int *construct_grid,int *nregrid_pts,double *regrid_width,double *regrid_xyzh){
    
    //TODO: assign external function
    *construct_grid = 0;
    *nregrid_pts = 0;
    *regrid_width = 0.0;
    //regrid_xyzh = NULL;
    
}
