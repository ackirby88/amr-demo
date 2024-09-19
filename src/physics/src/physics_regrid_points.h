/**
 * \file   physics_regrid_points.h
 * \author akirby
 *
 * Created on September 10, 2018, 8:31 AM
 */

#ifndef PHYSICS_REGRID_POINTS_H
#define PHYSICS_REGRID_POINTS_H


#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif


/** External solver function wrapper for user-specified point list regridding
 * 
 * @param [out] construct_grid      Flag indicator for constructing grid to new list of points
 * @param [out] nregrid_pts         Number of points in regrid point list
 * @param [out] regrid_width        Buffer width around each point for regridding 
 * @param [out] regrid_xyzh         Point coordinates and sizes array (x,y,z,h) 
 *                                  * See physics_regrid_points.c description for point ordering
 */
void physics_regrid_points(int *construct_grid,int *nregrid_pts,double *regrid_width,double *regrid_xyzh);


#endif /* PHYSICS_REGRID_POINTS_H */

