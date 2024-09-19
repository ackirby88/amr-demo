/**
 * \file   hpc_amr_regrid_points.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief AMR-Physics interface function for reading user-specified list of 
 *        points to regrid AMR grid.
 * 
 * Created on September 10, 2018, 8:27 AM
 */


/* header files */
#include "physics_regrid_points.h"


void hpc_amr_regrid_points( int *construct_grid,
                            int * nregrid_pts,
                            double *regrid_width,
                            double *regrid_xyzh){
    
    physics_regrid_points(construct_grid,nregrid_pts,regrid_width,regrid_xyzh);
    
}
