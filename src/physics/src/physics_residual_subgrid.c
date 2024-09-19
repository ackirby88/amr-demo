/**
 * \file   physics_residual_subgrid.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief External solver function wrapper for calculating the subgrid residual 
 *        in each AMR quadrant/octant
 * 
 * Created on August 27, 2018, 11:00 AM
 */

/* header files */
#include "physics_residual_subgrid.h"


void physics_residual_subgrid(int *dim,int* nfields,int *nelem_subgrid,int *type,int *level,double *cell_geom,double *soln,double *res){
    
    /* external solver kernel */
    solver_residual_subgrid(dim,nfields,nelem_subgrid,type,level,cell_geom,soln,res);
    
}
