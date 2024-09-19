/**
 * \file   physics_residual_bc_unstructured.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief External solver function wrapper for boundary condition residual 
 *        calculation on unstructured grids
 * 
 * Created on August 31, 2018, 11:53 AM
 */

/* header files */
#include "physics_residual_bc_unstructured.h"


void physics_residual_bc_unstructured(int *dim,int *nfields,int *nelem_subgrid,int *side,
                                      int *type,int *bc,double *geom,double *soln,double *res){

    /* external solver kernel */
    solver_residual_bc_unstructured(dim,nfields,nelem_subgrid,side,type,bc,geom,soln,res);
    
}
