/**
 * \file   physics_residual_bc_structured.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief External solver function wrapper for boundary condition residual 
 *        calculation on structured grids
 * 
 * Created on September 4, 2018, 9:34 AM
 */

/* header files */
#include "physics_residual_bc_structured.h"


void physics_residual_bc_structured(int *dim,int *nfields,int *nelem_subgrid,
                                    int *side,int *type,int *level,int *bc,
                                    double *h_base,double *xlo,double *ylo,double *zlo,
                                    double *soln,double *res){

    /* external solver kernel */
    solver_residual_bc_structured(dim,nfields,nelem_subgrid,side,type,level,bc,h_base,xlo,ylo,zlo,soln,res);
    
}
