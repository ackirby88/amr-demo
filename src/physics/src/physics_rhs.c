/**
 * \file   physics_rhs.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Physics right-hand-side (rhs) operator function interface.
 * 
 * Physics right-hand-side (rhs) operator interface which calls structured 
 * or unstructured rhs function.
 * 
 * Created on August 24, 2018, 3:50 PM
 */

/* header files */
#include "physics_rhs.h"


void physics_rhs(int *dim,int *unstructured,int *nfields,int *nelem_subgrid,
                 int *ncell,int *cell_info,double *geom,int *nface,
                 int *face_info,int *dof,double *soln,double *rhs){
    
    if(*unstructured){
        physics_rhs_unstructured(dim,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,soln,rhs);
    }else{
        physics_rhs_structured  (dim,nfields,nelem_subgrid,ncell,cell_info,geom,nface,face_info,dof,soln,rhs);
    }
    
}
