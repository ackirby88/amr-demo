/**
 * \file   physics_residual_face_full_structured.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief External solver function wrapper for calculating residual at 
 *        AMR quad face with no hanging elements (two elements sharing full face)
 * 
 * Created on August 28, 2018, 12:38 PM
 */

/* header files */
#include "physics_residual_face_full_structured.h"


void physics_residual_face_full_structured(int *dim,int *nfields,int *nelem_subgrid,int *side,int *level,
                                           int *type_l,int *type_r,int *indl,int *indr,int *dof,
                                           double *h_base,double *xlo,double *ylo,double *zlo,
                                           double *ql,double *qr,double *res){
    
    /* external solver kernel */
    solver_residual_face_full_structured(dim,nfields,nelem_subgrid,side,level,type_l,type_r,indl,indr,dof,h_base,xlo,ylo,zlo,ql,qr,res);
    
}
