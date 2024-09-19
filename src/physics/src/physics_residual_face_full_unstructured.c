/**
 * \file   physics_residual_face_full_unstructured.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief External solver function wrapper for calculating residual at 
 *        AMR quad face with no hanging elements (two elements sharing full face)
 * 
 * Created on August 31, 2018, 9:44 AM
 */

/* header files */
#include "physics_residual_face_full_unstructured.h"


void physics_residual_face_full_unstructured(int *dim,int *nfields,int *nelem_subgrid,int *side_l,int *side_r,
                                             int *type_l,int *type_r,int *indl,int *indr,int *dof,
                                             double *geom_l,double *geom_r,
                                             double *ql,double *qr,double *res){

    /* external solver kernel */
    solver_residual_face_full_unstructured(dim,nfields,nelem_subgrid,side_l,side_r,type_l,type_r,indl,indr,dof,geom_l,geom_r,ql,qr,res);
    
}

