/**
 * \file   physics_residual_face_hang_full_structured.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief External solver function wrapper for calculating residual on 
 *        hang-full face on structured grid (left element is one level finer)
 * 
 * Created on September 5, 2018, 2:33 PM
 */

/* header files */
#include "physics_residual_face_hang_full_structured.h"


void physics_residual_face_hang_full_structured(int *dim,int *nfields,int *nelem_subgrid,
                            int *side,int *level_l,int *level_r,int *type_l,int *type_r,
                            int *ind_l,int *ind_r,int *dof,double *h_base,
                            double *xlo,double *ylo,double *zlo,
                            double *soln,double *rhs){
    
    /* external solver kernel */
    if(*dim==2){
        
        solver_residual_face_hang_full_structured_2d(
                nfields,nelem_subgrid,side,
                level_l,level_r,
                &type_l[0],&type_l[1],&type_r[0],
                &ind_l[0],&ind_l[1],&ind_r[0],
                dof,h_base,xlo,ylo,zlo,
                &soln[ind_l[0]],
                &soln[ind_l[1]],
                &soln[ind_r[0]],
                rhs);
        
    }else{
        
        solver_residual_face_hang_full_structured_3d(
                nfields,nelem_subgrid,side,
                level_l,level_r,
                &type_l[0],&type_l[1],&type_l[2],&type_l[3],&type_r[0],
                &ind_l[0],&ind_l[1],&ind_l[2],&ind_l[3],&ind_r[0],
                dof,h_base,xlo,ylo,zlo,
                &soln[ind_l[0]],
                &soln[ind_l[1]],
                &soln[ind_l[2]],
                &soln[ind_l[3]],
                &soln[ind_r[0]],
                rhs);
        
    }
    
}

