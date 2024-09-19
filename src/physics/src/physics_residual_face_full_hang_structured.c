/**
 * \file   physics_residual_face_full_hang_structured.c
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief External solver function wrapper for calculating residual on 
 *        full-hang face on structured grid (right element is one level finer)
 *
 * Created on September 5, 2018, 1:36 PM
 */

/* header files */
#include "physics_residual_face_full_hang_structured.h"


void physics_residual_face_full_hang_structured(int *dim,int *nfields,int *nelem_subgrid,
                            int *side,int *level_l,int *level_r,int *type_l,int *type_r,
                            int *ind_l,int *ind_r,int *dof,double *h_base,
                            double *xlo,double *ylo,double *zlo,
                            double *soln,double *rhs){
    
    /* external solver kernel */
    if(*dim==2){
        
        solver_residual_face_full_hang_structured_2d(
                nfields,nelem_subgrid,side,
                level_l,level_r,
                &type_l[0],&type_r[0],&type_r[1],
                &ind_l[0],&ind_r[0],&ind_r[1],
                dof,h_base,xlo,ylo,zlo,
                &soln[ind_l[0]],
                &soln[ind_r[0]],
                &soln[ind_r[1]],
                rhs);
        
    }else{
        
        solver_residual_face_full_hang_structured_3d(
                nfields,nelem_subgrid,side,
                level_l,level_r,
                &type_l[0],&type_r[0],&type_r[1],&type_r[2],&type_r[3],
                &ind_l[0],&ind_r[0],&ind_r[1],&ind_r[2],&ind_r[3],
                dof,h_base,xlo,ylo,zlo,
                &soln[ind_l[0]],
                &soln[ind_r[0]],
                &soln[ind_r[1]],
                &soln[ind_r[2]],
                &soln[ind_r[3]],
                rhs);
        
    }
    
}
