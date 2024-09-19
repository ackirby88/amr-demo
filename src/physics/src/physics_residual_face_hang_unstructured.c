/**
 * \file   physics_residual_face_hang_unstructured.c
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief External solver function wrapper for calculating residual on 
 *        hanging face on unstructured grid
 *
 * Created on September 5, 2018, 3:51 PM
 */


/* header files */
#include "physics_residual_face_hang_unstructured.h"


void physics_residual_face_hang_unstructured(int *dim,int *nfields,int *nelem_subgrid,
                                             int *side_l,int *side_r,
                                             int *type_l,int *type_r,
                                             int *cell_l,int *cell_r,
                                             int *ind_l,int *ind_r,
                                             int *dof,double *geom,
                                             double *soln,double *rhs){    
    
    /* external solver kernel */
    if(*dim==2){
        
        solver_residual_face_hang_unstructured_2d(
                nfields,nelem_subgrid,
                &side_l[0],&side_r[0],&side_r[1],
                &type_l[0],&type_r[0],&type_r[1],
                &ind_l[0],
                &ind_r[0],
                &ind_r[1],
                dof,
                &geom[3*4*cell_l[0]],
                &geom[3*4*cell_r[0]],
                &geom[3*4*cell_r[1]],
                &soln[ind_l[0]],
                &soln[ind_r[0]],
                &soln[ind_r[1]],
                rhs);
        
    }else{
        
        solver_residual_face_hang_unstructured_3d(
                nfields,nelem_subgrid,
                &side_l[0],&side_r[0],&side_r[1],&side_r[2],&side_r[3],
                &type_l[0],&type_r[0],&type_r[1],&type_r[2],&type_r[3],
                &ind_l[0],
                &ind_r[0],
                &ind_r[1],
                &ind_r[2],
                &ind_r[3],
                dof,
                &geom[3*8*cell_l[0]],
                &geom[3*8*cell_r[0]],
                &geom[3*8*cell_r[1]],
                &geom[3*8*cell_r[2]],
                &geom[3*8*cell_r[3]],
                &soln[ind_l[0]],
                &soln[ind_r[0]],
                &soln[ind_r[1]],
                &soln[ind_r[2]],
                &soln[ind_r[3]],
                rhs);
        
    }
    
    
}
