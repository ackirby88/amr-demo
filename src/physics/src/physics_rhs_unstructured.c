/**
 * \file   physics_rhs_unstructured.c
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief Unstructured grid right-hand-side operator.
 *
 * This function demonstrates how to parse the AMR data structures 
 * (cell_info, face_info, edge_info, node_info) to calculate discretized physics 
 * using neighboring cell information on unstructured quad/hex grids.
 * 
 * Created on August 29, 2018, 3:54 PM
 */

/* header files */
#include "physics_rhs_unstructured.h"


void physics_rhs_unstructured(int *dim,int *nfields,int *nelem_subgrid,int *ncell,
                              int *cell_info,double *geom,
                              int *nface,int *face_info,int *dof,
                              double *soln,double *rhs){
    
    int i;
    int cell_id;
    int face_id;
    int cell_index;
    int face_index;
    
    int ncorners;
    int facetype;
    
    int side_l;
    int side_r[4];
    
    int cell;
    int cell_l;
    int cell_r[4];
    
    int type;
    int type_l;
    int type_r[4];
    
    int ind_l;
    int ind_r[4];
    
    int level;
    int bc;
    
    
    
    ncorners = 4;
    if(*dim==3){
        ncorners = 8;
    }
    ncorners *= 3; // multiply in 3 coordinates per node (2D also)
        
    
    
    /* zero the rhs */
    for(i = 0; i < *dof; i++){
        rhs[i] = 0.0;
    }
    
    
    /* face loop - compute face residual */
    for(face_id = 0; face_id < *nface; face_id++){
        
        /* this cell's information */
        face_index = d_physics_interface.extern_face_size*face_id;
        facetype    = face_info[face_index+0];
        side_l      = face_info[face_index+1];
        cell_l      = face_info[face_index+2];
        
        cell_index  = d_physics_interface.extern_cell_size*cell_l;
        type_l      = cell_info[cell_index+0];
        ind_l       = cell_info[cell_index+2];
        
        
        switch (facetype){
            
            case 1:
                /* mesh boundary --> cell is real */
                
                //TODO: assign bc type
                //TODO: Im just using the element type right now since I assigned in
                //TODO: the solver initialization
                bc = type_l;
                
                physics_residual_bc_unstructured(dim,nfields,nelem_subgrid,
                                                 &side_l,&type_l,&bc,
                                                 &geom[ncorners*cell_l],
                                                 &soln[ind_l],&rhs[ind_l]);
                
                break;
                
            case 2: 
                /* full-full --> check both if real or ghost */
                
                /* full neighbor's information */
                cell_r[0]   = face_info[face_index+3];
                side_r[0]   = face_info[face_index+7];
                type_r[0]   = cell_info[d_physics_interface.extern_cell_size*cell_r[0]+0];
                ind_r[0]    = cell_info[d_physics_interface.extern_cell_size*cell_r[0]+2];
                
                
                physics_residual_face_full_unstructured(
                    dim,nfields,nelem_subgrid,
                    &side_l,&side_r[0],
                    &type_l,&type_r[0],
                    &ind_l,&ind_r[0],dof,
                    &geom[ncorners*cell_l],
                    &geom[ncorners*cell_r[0]],
                    &soln[ind_l],&soln[ind_r[0]],
                    rhs);
                
                break;
                
            case 3:
                /* hanging 2D --> check all if real or ghost: ind > dof = ghost */
                
                /* hanging neighbors' information */
                cell_r[0]   = face_info[face_index+3];
                cell_r[1]   = face_info[face_index+4];
                side_r[0]   = face_info[face_index+7];
                side_r[1]   = face_info[face_index+8];
                type_r[0]   = cell_info[d_physics_interface.extern_cell_size*cell_r[0]+0];
                type_r[1]   = cell_info[d_physics_interface.extern_cell_size*cell_r[1]+0];
                ind_r[0]    = cell_info[d_physics_interface.extern_cell_size*cell_r[0]+2];
                ind_r[1]    = cell_info[d_physics_interface.extern_cell_size*cell_r[1]+2];
                
                physics_residual_face_hang_unstructured(
                    dim,nfields,nelem_subgrid,
                    &side_l,&side_r[0],
                    &type_l,&type_r[0],
                    &cell_l,&cell_r[0],
                    &ind_l,&ind_r[0],
                    dof,geom,soln,rhs);
                
                break;
                
            case 5:
                /* hanging 3D --> check all if real or ghost: ind > dof = ghost */
                
                /* hanging neighbors' information */
                cell_r[0]   = face_info[face_index+3];
                cell_r[1]   = face_info[face_index+4];
                cell_r[2]   = face_info[face_index+5];
                cell_r[3]   = face_info[face_index+6];
                side_r[0]   = face_info[face_index+7];
                side_r[1]   = face_info[face_index+8];
                side_r[2]   = face_info[face_index+9];
                side_r[3]   = face_info[face_index+10];
                type_r[0]   = cell_info[d_physics_interface.extern_cell_size*cell_r[0]+0];
                type_r[1]   = cell_info[d_physics_interface.extern_cell_size*cell_r[1]+0];
                type_r[2]   = cell_info[d_physics_interface.extern_cell_size*cell_r[2]+0];
                type_r[3]   = cell_info[d_physics_interface.extern_cell_size*cell_r[3]+0];
                ind_r[0]    = cell_info[d_physics_interface.extern_cell_size*cell_r[0]+2];
                ind_r[1]    = cell_info[d_physics_interface.extern_cell_size*cell_r[1]+2];
                ind_r[2]    = cell_info[d_physics_interface.extern_cell_size*cell_r[2]+2];
                ind_r[3]    = cell_info[d_physics_interface.extern_cell_size*cell_r[3]+2];
                
                physics_residual_face_hang_unstructured(
                    dim,nfields,nelem_subgrid,
                    &side_l,&side_r[0],
                    &type_l,&type_r[0],
                    &cell_l,&cell_r[0],
                    &ind_l,&ind_r[0],
                    dof,geom,soln,rhs);
                
                break;
                
        }
        
    }
    
    
    /* element loop: compute cell subgrid face, volume, source term residuals */
    for(cell_id = 0; cell_id < *ncell; cell_id++){
        
        cell_index = d_physics_interface.extern_cell_size*cell_id;
        type    = cell_info[cell_index+0];
        level   = cell_info[cell_index+1];
        cell    = cell_info[cell_index+2];
        
//TODO: physics_residual_source (dim,nfields,nelem_subgrid,&type_l,&geom[ncorners*cell_id],&soln[cell],&rhs[cell]);
        physics_residual_subgrid(dim,nfields,nelem_subgrid,&type,&level,&geom[ncorners*cell_id],&soln[cell],&rhs[cell]);
        
    }
    
    /* rhs = -R */
    for(i = 0; i < *dof; i++){
        rhs[i] *= -1.0;
    }
    
}

