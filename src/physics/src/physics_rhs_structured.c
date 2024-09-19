/**
 * \file   physics_rhs_structured.c
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief Structured grid right-hand-side operator.
 * 
 * This function demonstrates how to parse the AMR data structures 
 * (cell_info, face_info, edge_info, node_info) to calculate discretized physics 
 * using neighboring cell information on structured quad/hex grids.
 *
 * Created on August 29, 2018, 3:50 PM
 */

/* header files */
#include "physics_rhs_structured.h"


void physics_rhs_structured(int *dim,int *nfields,int *nelem_subgrid,int *ncell,
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
    
    int side;
    int side_plus;
    
    int cell;
    int ql[4];
    int qr[4];
    
    int type;
    int type_l[4];
    int type_r[4];
    
    int ind;
    int ind_l[4];
    int ind_r[4];
    
    int level;
    int level_l;
    int level_r;
    
    double h_base;
    double xlo[2] = {0.0,0.0};
    double ylo[2] = {0.0,0.0};
    double zlo[2] = {0.0,0.0};
    
    
    ncorners = 4;
    if(*dim==3){
        ncorners = 8;
    }
    ncorners *= 3; // multiple in 3 coordinates per node (2D also)
        
    
    
    /* zero the rhs */
    for(i = 0; i < *dof; i++){
        rhs[i] = 0.0;
    }
    
    
    /* face loop - compute face residual */
    for(face_id = 0; face_id < *nface; face_id++){
        
        /* this cell's information */
        face_index = d_physics_interface.extern_face_size*face_id;
        facetype= face_info[face_index+0];
        side    = face_info[face_index+1];
        cell    = face_info[face_index+2];
        
        cell_index = d_physics_interface.extern_cell_size*cell;
        type    = cell_info[cell_index+0];
        level   = cell_info[cell_index+1];
        ind     = cell_info[cell_index+2];
        
        
        if(side%2==0){
            /** lo face */
            // ql is cell_r
            // qr is cell_l
            
            ql[0] = face_info[face_index+3];
            ql[1] = face_info[face_index+4];
            ql[2] = face_info[face_index+5];
            ql[3] = face_info[face_index+6];
            
            type_l[0] = cell_info[d_physics_interface.extern_cell_size*ql[0]+0];
            type_l[1] = cell_info[d_physics_interface.extern_cell_size*ql[1]+0];
            type_l[2] = cell_info[d_physics_interface.extern_cell_size*ql[2]+0];
            type_l[3] = cell_info[d_physics_interface.extern_cell_size*ql[3]+0];
            
            ind_l[0] = cell_info[d_physics_interface.extern_cell_size*ql[0]+2];
            ind_l[1] = cell_info[d_physics_interface.extern_cell_size*ql[1]+2];
            ind_l[2] = cell_info[d_physics_interface.extern_cell_size*ql[2]+2];
            ind_l[3] = cell_info[d_physics_interface.extern_cell_size*ql[3]+2];
            
            qr[0]     = cell;
            type_r[0] = cell_info[d_physics_interface.extern_cell_size*qr[0]+0];
            ind_r[0]  = cell_info[d_physics_interface.extern_cell_size*qr[0]+2];
            
        }else{
            /** hi face */
            // ql is cell_l
            // qr is cell_r
            
            qr[0] = face_info[face_index+3];
            qr[1] = face_info[face_index+4];
            qr[2] = face_info[face_index+5];
            qr[3] = face_info[face_index+6];
            
            type_r[0] = cell_info[d_physics_interface.extern_cell_size*qr[0]+0];
            type_r[1] = cell_info[d_physics_interface.extern_cell_size*qr[1]+0];
            type_r[2] = cell_info[d_physics_interface.extern_cell_size*qr[2]+0];
            type_r[3] = cell_info[d_physics_interface.extern_cell_size*qr[3]+0];
            
            ind_r[0] = cell_info[d_physics_interface.extern_cell_size*qr[0]+2];
            ind_r[1] = cell_info[d_physics_interface.extern_cell_size*qr[1]+2];
            ind_r[2] = cell_info[d_physics_interface.extern_cell_size*qr[2]+2];
            ind_r[3] = cell_info[d_physics_interface.extern_cell_size*qr[3]+2];
            
            ql[0]     = cell;
            type_l[0] = cell_info[d_physics_interface.extern_cell_size*ql[0]+0];
            ind_l[0]  = cell_info[d_physics_interface.extern_cell_size*ql[0]+2];
            
        }
        
        level_l = cell_info[3*ql[0]+1];
        level_r = cell_info[3*qr[0]+1];
        
        switch (facetype){
            
            case 1:
                /* boundary case --> cell is real */
                
                xlo[0] = geom[ncorners*cell+0];
                ylo[0] = geom[ncorners*cell+1];
                zlo[0] = geom[ncorners*cell+2];
                h_base = d_physics_interface.coarse_grid_h;
                
                physics_residual_bc_structured(dim,nfields,nelem_subgrid,
                                               &side,&type,&level,&d_physics_fluid.bc[side],
                                               &h_base,&xlo[0],&ylo[0],&zlo[0],
                                               &soln[ind],&rhs[ind]);
                
                break;
                
            case 2: 
                /* full-full case --> check both if real or ghost: ind > dof = ghost */
                
                xlo[0] = geom[ncorners*ql[0]+0];
                ylo[0] = geom[ncorners*ql[0]+1];
                zlo[0] = geom[ncorners*ql[0]+2];
                xlo[1] = geom[ncorners*qr[0]+0];
                ylo[1] = geom[ncorners*qr[0]+1];
                zlo[1] = geom[ncorners*qr[0]+2];
                h_base = d_physics_interface.coarse_grid_h;
                
                
                side_plus = (side%2==0) ? side + 1 : side;
                physics_residual_face_full_structured(dim,nfields,nelem_subgrid,
                                        &side_plus,&level_l,&type_l[0],&type_r[0],
                                        &ind_l[0],&ind_r[0],dof,
                                        &h_base,&xlo[0],&ylo[0],&zlo[0],
                                        &soln[ind_l[0]],&soln[ind_r[0]],rhs);
                
                break;
                
            case 3:
                /* hanging 2D --> check all if real or ghost: ind > dof = ghost */
                
                if(side%2==0){
                    /* lo face */
                    //ql = hanging
                    //qr = full
                    
                    xlo[0]  = geom[ncorners*ql[1]+0];
                    ylo[0]  = geom[ncorners*ql[1]+1];
                    zlo[0]  = geom[ncorners*ql[1]+2];
                    xlo[1]  = geom[ncorners*qr[0]+0];
                    ylo[1]  = geom[ncorners*qr[0]+1];
                    zlo[1]  = geom[ncorners*qr[0]+2];
                    h_base = d_physics_interface.coarse_grid_h;
                    
                    side_plus = side + 1;
                    physics_residual_face_hang_full_structured(dim,nfields,nelem_subgrid,
                            &side_plus,&level_l,&level_r,type_l,type_r,
                            ind_l,ind_r,dof,&h_base,xlo,ylo,zlo,soln,rhs);
                    
                }else{
                    /* hi face */
                    //ql = full
                    //qr = hanging
                    
                    xlo[0]  = geom[ncorners*ql[0]+0];
                    ylo[0]  = geom[ncorners*ql[0]+1];
                    zlo[0]  = geom[ncorners*ql[0]+2];
                    xlo[1]  = geom[ncorners*qr[1]+0];
                    ylo[1]  = geom[ncorners*qr[1]+1];
                    zlo[1]  = geom[ncorners*qr[1]+2];
                    h_base = d_physics_interface.coarse_grid_h;
                    
                    side_plus = side;
                    physics_residual_face_full_hang_structured(dim,nfields,nelem_subgrid,
                            &side_plus,&level_l,&level_r,type_l,type_r,
                            ind_l,ind_r,dof,&h_base,xlo,ylo,zlo,soln,rhs);
                    
                }
                break;
                
            case 5:
                /* hanging --> check all if real or ghost: ind > dof = ghost */
                
                if(side%2==0){
                    /* lo face */
                    // ql = hanging
                    // qr = full
                    
                    xlo[0]  = geom[ncorners*ql[3]+0];
                    ylo[0]  = geom[ncorners*ql[3]+1];
                    zlo[0]  = geom[ncorners*ql[3]+2];
                    xlo[1]  = geom[ncorners*qr[0]+0];
                    ylo[1]  = geom[ncorners*qr[0]+1];
                    zlo[1]  = geom[ncorners*qr[0]+2];
                    
                    side_plus = side + 1;
                    physics_residual_face_hang_full_structured(dim,nfields,nelem_subgrid,
                            &side_plus,&level_l,&level_r,type_l,type_r,
                            ind_l,ind_r,dof,&h_base,xlo,ylo,zlo,soln,rhs);
                    
                }else{
                    /* hi face */
                    // ql = full
                    // qr = hanging
                    
                    xlo[0]  = geom[ncorners*ql[0]+0];
                    ylo[0]  = geom[ncorners*ql[0]+1];
                    zlo[0]  = geom[ncorners*ql[0]+2];
                    xlo[1]  = geom[ncorners*qr[3]+0];
                    ylo[1]  = geom[ncorners*qr[3]+1];
                    zlo[1]  = geom[ncorners*qr[3]+2];
                    
                    side_plus = side;
                    physics_residual_face_full_hang_structured(dim,nfields,nelem_subgrid,
                            &side_plus,&level_l,&level_r,type_l,type_r,
                            ind_l,ind_r,dof,&h_base,xlo,ylo,zlo,soln,rhs);
                            
                }
                break;
                
        }
        
    }
    
    
    /* element loop: compute cell subgrid face, volume, source term residuals */
    for(cell_id = 0; cell_id < *ncell; cell_id++){
        
        cell_index = d_physics_interface.extern_cell_size*cell_id;
        type    = cell_info[cell_index+0];
        level   = cell_info[cell_index+1];
        cell    = cell_info[cell_index+2];
        
//TODO: physics_residual_source (dim,nfields,nelem_subgrid,&type,&level,&geom[ncorners*cell_id],&soln[cell],&rhs[cell]);
        physics_residual_subgrid(dim,nfields,nelem_subgrid,&type,&level,&geom[ncorners*cell_id],&soln[cell],&rhs[cell]);
        
    }
    
    /* rhs = -R */
    for(i = 0; i < *dof; i++){
        rhs[i] *= -1.0;
    }
    
}

