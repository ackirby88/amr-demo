/**
 * \file    physics_time_step.c
 * \ingroup physics_group
 * \author  akirby
 * 
 * \brief External solver function wrapper for calculating the maximum stable 
 *        time step allowable for the cells on this mpi rank.
 *
 * Created on September 19, 2018, 11:10 AM
 */

/* header files */
#include "physics_time_step.h"


void physics_time_step(int *dim,int *nfields,int *nelem_subgrid,
                       int *ncell,int *cell_info,double *geom,double *soln, 
                       double *dt){

    int cell_id;
    int ncorners;
    double cfl;
    
    int cell_index;
    int type;
    int level;
    int cell;
    
    
    ncorners = 4;
    if(*dim==3){
        ncorners = 8;
    }
    ncorners *= 3; // multiply in 3 coordinates per node (2D also)
    
    cfl = d_physics_fluid.cfl;
    
    for(cell_id = 0; cell_id < *ncell; cell_id++){
        
        cell_index = d_physics_interface.extern_cell_size*cell_id;
        type    = cell_info[cell_index+0];
        level   = cell_info[cell_index+1];
        cell    = cell_info[cell_index+2];
        
        solver_time_step(dim,nfields,nelem_subgrid,&type,&level,&geom[ncorners*cell_id],&soln[cell],&cfl,dt);
        
    }
    
}

