/**
 * \file   amr_vtk.c
 * \ingroup amr_group
 * \author akirby
 *
 * \brief Visualization functions for the AMR code module.
 * 
 * Created on August 20, 2018, 5:25 PM
 */


/* header files */
#include "amr_vtk.h"


void vtk_visualize_grid(p4est_t *p4est){
    
    p4est_vtk_write_file(p4est,NULL,P4EST_STRING "_grid");
    
}


void vtk_write_solution(p4est_t *p4est,int timestep){
    
    char filename[BUFF_SIZE] = { '\0' };
    p4est_locidx_t numquads;
    snprintf(filename,13,"soln_%07d",timestep);
    
    numquads = p4est->local_num_quadrants;
    
    
    vtk_plot_fluid_t fluid_data;
    fluid_data.rho  = sc_array_new_size(sizeof(double),numquads*P4EST_CHILDREN);
    fluid_data.rhou = sc_array_new_size(sizeof(double),numquads*P4EST_CHILDREN);
    fluid_data.rhov = sc_array_new_size(sizeof(double),numquads*P4EST_CHILDREN);
    fluid_data.rhoe = sc_array_new_size(sizeof(double),numquads*P4EST_CHILDREN);
#ifdef P4_TO_P8
    fluid_data.rhow = sc_array_new_size(sizeof(double),numquads*P4EST_CHILDREN);
#endif
    
    
    
    p4est_iterate(p4est,NULL,(void *) &fluid_data,interpolate_solution,
                 NULL,
#ifdef P4_TO_P8
                 NULL,
#endif
                 NULL);
    

    
    /* create VTK output context and set its parameters */
    p4est_vtk_context_t *context = p4est_vtk_context_new(p4est,filename);
    p4est_vtk_context_set_scale(context,0.99999);

    
    /* begin writing the output files */
    context = p4est_vtk_write_header(context);
    SC_CHECK_ABORT(context != NULL,P4EST_STRING "_vtk: Error writing vtk header");
    

    context = p4est_vtk_write_cell_dataf(   context, 
                                            1,  /* write tree */
                                            1,  /* write level */
                                            1,  /* write mpi rank */
                                            0,  /* do not wrap the mpi rank */
                                            0,  /* num custom cell scalar data. */
                                            0,  /* num custom cell vector data. */
                                            context); /* mark the end of the variable cell data. */
    
    SC_CHECK_ABORT(context != NULL,P4EST_STRING "_vtk: Error writing cell data");

    
    
    /* write one scalar field: the solution value */
    context = p4est_vtk_write_point_dataf(  context, 
                                            NFIELDS,/* number of scalar fields  */
                                            0,      /* number of vector fields  */
                                            "rho", fluid_data.rho,
                                            "rhou",fluid_data.rhou,
                                            "rhov",fluid_data.rhov,
                                            "rhoe",fluid_data.rhoe,
#ifdef P4_TO_P8
                                            "rhow",fluid_data.rhow,
#endif
                                            context);
    
    
    SC_CHECK_ABORT(context != NULL,P4EST_STRING "_vtk: Error writing cell data");
    
    
    const int retval = p4est_vtk_write_footer (context);
    SC_CHECK_ABORT(!retval, P4EST_STRING "_vtk: Error writing footer");
    
    
    
    sc_array_destroy(fluid_data.rho);
    sc_array_destroy(fluid_data.rhou);
    sc_array_destroy(fluid_data.rhov);
    sc_array_destroy(fluid_data.rhoe);
#ifdef P4_TO_P8
    sc_array_destroy(fluid_data.rhow);
#endif
    
  
}


void interpolate_solution(p4est_iter_volume_info_t *info, void *user_data){
    
    vtk_plot_fluid_t   *soln_viz = (vtk_plot_fluid_t *) user_data;
    p4est_t            *p4est = info->p4est;
    //ctx_t              *d_ctx = (ctx_t *) p4est->user_pointer;
    p4est_quadrant_t   *q = info->quad;
    p4est_topidx_t      which_tree = info->treeid;
    p4est_locidx_t      local_id = info->quadid;
    p4est_tree_t       *tree;
    quad_data_t        *data = (quad_data_t *) q->p.user_data;
    p4est_locidx_t      arrayoffset;
    double             *this_u_ptr;
    
    
    int ind,ngp;
    int i,j;
    
    double oneO_total;
    double u_ave[NFIELDS] = {0.0};
    
    //double h;
    //h = d_ctx->d_grid.scale * (double) P4EST_QUADRANT_LEN(q->level)/(double) P4EST_ROOT_LEN;
    
    
    tree = p4est_tree_array_index(p4est->trees, which_tree);
    
    /* now the id is relative to the MPI process */
    local_id += tree->quadrants_offset;
    
    /* each local quadrant has 2^d (P4EST_CHILDREN) values in soln_viz */
    arrayoffset = P4EST_CHILDREN * local_id;
    
    
    //TODO: call external callback function to interpolate the data to corners
    
    // Average the subgrid fields
#ifdef P4_TO_P8
    int k;
    for(k=0;k<NXPATCH;k++){
        for(j=0;j<NXPATCH;j++){
            for(i=0;i<NXPATCH;i++){
                ind = NFIELDS*NXPATCH*NXPATCH*k + NFIELDS*NXPATCH*j + NFIELDS*i;
                u_ave[0] += data->soln[ind+0];
                u_ave[1] += data->soln[ind+1];
                u_ave[2] += data->soln[ind+2];
                u_ave[3] += data->soln[ind+3];
                u_ave[4] += data->soln[ind+4];
            }
        }
    }
    oneO_total = 1.0 / (double) (NXPATCH*NXPATCH*NXPATCH);
    u_ave[0] *= oneO_total;
    u_ave[1] *= oneO_total;
    u_ave[2] *= oneO_total;
    u_ave[3] *= oneO_total;
    u_ave[4] *= oneO_total;
#else
    for(j=0;j<NXPATCH;j++){
        for(i=0;i<NXPATCH;i++){
            ind = NFIELDS*NXPATCH*j + NFIELDS*i;
            u_ave[0] += data->soln[ind+0];
            u_ave[1] += data->soln[ind+1];
            u_ave[2] += data->soln[ind+2];
            u_ave[3] += data->soln[ind+3];
        }
    }
    oneO_total = 1.0 / (double) (NXPATCH*NXPATCH);
    u_ave[0] *= oneO_total;
    u_ave[1] *= oneO_total;
    u_ave[2] *= oneO_total;
    u_ave[3] *= oneO_total;
#endif
    
  
    ngp = 2;
#ifdef P4_TO_P8
    for(k=0;k<2;k++){
        for(j=0;j<2;j++){
            for(i=0;i<2;i++){
                this_u_ptr = (double *) sc_array_index(soln_viz->rho ,arrayoffset+ngp*ngp*k+ngp*j+i);
                this_u_ptr[0] = u_ave[0];
                this_u_ptr = (double *) sc_array_index(soln_viz->rhou,arrayoffset+ngp*ngp*k+ngp*j+i);
                this_u_ptr[0] = u_ave[1];
                this_u_ptr = (double *) sc_array_index(soln_viz->rhov,arrayoffset+ngp*ngp*k+ngp*j+i);
                this_u_ptr[0] = u_ave[2];
                this_u_ptr = (double *) sc_array_index(soln_viz->rhow,arrayoffset+ngp*ngp*k+ngp*j+i);
                this_u_ptr[0] = u_ave[3];
                this_u_ptr = (double *) sc_array_index(soln_viz->rhoe,arrayoffset+ngp*ngp*k+ngp*j+i);
                this_u_ptr[0] = u_ave[4];
            }
        }
    }
#else
    for(j=0;j<2;j++){
        for(i=0;i<2;i++){
            this_u_ptr = (double *) sc_array_index(soln_viz->rho ,arrayoffset+ngp*j+i);
            this_u_ptr[0] = u_ave[0];
            this_u_ptr = (double *) sc_array_index(soln_viz->rhou,arrayoffset+ngp*j+i);
            this_u_ptr[0] = u_ave[1];
            this_u_ptr = (double *) sc_array_index(soln_viz->rhov,arrayoffset+ngp*j+i);
            this_u_ptr[0] = u_ave[2];
            this_u_ptr = (double *) sc_array_index(soln_viz->rhoe,arrayoffset+ngp*j+i);
            this_u_ptr[0] = u_ave[3];
        }
    }
#endif
    
    
}
