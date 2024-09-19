/**
 * \file   amr_regrid.c
 * \ingroup amr_group
 * \author akirby
 * 
 * \brief Regridding functions for the AMR code module.
 *
 * Created on September 4, 2018, 12:29 PM
 */

/* header files */
#include "amr_regrid.h"


void regrid_solution(ctx_t *d_ctx,p4est_t *p4est,int initial,
                     int point_regrid,int feature_regrid,int allow_coarsening){
    
    double total_t1,total_t2;
    double t1,t2;
    
    double max_time;
    double max_point_time = 0.0;
    double max_box_time = 0.0;
    double max_feature_time = 0.0;
    double coarsen_time = 0.0;
    double balance_time = 0.0;
    double partition_time = 0.0;
    
    
    if(d_ctx->d_mpi.rank==0) printf("\n[regrid] Into regridding...\n");
    
    
    utilities_timer(&total_t1);   /* start timer */
        tag_reset(p4est,no_tag);
        
        /* regrid points */
        if(point_regrid){
            utilities_timer(&t1);
                regrid_points(d_ctx,p4est,initial);
            utilities_timer(&t2);
            max_point_time = utilities_mpireducemax_double(d_ctx->d_mpi.comm,t2-t1);
        }
        
        /* regrid box region */
        if(d_ctx->d_grid.refined_box){
            utilities_timer(&t1);
                regrid_box(d_ctx,p4est,initial);
            utilities_timer(&t2);
            max_box_time = utilities_mpireducemax_double(d_ctx->d_mpi.comm,t2-t1);
        }
        
        /* regrid features */
        if(feature_regrid){
            utilities_timer(&t1);
                regrid_feature(d_ctx,p4est,initial);
            utilities_timer(&t2);
            max_feature_time = utilities_mpireducemax_double(d_ctx->d_mpi.comm,t2-t1);
        }
        
        /* coarsen grid */
        if(allow_coarsening){
            regrid_coarsen(d_ctx,p4est,initial,&coarsen_time);
        }
        
        /* balance and partition */
        regrid_balance_partition(d_ctx,p4est,initial,&balance_time,&partition_time);
        
    utilities_timer(&total_t2);   /* stop timer*/
    
    
    max_time = utilities_mpireducemax_double(d_ctx->d_mpi.comm,total_t2-total_t1);
    
    p4est_utilities_mesh_stats(p4est);
    utilities_write_amr_regrid(d_ctx->d_mpi.rank,max_time,
                               balance_time,partition_time,
                               max_point_time,max_box_time,
                               max_feature_time);
    
    if(d_ctx->d_mpi.rank==0) printf("[regrid] Done regridding...\n\n");
    
}


void regrid_coarsen(ctx_t *d_ctx,p4est_t *p4est,int initial,double *coarsen_time){
    
    int recursive = 0;
    int callbackorphans = 0;
    double t1,t2;
    
    /* coarsen if not tagged*/
    utilities_timer(&t1);
        if(initial==0){
            p4est_coarsen_ext(p4est,recursive,callbackorphans,tag_coarsen,NULL,regrid_replace_quads);
        }
    utilities_timer(&t2);
    *coarsen_time = utilities_mpireducemax_double(d_ctx->d_mpi.comm,t2-t1);
    
}


void regrid_balance_partition(ctx_t *d_ctx,p4est_t *p4est,int initial,double *bal_time,double *part_time){
    
    int partition_for_coarsening = 1;
    double t1,t2;
    
    
    /* balance mesh 2:1 element ratio */
    utilities_timer(&t1);
        p4est_balance_ext(p4est,P4EST_CONNECT_FULL,NULL,regrid_replace_quads);
        /* Do not use initialize quadrant data -- 
         *   The coarsen operator destroys the 2:1 balance cells because they 
         *   are not tagged. Then we rebalance, it uses regrid_replace_quads.
         *   Therefore the grid changes even if we don't change the solution.
         */
        
    utilities_timer(&t2);
    *bal_time = utilities_mpireducemax_double(d_ctx->d_mpi.comm,t2-t1);
    
    /* partition mesh based on weighting */
    utilities_timer(&t1);
        p4est_partition(p4est,partition_for_coarsening,regrid_partition_weight);
    utilities_timer(&t2);
    *part_time = utilities_mpireducemax_double(d_ctx->d_mpi.comm,t2-t1);
    
}


void regrid_points(ctx_t *d_ctx,p4est_t *p4est,int initial){
    
    int recursive;
    int allowcoarsening;
    int i,j;
    
    double *p;
    sc_array_t *point_data;
    
    external_regrid_points(&d_ctx->d_grid.construct_grid,
                           &d_ctx->d_grid.nregrid_pts,
                           &d_ctx->d_grid.regrid_width,
                            d_ctx->d_grid.regrid_xyzh);
    
    
    if(d_ctx->d_grid.nregrid_pts>0){
        point_data = sc_array_new_size(4*sizeof(double),d_ctx->d_grid.nregrid_pts);
        
        for(j = 0; j < d_ctx->d_grid.nregrid_pts; j++){
            p = sc_array_index(point_data,j);

            for(i = 0; i < 4; i++){
                p[i] = d_ctx->d_grid.regrid_xyzh[4*j+i];
            }
        }
        
        if(d_ctx->d_grid.construct_grid){
            
            /* build one level at a time */
            recursive = 0;
            allowcoarsening = 0;
            for(i = d_ctx->d_grid.min_level; i <= d_ctx->d_grid.max_level; i++){
                
                /* reset regrid tags */
                tag_reset(p4est,no_tag);
                
                /* search the point list */
                p4est_search(p4est,NULL,tag_point_search,point_data);
                
                /* refine tagged quadrants */
                p4est_refine_ext(p4est,recursive,i,tag_point,initialize_quadrant_data,NULL);
                
                /* partition cells for load balancing as the grid is built */
                p4est_partition(p4est,allowcoarsening,regrid_partition_weight);

            }
            d_ctx->d_grid.construct_grid = 0;
            
        }else{
                
            /* search the point list */
            p4est_search(p4est,NULL,tag_point_search,point_data);
            
            /* refine tagged quadrants */
            recursive = 1;
            p4est_refine_ext(p4est,recursive,d_ctx->d_grid.max_level,tag_point,NULL,regrid_replace_quads);
            
        }
        sc_array_destroy(point_data);
        
    }
    
}


void regrid_box(ctx_t *d_ctx,p4est_t *p4est,int initial){
    
    int l;
    int recursive;
    int allowcoarsening;
    
    if(initial){
    
        /* build one level at a time */
        recursive = 0;
        allowcoarsening = 1;
        
        for(l = d_ctx->d_grid.min_level; l <= d_ctx->d_grid.max_level; l++){
            tag_reset(p4est,no_tag);
            p4est_refine_ext(p4est,recursive,d_ctx->d_grid.max_level,tag_box,initialize_quadrant_data,NULL);
            p4est_partition(p4est,allowcoarsening,regrid_partition_weight);
        }
        
    }else{
        
        recursive = 1;
        p4est_refine_ext(p4est,recursive,d_ctx->d_grid.max_level,tag_box,NULL,regrid_replace_quads);
        
    }
    
}


void regrid_feature(ctx_t *d_ctx,p4est_t *p4est,int initial){
    
    int l;
    int recursive;
    int allowcoarsening;
    
    if(initial){
    
        /* build one level at a time */
        recursive = 1;
        allowcoarsening = 0;
        
        for(l = d_ctx->d_grid.min_level; l <= d_ctx->d_grid.max_level; l++){
            tag_reset(p4est,no_tag);
            p4est_refine_ext(p4est,recursive,d_ctx->d_grid.max_level,tag_feature,initialize_quadrant_data,NULL);
            p4est_partition(p4est,allowcoarsening,regrid_partition_weight);
        }
        
    }else{
        
        recursive = 0;
        p4est_refine_ext(p4est,recursive,d_ctx->d_grid.max_level,tag_feature,NULL,regrid_replace_quads);
        
    }
    
}


void regrid_replace_quads(p4est_t *p4est,
                          p4est_topidx_t which_tree,
                          int num_outgoing,
                          p4est_quadrant_t *outgoing[],
                          int num_incoming,
                          p4est_quadrant_t *incoming[]){
    
    
    double h;
    double xyz[3*P4EST_CHILDREN];
    double *uc[P4EST_CHILDREN];
    ctx_t  *d_ctx = (ctx_t *) p4est->user_pointer;
    quad_data_t *parent_data;
    quad_data_t *child_data;
    int display_quad = 0;
    
    int i;
    
    
    if (num_outgoing > 1) {
        
        /* coarsen */
        parent_data = (quad_data_t *) incoming[0]->p.user_data;        
        //TODO: set the parent type based on the children
        //TODO: for now, its sets to the last child's type
        
        for (i = 0; i < P4EST_CHILDREN; i++) {
          child_data = (quad_data_t *) outgoing[i]->p.user_data;
          parent_data->type = child_data->type;
          uc[i] = child_data->soln;
        }
        
        /* element size */
        h = d_ctx->d_grid.scale * P4EST_QUADRANT_LEN(incoming[0]->level) / P4EST_ROOT_LEN;
        p4est_utilities_quad_coordinates(p4est,which_tree,incoming[0],xyz,display_quad);
        external_coarsen_operator(d_ctx->d_grid.dim,d_ctx->d_grid.unstructured_flag,NFIELDS,NXPATCH,&parent_data->type,h,xyz,uc,parent_data->soln);
        
        
    } else if(num_incoming > 1) {
        
        /* refine */
        parent_data = (quad_data_t *) outgoing[0]->p.user_data;
        
        for(i=0;i<P4EST_CHILDREN;i++){
          child_data = (quad_data_t *) incoming[i]->p.user_data;
          child_data->type = parent_data->type;
          uc[i] = child_data->soln;
        }
        
        /* element size */
        h = d_ctx->d_grid.scale * P4EST_QUADRANT_LEN(outgoing[0]->level) / P4EST_ROOT_LEN;
        p4est_utilities_quad_coordinates(p4est,which_tree,outgoing[0],xyz,display_quad);
        external_refine_operator(d_ctx->d_grid.dim,d_ctx->d_grid.unstructured_flag,NFIELDS,NXPATCH,&parent_data->type,h,xyz,parent_data->soln,uc);

    }

}


int regrid_partition_weight(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant){

    //TODO: Callback function with cell type for partition weighting
    //quad_data_t *data;
    //int cell_type;
    //
    //data = (quad_data_t *) quadrant->p.user_data;
    //cell_type = data->type;
    
    return 1;

}
