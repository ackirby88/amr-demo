/**
 * \file   amr_tag.c
 * \ingroup amr_group
 * \author akirby
 *
 * \brief Tagging functions for regridding for the AMR code module.
 * 
 * Created on September 7, 2018, 9:15 AM
 */

/* header files */
#include "amr_tag.h"


int tag_reset(p4est_t *p4est,int tag_default){
    
    p4est_topidx_t t;
    p4est_topidx_t first_local_tree = p4est->first_local_tree;
    p4est_topidx_t last_local_tree = p4est->last_local_tree;
    sc_array_t *trees = p4est->trees;
    sc_array_t *quadrants;
    p4est_tree_t *tree;
    
    p4est_quadrant_t *quad;
    quad_data_t *data;
    
    size_t n_quads;
    size_t si;
    

    for(t = first_local_tree; t <= last_local_tree; t++){
        tree = p4est_tree_array_index(trees,t);
        quadrants = &(tree->quadrants);
        n_quads = quadrants->elem_count;
        
        for(si = 0; si < n_quads; si++){
            quad = p4est_quadrant_array_index(quadrants,si);
            data = (quad_data_t *) quad->p.user_data;
            data->tag = tag_default;
        }
        
    }
    return 0;    
}


int tag_all(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant){
    
    quad_data_t *data;
    data = (quad_data_t *) quadrant->p.user_data;
    data->tag = all_tag;
            
    return 1;
    
}


int tag_point(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant){
    
    quad_data_t *data = (quad_data_t *) quadrant->p.user_data;
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;
    
    
    if(data->tag == no_tag)return 0;
    if(quadrant->level == d_ctx->d_grid.max_level) return 0;
    if(data->tag == point_tag) return 1;
    
    return 0;
    
}


int tag_point_search(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant,p4est_locidx_t local_num,void *points_in){
  
    
    quad_data_t *data = (quad_data_t *) quadrant->p.user_data;
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;
    double *point = (double *) points_in;
    
    double xlo[3];
    double dx;
    double h;
    
  
    
    /* quadrant mesh size */
    h = d_ctx->d_grid.scale * (double) P4EST_QUADRANT_LEN (quadrant->level) / (double) P4EST_ROOT_LEN;

    /* lower corner of quadrant */
    p4est_utilities_get_xlo(p4est,which_tree,quadrant,xlo);

    /* point mesh size */
    dx = d_ctx->d_grid.regrid_width*point[3];

    /* check intersection with box around each point */   
    if(xlo[0] <= point[0]+dx && xlo[0]+h >= point[0]-dx &&
       xlo[1] <= point[1]+dx && xlo[1]+h >= point[1]-dx 
#ifdef P4_TO_P8
                                                        &&
       xlo[2] <= point[2]+dx && xlo[2]+h >= point[2]-dx 
#endif
                                                        ){
        
        /* quadrant contains point but not a leaf */
        if(local_num == -1) return 1;
        
        data->tag = point_tag;
        return 1;
        
    }
    
    return 0;

}



int tag_box(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant){
    
    double xlo[3];
    double h;
    
    /* access our simulation context data */
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;
    
    /* access our per quadrant data pointer */
    quad_data_t *data = (quad_data_t *) quadrant->p.user_data;
    
    /* check if already tagged */
    if(data->tag > 0) return 0;
    
    
    
    /* quadrant mesh size */
    h = d_ctx->d_grid.scale * (double) P4EST_QUADRANT_LEN (quadrant->level) / (double) P4EST_ROOT_LEN;
    
    /* lower corner of quadrant */
    p4est_utilities_get_xlo(p4est,which_tree,quadrant,xlo);
    
    
    
    if(xlo[0]+0.001*h >= d_ctx->d_grid.box_xlo[0] && 
       xlo[0]+0.999*h <= d_ctx->d_grid.box_xhi[0] &&
       xlo[1]+0.001*h >= d_ctx->d_grid.box_xlo[1] &&
       xlo[1]+0.999*h <= d_ctx->d_grid.box_xhi[1]
#ifdef P4_TO_P8
                                                  &&
       xlo[2]+0.001*h >= d_ctx->d_grid.box_xlo[2] && 
       xlo[2]+0.999*h <= d_ctx->d_grid.box_xhi[2] 
#endif
    ){
        data->tag = box_tag;
        return 1;
    }
    
    
    /* check overlapping */
    if(check_overlap(xlo[0],xlo[0],d_ctx->d_grid.box_xlo[0],d_ctx->d_grid.box_xhi[0]) &&
       check_overlap(xlo[1],xlo[1],d_ctx->d_grid.box_xlo[1],d_ctx->d_grid.box_xhi[1])     
#ifdef P4_TO_P8
                                                                                      &&
       check_overlap(xlo[2],xlo[2],d_ctx->d_grid.box_xlo[2],d_ctx->d_grid.box_xhi[2])
#endif
                                                                                      ){
        data->tag = box_tag;
        return 1;
    }
    
    /* otherwise */
    return 0;
    
}


int tag_feature(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant){
    
    double h;
    double xyz[3*P4EST_CHILDREN];
    int display_quad = 0;
    int tag;
    
    /* access our simulation context data */
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;
    
    /* access our per quadrant data pointer */
    quad_data_t *data = (quad_data_t *) quadrant->p.user_data;
    
    /* check if already tagged */
    if(data->tag > 0) return 0;
    
    /* element size */
    h = d_ctx->d_grid.scale * P4EST_QUADRANT_LEN(quadrant->level) / P4EST_ROOT_LEN;
    p4est_utilities_quad_coordinates(p4est,which_tree,quadrant,xyz,display_quad);
    
    /* initialize the solution data in the quadrant via callback function */
    external_tag_feature(d_ctx->d_grid.dim,d_ctx->d_grid.unstructured_flag,NFIELDS,NXPATCH,&data->type,h,xyz,data->soln,&tag);
    if(tag) data->tag = feature_tag;
    return tag;
    
}

int tag_coarsen(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *children[]){
    
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;
    quad_data_t *data;
    int i;
  
    /* loop over the children and if any children igbp tagged do not coarsen */
    for (i = 0; i < P4EST_CHILDREN; i++) {
        if(children[i]){
            data = (quad_data_t *) children[i]->p.user_data;
            if(data->tag > 0) return 0;
            if(children[i]->level == d_ctx->d_grid.min_level) return 0;
        }
    }
  
    /* all children are tagged 0 */
    return 1;
  
}


int check_overlap(double x1,double x2,double xx1,double xx2){
 
    if(x1 >= xx1 && x1 <= xx2) return 1;
    if(x2 >= xx1 && x2 <= xx2) return 1;
    if(xx1 >= x1 && xx1 <= x2) return 1;
    if(xx2 >= x1 && xx2 <= x2) return 1;
    
    return 0;
}
