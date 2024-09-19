/**
 * \file   amr_p4est_utilities.c
 * \ingroup amr_group
 * \author akirby
 *
 * \brief Utility functions for accessing info about p4est data.
 * 
 * Created on August 17, 2018, 2:45 PM
 */

/* header files */
#include "amr_p4est_utilities.h"


void p4est_utilities_quad_coordinates(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *q, double *xyz, int display_quad){
    
    
    p4est_topidx_t vt[P4EST_CHILDREN];
    p4est_connectivity_t *connectivity;
    const double *v;
    const p4est_topidx_t *tree_to_vertex;
    const double intsize = 1.0 / (double) P4EST_ROOT_LEN;
    
    int j,k;
    int xi,yi;
    double h2,eta_x,eta_y,eta_z;
    
    
    
    connectivity = p4est->connectivity;
    v = connectivity->vertices;
    tree_to_vertex = connectivity->tree_to_vertex;
    
    
    /* fill in element size */
    h2 = 0.5 * intsize * P4EST_QUADRANT_LEN(q->level);;
    
    
    for(k = 0; k < P4EST_CHILDREN; ++k) {
      vt[k] = tree_to_vertex[which_tree * P4EST_CHILDREN + k];
    }
    
    
    eta_z = 0.0;
    /* fill in quadrant node coordinates */
    k = 0;
#ifdef P4_TO_P8
    int zi;
    for(zi = 0; zi < 2; ++zi){
        eta_z = intsize * q->z + h2 * (1.0 + (2*zi - 1));
#endif
        for(yi = 0; yi < 2; ++yi){
            eta_y = intsize * q->y + h2 * (1.0 + (2*yi - 1));      
            for(xi = 0; xi < 2; ++xi){
                eta_x = intsize * q->x + h2 * (1.0 + (2*xi - 1));
                
                for(j = 0; j < 3; ++j){
                    xyz[3*k+j] =
                        ((1. - eta_z) * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[0] + j] +
                                                               eta_x  * v[3 * vt[1] + j]) +
                                               eta_y  * ((1. - eta_x) * v[3 * vt[2] + j] +
                                                               eta_x  * v[3 * vt[3] + j]))
#ifdef P4_TO_P8
                         +     eta_z  * ((1. - eta_y) * ((1. - eta_x) * v[3 * vt[4] + j] +
                                                               eta_x  * v[3 * vt[5] + j]) +
                                               eta_y  * ((1. - eta_x) * v[3 * vt[6] + j] +
                                                               eta_x  * v[3 * vt[7] + j]))
#endif
                        );
                }
                if(display_quad){
                    printf("Quad geometry: %f %f %f\n",xyz[3*k+0],xyz[3*k+1],xyz[3*k+2]);
                }
                k++;
            }
        }
#ifdef P4_TO_P8
    }
#endif
    
}


void p4est_utilities_get_xlo(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *q,double xyz[3]){
  
    p4est_qcoord_to_vertex (p4est->connectivity,which_tree,q->x,q->y,
#ifdef P4_TO_P8
        q->z,
#endif
        xyz);

}


p4est_locidx_t p4est_utilities_get_local_id_volume(p4est_iter_volume_info_t *info){

    p4est_t            *p4est;
    p4est_topidx_t     which_tree;
    p4est_locidx_t     local_id;
    p4est_tree_t       *tree;

    p4est = info->p4est;
    which_tree = info->treeid;
    tree = p4est_tree_array_index (p4est->trees, which_tree);
    local_id = info->quadid + tree->quadrants_offset; 

    return local_id;
}


p4est_locidx_t p4est_utilities_get_local_id_face_hang(p4est_iter_face_info_t *info, int which_side, int which_hang){

    p4est_t                 *p4est;
    p4est_iter_face_side_t  *side;
    sc_array_t              *sides;
    p4est_topidx_t          which_tree;
    p4est_locidx_t          local_id;
    p4est_tree_t            *tree;

    p4est = info->p4est;
    sides = &(info->sides);
    side = p4est_iter_fside_array_index_int (sides, which_side);
    which_tree = side->treeid;
    tree = p4est_tree_array_index (p4est->trees, which_tree); 
    local_id = side->is.hanging.quadid[which_hang] + tree->quadrants_offset; 

    return local_id;
}


p4est_locidx_t p4est_utilities_get_local_id_face_full(p4est_iter_face_info_t *info, int which_side){
    
    p4est_t                 *p4est;
    p4est_iter_face_side_t  *side;
    sc_array_t              *sides;
    p4est_topidx_t          which_tree;
    p4est_locidx_t          local_id;
    p4est_tree_t            *tree;

    p4est = info->p4est;
    sides = &(info->sides);
    side = p4est_iter_fside_array_index_int (sides, which_side);
    which_tree = side->treeid;
    tree = p4est_tree_array_index (p4est->trees, which_tree); 
    local_id = side->is.full.quadid + tree->quadrants_offset;

    return local_id;
}


void p4est_utilities_mesh_stats(p4est_t *p4est){

    ctx_t              *d_ctx = (ctx_t *) p4est->user_pointer;
    p4est_topidx_t      t;
    p4est_topidx_t      first_local_tree = p4est->first_local_tree;
    p4est_topidx_t      last_local_tree = p4est->last_local_tree;
    sc_array_t         *trees = p4est->trees;
    p4est_tree_t       *tree;
    sc_array_t         *quadrants;
    p4est_quadrant_t   *quad;
    //quad_data_t        *data;
    p4est_locidx_t      si,n_quads;
    
    int l;

    int cells_per_level[d_ctx->d_grid.max_level+1];
    int cells_per_level_sum[d_ctx->d_grid.max_level+1];
    int cell_total = 0;
    //int cell_type;
    
    
    for(l=0;l<=d_ctx->d_grid.max_level;l++){
        cells_per_level[l] = 0;
    }
  
    
    for (t = first_local_tree; t <= last_local_tree; t++) {
        tree = p4est_tree_array_index(trees,t);
        quadrants = &(tree->quadrants);
        n_quads = (p4est_locidx_t) quadrants->elem_count;
        
        for (si = 0; si < n_quads; si++) {
            quad = p4est_quadrant_array_index(quadrants,si);
            
            //TODO: count cell types
            //TODO: data = (quad_data_t *) quad->p.user_data;
            //TODO: cell_type = data->type;
            ++cells_per_level[quad->level];
        }
    }


    MPI_Reduce(&cells_per_level,&cells_per_level_sum,(d_ctx->d_grid.max_level+1),MPI_INT,MPI_SUM,0,p4est->mpicomm);

    for(l=0;l<=d_ctx->d_grid.max_level;l++){
        cell_total +=cells_per_level_sum[l];
    }
        
    
    if(d_ctx->d_mpi.rank==0){
        printf("[regrid] Cell Counts: %d\n",cell_total);
        printf("         ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n");
        for(l=0;l<=d_ctx->d_grid.max_level;l++){
            if(cells_per_level_sum[l]){
                printf("         Level[%02d]: %d\n",l,cells_per_level_sum[l]);
            }
        }
        printf("         ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n");
    }
    
}
