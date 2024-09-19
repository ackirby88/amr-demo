/**
 * \file   amr_external.c
 * \ingroup amr_group
 * \author akirby
 *
 * \brief External wrapper functions to call PHYSICS code module.
 *
 * Created on August 20, 2018, 11:26 AM
 */

/* header files */
#include "amr_external.h"

void external_initialize_solver(ctx_t *d_ctx){

    int nfields = NFIELDS;
    int nfringe = NFRINGE;
    int nelem_subgrid = NXPATCH;
    int extern_comm_cell_size = EXTERN_COMM_CELL_SIZE;
    int extern_comm_face_size = EXTERN_COMM_FACE_SIZE;
    int extern_comm_edge_size = EXTERN_COMM_EDGE_SIZE;
    int extern_comm_node_size = EXTERN_COMM_NODE_SIZE;
    int extern_comm_mpi_size  = EXTERN_COMM_MPI_SIZE;

    /* external physics kernel initialize callback function */
    hpc_amr_initialize_solver(  &d_ctx->d_mpi.rank,
                                &d_ctx->d_grid.dim,
                                &d_ctx->d_grid.scale,
                                &d_ctx->d_grid.max_level,
                                &nfields,&nfringe,&nelem_subgrid,
                                &extern_comm_cell_size,&extern_comm_face_size,
                                &extern_comm_edge_size,&extern_comm_node_size,
                                &extern_comm_mpi_size);

}

void external_initialize_quadrant(int dim,int unstructured,int numfields,int nsub_elem,int *cell_type,double elem_h,double *elem_xyz,double *elem_soln){

    /* external physics kernel quadrant initialize callback function */
    hpc_amr_initialize_quadrant(&dim,&unstructured,&numfields,&nsub_elem,cell_type,&elem_h,elem_xyz,elem_soln);

}

void external_regrid_points(int *construct_grid,int *nregrid_pts,double *regrid_width,double *regrid_xyzh){

    /* external physics kernel regrid points callback function */
    hpc_amr_regrid_points(construct_grid,nregrid_pts,regrid_width,regrid_xyzh);

}

void external_tag_feature(int dim,int unstructured,int numfields,int nsub_elem,int *cell_type,double elem_h,double *elem_xyz,double *elem_soln,int *tag){

    /* external physics kernel quadrant initialize callback function */
    hpc_amr_tag_feature(&dim,&unstructured,&numfields,&nsub_elem,cell_type,&elem_h,elem_xyz,elem_soln,tag);

}

void external_evolve(ctx_t *d_ctx,p4est_t **p4est,int *istep,int *regrid,int *visualize,int *checkpoint,int *evolve_solution){

    p4est_ghost_t *ghost;
    quad_data_t   *ghost_data;
    int nfields = NFIELDS;

    /* initialize simulation parameters */
    *regrid = 0;
    *visualize = 0;
    *checkpoint = 0;

    /* create the ghost quadrants */
    ghost = p4est_ghost_new(*p4est,P4EST_CONNECT_FULL);

    /* allocate ghost data to exchange cell types */
    ghost_data = P4EST_ALLOC(quad_data_t,ghost->ghosts.elem_count);

    /* exchange ghost data to synchronize */
    p4est_ghost_exchange_data(*p4est,ghost,ghost_data);

    /* setup external solver buffers and mpi information */
    external_setup(*p4est,ghost,ghost_data);

    /* copy p4est solution data to the external solution vector */
    external_p4est_to_soln(*p4est);

    /* call external physics kernel evolve callback function */
    hpc_amr_evolve_solver( &d_ctx->d_grid.dim,
                           &d_ctx->d_grid.unstructured_flag,
                           &nfields,
                           &d_ctx->d_grid.nelem_subgrid,
                           &d_ctx->d_external.ncell,
                           &d_ctx->d_external.nghost,
                           &d_ctx->d_external.nface,
                           &d_ctx->d_external.dof,
                            d_ctx->d_external.cell_info,
                            d_ctx->d_external.face_info,
                           &d_ctx->d_external.nsend,
                           &d_ctx->d_external.nrecv,
                            d_ctx->d_external.send,
                            d_ctx->d_external.recv,
                            d_ctx->d_external.sizes,
                            d_ctx->d_external.cell_geom,
                            d_ctx->d_external.soln,
                            istep,regrid,visualize,checkpoint,evolve_solution,
                            d_ctx->d_mpi.comm,
                           &d_ctx->total_time,
                           &d_ctx->compute_time,
                           &d_ctx->residual_time);

    /* copy external solution vector to the p4est solution data */
    external_soln_to_p4est(*p4est);

    /* reset the p4est ghost data structure */
    P4EST_FREE(ghost_data);
    p4est_ghost_destroy(ghost);
}

void external_setup(p4est_t *p4est,p4est_ghost_t *ghost,quad_data_t *ghost_data){

    /* allocate external soln vector */
    external_allocate_soln(p4est,ghost,ghost_data);

    /* fill mpi exchange information */
    external_setup_ghost_exchange_data(p4est,ghost);

    /* fill external cell_info and face_info through p4est iterate */
    p4est_iterate(p4est,ghost,NULL,
            external_cellinfo_callback,
            external_faceinfo_callback,
#ifdef P4_TO_P8
            NULL,
#endif
            NULL);
}

void external_p4est_to_soln(p4est_t *p4est){
    p4est_topidx_t     t;
    ctx_t             *d_ctx = (ctx_t *) p4est->user_pointer;
    p4est_topidx_t     first_local_tree = p4est->first_local_tree;
    p4est_topidx_t     last_local_tree = p4est->last_local_tree;
    sc_array_t        *trees = p4est->trees;
    p4est_tree_t      *tree;
    sc_array_t        *quadrants;
    p4est_quadrant_t  *quad;
    quad_data_t       *data;
    p4est_locidx_t     si,n_quads;
    int                i,ind,nftm;

#ifndef P4_TO_P8
    nftm = NFIELDS*NXPATCH*NXPATCH;
#else
    nftm = NFIELDS*NXPATCH*NXPATCH*NXPATCH;
#endif

    /* copy quad_data->soln solution to d_ctx->d_external.soln  */
    ind = 0;
    for(t = first_local_tree; t <= last_local_tree; t++) {
        tree = p4est_tree_array_index (trees, t);
        quadrants = &(tree->quadrants);
        n_quads = (p4est_locidx_t) quadrants->elem_count;
        for (si = 0; si < n_quads; si++) {

            quad = p4est_quadrant_array_index (quadrants, si);
            data = (quad_data_t *) quad->p.user_data;

            for(i = 0; i < nftm; i++){
                d_ctx->d_external.soln[ind+i] = data->soln[i];
            }
            ind += nftm;
        }
    }
}

void external_soln_to_p4est(p4est_t *p4est){

    p4est_topidx_t     t;
    ctx_t             *d_ctx = (ctx_t *) p4est->user_pointer;
    p4est_topidx_t     first_local_tree = p4est->first_local_tree;
    p4est_topidx_t     last_local_tree = p4est->last_local_tree;
    sc_array_t        *trees = p4est->trees;
    p4est_tree_t      *tree;
    sc_array_t        *quadrants;
    p4est_quadrant_t  *quad;
    quad_data_t       *data;
    p4est_locidx_t     si,n_quads;
    int                i,ind,nftm;

#ifndef P4_TO_P8
    nftm = NFIELDS*NXPATCH*NXPATCH;
#else
    nftm = NFIELDS*NXPATCH*NXPATCH*NXPATCH;
#endif

    /* copy d_ctx->d_external.soln to quad_data->soln */
    ind = 0;
    for(t = first_local_tree; t <= last_local_tree; t++) {
        tree = p4est_tree_array_index (trees, t);
        quadrants = &(tree->quadrants);
        n_quads = (p4est_locidx_t) quadrants->elem_count;
        for (si = 0; si < n_quads; si++) {

            quad = p4est_quadrant_array_index (quadrants, si);
            data = (quad_data_t *) quad->p.user_data;

            for(i = 0; i < nftm; i++){
                data->soln[i] = d_ctx->d_external.soln[ind+i];
            }
            ind += nftm;
        }
    }
}

void external_allocate_soln(p4est_t *p4est,p4est_ghost_t *ghost,quad_data_t *ghost_data){
    int i,j,dof,loc_nfields;
    quad_data_t *data;
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;

    /* free cell info, face info, and solution */
    external_deallocate_soln(p4est);

    /* count total cells real and ghost */
    const int ncell = p4est->local_num_quadrants + (int) ghost->ghosts.elem_count;

    d_ctx->d_external.cell_info = P4EST_ALLOC(int,EXTERN_COMM_CELL_SIZE*ncell);
    d_ctx->d_external.cell_geom = P4EST_ALLOC(double,3*P4EST_CHILDREN*ncell);
    d_ctx->d_external.face_info = P4EST_ALLOC(int,EXTERN_COMM_FACE_SIZE*6*ncell);

    d_ctx->d_external.nsend = (int) ghost->mirror_proc_offsets[p4est->mpisize];
    d_ctx->d_external.nrecv = (int) ghost->ghosts.elem_count;

    d_ctx->d_external.send = P4EST_ALLOC(int,EXTERN_COMM_MPI_SIZE*d_ctx->d_external.nsend);
    d_ctx->d_external.recv = P4EST_ALLOC(int,EXTERN_COMM_MPI_SIZE*d_ctx->d_external.nrecv);
    d_ctx->d_external.sizes = P4EST_ALLOC(int,ncell);

    /* initialize all info */
    for(i=0;i<ncell;i++){
        for(j=0;j<EXTERN_COMM_CELL_SIZE;j++){
            d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*i+j] = -1;
        }
        for(j=0;j<3*P4EST_CHILDREN;j++){
            d_ctx->d_external.cell_geom[3*P4EST_CHILDREN*i+j] = 0.0;
        }
    }

    for(i=0;i<6*ncell;i++){
        for(j=0;j<EXTERN_COMM_FACE_SIZE;j++){
            d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*i+j] = -1;
        }
    }

#ifndef P4_TO_P8
    loc_nfields = NFIELDS*NXPATCH*NXPATCH;
#else
    loc_nfields = NFIELDS*NXPATCH*NXPATCH*NXPATCH;
#endif
    dof = loc_nfields*p4est->local_num_quadrants;
    d_ctx->d_external.dof = dof;


    for(i=0;i<ghost->ghosts.elem_count;++i){
        data = (quad_data_t *) & ghost_data[i];

        d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*(p4est->local_num_quadrants+i)+0] = data->type;
        d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*(p4est->local_num_quadrants+i)+2] = dof;
        d_ctx->d_external.sizes[p4est->local_num_quadrants+i] = loc_nfields;
        dof += loc_nfields;

    }

    d_ctx->d_external.soln = P4EST_ALLOC(double,dof);
    d_ctx->d_external.ncell = p4est->local_num_quadrants;
    d_ctx->d_external.nghost = (int) ghost->ghosts.elem_count;
    d_ctx->d_external.nface = 0;

    //printf("Rank[%d] has # ghosts: %d\n",d_ctx->d_mpi.rank,ghost->ghosts.elem_count);
}

void external_setup_ghost_exchange_data(p4est_t *p4est, p4est_ghost_t *ghost){
    p4est_quadrant_t   *mirror;
    ctx_t              *d_ctx = (ctx_t *) p4est->user_pointer;
    const int           num_procs = p4est->mpisize;
    int                 q,g;
    p4est_locidx_t      ng_excl,ng_incl,ng,theg;
    p4est_locidx_t      mirr;

    int nrecv = 0;
    int nsend = 0;

    /* receive data from other processors */
    ng_excl = 0;
    for (q = 0; q < num_procs; ++q) {
        ng_incl = ghost->proc_offsets[q+1];
        ng = ng_incl - ng_excl;
        P4EST_ASSERT(ng >= 0);
        if(ng > 0){
            for(g=0;g<ng;++g){
                d_ctx->d_external.recv[EXTERN_COMM_MPI_SIZE*nrecv+0] = p4est->local_num_quadrants + ng_excl + g; //cell_index
                d_ctx->d_external.recv[EXTERN_COMM_MPI_SIZE*nrecv+1] = q; //proc_id
                ++nrecv;
            }
            ng_excl = ng_incl;
        }
    }
    P4EST_ASSERT (ng_excl == (p4est_locidx_t) ghost->ghosts.elem_count);

    /* send data to other processors */
    ng_excl = 0;
    for (q = 0; q < num_procs; ++q) {
        ng_incl = ghost->mirror_proc_offsets[q + 1];
        ng = ng_incl - ng_excl;
        P4EST_ASSERT (ng >= 0);
        if (ng > 0) {
            for (theg = 0; theg < ng; ++theg) {
                mirr = ghost->mirror_proc_mirrors[ng_excl + theg];
                mirror = p4est_quadrant_array_index (&ghost->mirrors, mirr);
                d_ctx->d_external.send[EXTERN_COMM_MPI_SIZE*nsend+0] = mirror->p.piggy3.local_num; //cell_index
                d_ctx->d_external.send[EXTERN_COMM_MPI_SIZE*nsend+1] = q; //proc_id
                ++nsend;
            }
            ng_excl = ng_incl;
        }
    }
    //printf("Rank[%d]: nsend[%d] recv[%d]\n",d_ctx->d_mpi.rank,nsend,nrecv);
}

void external_deallocate_soln(p4est_t *p4est){
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;

    if(d_ctx->d_external.cell_info) P4EST_FREE(d_ctx->d_external.cell_info);
    if(d_ctx->d_external.cell_geom) P4EST_FREE(d_ctx->d_external.cell_geom);
    if(d_ctx->d_external.face_info) P4EST_FREE(d_ctx->d_external.face_info);
    if(d_ctx->d_external.send)      P4EST_FREE(d_ctx->d_external.send);
    if(d_ctx->d_external.recv)      P4EST_FREE(d_ctx->d_external.recv);
    if(d_ctx->d_external.sizes)     P4EST_FREE(d_ctx->d_external.sizes);
    if(d_ctx->d_external.soln)      P4EST_FREE(d_ctx->d_external.soln);
}

/* volume callback function to fill cell info */
void external_cellinfo_callback(p4est_iter_volume_info_t *info, void *user_data){
    p4est_t          *p4est = info->p4est;
    p4est_quadrant_t *q = info->quad;
    quad_data_t      *data = (quad_data_t *) q->p.user_data;
    ctx_t            *d_ctx = (ctx_t *) p4est->user_pointer;
    p4est_locidx_t    local_id;
    p4est_topidx_t    which_tree = info->treeid;

    int loc_nfields;
    int display_quad = 0;

    local_id = p4est_utilities_get_local_id_volume(info);
    which_tree = info->treeid;

#ifndef P4_TO_P8
    loc_nfields = NFIELDS*NXPATCH*NXPATCH;
#else
    loc_nfields = NFIELDS*NXPATCH*NXPATCH*NXPATCH;
#endif
    d_ctx->d_external.sizes[local_id] = loc_nfields;

    d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*local_id+0] = data->type;
    d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*local_id+1] = q->level;
    d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*local_id+2] = loc_nfields*local_id; // soln index

    /* find the quadrant by searching over trees */
    p4est_utilities_quad_coordinates(p4est,which_tree,q,&d_ctx->d_external.cell_geom[3*P4EST_CHILDREN*local_id],display_quad);
}

/* face callback function to fill face info */
void external_faceinfo_callback(p4est_iter_face_info_t * info, void *user_data){

    p4est_t                 *p4est = info->p4est;
    ctx_t                   *d_ctx = (ctx_t *) p4est->user_pointer;
    p4est_iter_face_side_t  *side[2];
    sc_array_t              *sides = &(info->sides);
    p4est_locidx_t          local_id;
    p4est_topidx_t          which_tree;
    int                     j,hangside,fullside,hanging,boundary;
    int display_quad = 0;

    if(sides->elem_count > 1){
        boundary = 0;
    } else {
        boundary = 1;
    }

    side[0] = p4est_iter_fside_array_index_int(sides,0);
    if(!boundary){
        side[1] = p4est_iter_fside_array_index_int(sides,1);
    }

    /* face_info description */
    // face_info[0]-> face type: 1-boundary, 2-full-full, 3/5-full-hang
    // face_info[1]-> this cell side: 0-xlo, 1-xhi, ...
    // face_info[2]-> this cell id
    // face_info[3]-> nhbr cell id 1 (if present)
    // face_info[4]-> nhbr cell id 2 (if present)
    // face_info[5]-> nhbr cell id 3 (if present)
    // face_info[6]-> nhbr cell id 4 (if present)
    // face_info[7]-> nhbr cell 1 side (if present)
    // face_info[8]-> nhbr cell 2 side (if present)
    // face_info[9]-> nhbr cell 3 side (if present)
    // face_info[10]->nhbr cell 4 side (if present)

    hangside = 0;
    fullside = 0;
    hanging = 0;

    /* detect if a face is hanging */
    if(side[0]->is_hanging) {
        hangside = 0;
        fullside = 1;
        hanging = 1;
    }
    if(!boundary){
        if(side[1]->is_hanging){
            hangside = 1;
            fullside = 0;
            hanging = 1;
        }
    }

    /* if a boundary tag with 1, if full-full tag 2, and if hanging tag 1+2 or 1+4 */
    if(boundary){
        d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+0] = 1;
    }else if(hanging){
        d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+0] = 1+P4EST_HALF;
    }else{
        d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+0] = 2;
    }

    if(hanging){
        if(side[fullside]->is.full.is_ghost){
            local_id = side[fullside]->is.full.quadid+p4est->local_num_quadrants;
            which_tree = side[fullside]->treeid;
        } else {
            local_id = p4est_utilities_get_local_id_face_full(info,fullside);
            which_tree = p4est_iter_fside_array_index_int(sides,fullside)->treeid;
        }

        if(local_id > d_ctx->d_external.ncell+d_ctx->d_external.nghost){
            printf("Local id is out of bounds %d\n",local_id);
        }

        d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+1] = side[fullside]->face;
        d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+2] = local_id;

        d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*local_id+1] = side[fullside]->is.full.quad->level;

        p4est_utilities_quad_coordinates(p4est,which_tree,side[fullside]->is.full.quad,&d_ctx->d_external.cell_geom[3*P4EST_CHILDREN*local_id],display_quad);

        for (j = 0; j < P4EST_HALF; j++) {
            if(side[hangside]->is.hanging.is_ghost[j]){
                local_id = side[hangside]->is.hanging.quadid[j]+p4est->local_num_quadrants;
                which_tree = side[hangside]->treeid;
            }else{
                local_id = p4est_utilities_get_local_id_face_hang(info,hangside,j);
                which_tree = p4est_iter_fside_array_index_int(sides,hangside)->treeid;
            }

            if(local_id > d_ctx->d_external.ncell+d_ctx->d_external.nghost){
                printf("Local id is out of bounds %d\n",local_id);
            }

            d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+3+j] = local_id;
            d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+7+j] = side[hangside]->face;
            d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*local_id+1] = side[hangside]->is.hanging.quad[j]->level;

            p4est_utilities_quad_coordinates(p4est,which_tree,side[hangside]->is.hanging.quad[j],&d_ctx->d_external.cell_geom[3*P4EST_CHILDREN*local_id],display_quad);
        }
    }else{
        if (side[0]->is.full.is_ghost) {
            local_id = side[0]->is.full.quadid+p4est->local_num_quadrants;
            which_tree = side[0]->treeid;
        } else {
            local_id = p4est_utilities_get_local_id_face_full(info, 0);
            which_tree = p4est_iter_fside_array_index_int (sides, 0)->treeid;
        }

        if(local_id > d_ctx->d_external.ncell+d_ctx->d_external.nghost){
            printf("Local id is out of bounds %d\n",local_id);
        }

        d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+1] = side[0]->face;
        d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+2] = local_id;
        d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*local_id+1] = side[0]->is.full.quad->level;

        p4est_utilities_quad_coordinates(p4est,which_tree,side[0]->is.full.quad,&d_ctx->d_external.cell_geom[3*P4EST_CHILDREN*local_id],display_quad);

        if(!boundary){
            if (side[1]->is.full.is_ghost) {
                local_id = side[1]->is.full.quadid+p4est->local_num_quadrants;
                which_tree = side[1]->treeid;
            }else{
                local_id = p4est_utilities_get_local_id_face_full(info, 1);
                which_tree = p4est_iter_fside_array_index_int(sides, 1)->treeid;
            }

            if(local_id > d_ctx->d_external.ncell+d_ctx->d_external.nghost){
                printf("Local id is out of bounds %d\n",local_id);
            }

            d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+3] = local_id;
            d_ctx->d_external.face_info[EXTERN_COMM_FACE_SIZE*d_ctx->d_external.nface+7] = side[1]->face;
            d_ctx->d_external.cell_info[EXTERN_COMM_CELL_SIZE*local_id+1] = side[1]->is.full.quad->level;

            p4est_utilities_quad_coordinates(p4est,which_tree,side[1]->is.full.quad,&d_ctx->d_external.cell_geom[3*P4EST_CHILDREN*local_id],display_quad);
        }
    }
    ++d_ctx->d_external.nface;
}

void external_coarsen_operator(int dim,int unstructured,int numfields,int nsub_elem,int *cell_type,double elem_h,double *xyz,double **uc,double *soln){
    if(dim==2){
        hpc_amr_coarsen_operator_2d(&unstructured,&numfields,&nsub_elem,cell_type,&elem_h,xyz,uc[0],uc[1],uc[2],uc[3],soln);
    }else{
        hpc_amr_coarsen_operator_3d(&unstructured,&numfields,&nsub_elem,cell_type,&elem_h,xyz,uc[0],uc[1],uc[2],uc[3],uc[4],uc[5],uc[6],uc[7],soln);
    }
}

void external_refine_operator(int dim,int unstructured,int numfields,int nsub_elem,int *cell_type,double elem_h,double *xyz,double *soln,double **uc){
    if(dim==2){
        hpc_amr_refine_operator_2d(&unstructured,&numfields,&nsub_elem,cell_type,&elem_h,xyz,soln,uc[0],uc[1],uc[2],uc[3]);
    }else{
        hpc_amr_refine_operator_3d(&unstructured,&numfields,&nsub_elem,cell_type,&elem_h,xyz,soln,uc[0],uc[1],uc[2],uc[3],uc[4],uc[5],uc[6],uc[7]);
    }
}