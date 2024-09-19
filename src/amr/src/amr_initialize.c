/**
 * \file   amr_initialize.c
 * \ingroup amr_group
 * \author akirby
 *
 * \brief Initialization functions for the AMR code module
 *
 * Created on August 15, 2018, 1:59 PM
 */

/* header files */
#include "amr_initialize.h"

void initialize(int argc,char** argv,ctx_t *d_ctx,p4est_t **p4est,p4est_connectivity_t **conn){
    int noinput = 0;

    /* read the command line */
    if(argc < 2){
        noinput = 1;
        if(d_ctx->d_mpi.rank==0){
            printf("\n USAGE:\n");
            printf("   hpc_driver <input.file>\n\n");
            printf("   Please provide an input file\n "
                    "  Sample provided: sample.input.amr\n");
        }
    }else{
        if(d_ctx->d_mpi.rank==0){
            printf("[ AMR ] Input file: %s\n",argv[1]);
        }
        strcpy(d_ctx->d_initialize.input_file,argv[1]);
    }

    /* build wrk directories */
    utilities_create_directories(d_ctx->d_mpi.rank,d_ctx->d_mpi.comm);

    /* read the input file */
    initialize_inputs(d_ctx->d_initialize.input_file,d_ctx,noinput);

    /* construct the grid  */
    initialize_solver(d_ctx,p4est,conn);

    /* save inputs to file */
    initialize_save_inputs(d_ctx->d_initialize.input_file,d_ctx,noinput);
}

void initialize_inputs(char *filename,ctx_t *d_ctx,int noinput){
    /* default all inputs */
    d_ctx->d_initialize.restart_flag = 0;
    strcpy(d_ctx->d_initialize.restart_file,"amr_soln_0000000.bin");

    d_ctx->d_grid.scale = 1.0;
    d_ctx->d_grid.min_dx = 1.0;
    d_ctx->d_grid.max_level = 0;
    d_ctx->d_grid.min_level = 0;
    d_ctx->d_grid.nelem_subgrid = NXPATCH;
    d_ctx->d_grid.xlo[0] = 0.0;
    d_ctx->d_grid.xlo[1] = 0.0;
    d_ctx->d_grid.xlo[2] = 0.0;
    d_ctx->d_grid.xhi[0] = 1.0;
    d_ctx->d_grid.xhi[1] = 1.0;
    d_ctx->d_grid.xhi[2] = 1.0;
    d_ctx->d_grid.periodic[0] = 0;
    d_ctx->d_grid.periodic[1] = 0;
    d_ctx->d_grid.periodic[2] = 0;
    d_ctx->d_grid.regrid_interval = 0;
    d_ctx->d_grid.unstructured_flag = 0;

    d_ctx->d_grid.refined_box = 0;
    d_ctx->d_grid.box_xlo[0] = 0.0;
    d_ctx->d_grid.box_xlo[1] = 0.0;
    d_ctx->d_grid.box_xlo[2] = 0.0;
    d_ctx->d_grid.box_xhi[0] = 0.0;
    d_ctx->d_grid.box_xhi[1] = 0.0;
    d_ctx->d_grid.box_xhi[2] = 0.0;

    d_ctx->d_grid.construct_grid = 1;
    d_ctx->d_grid.nregrid_pts = 0;
    d_ctx->d_grid.regrid_width = 1.0;
    d_ctx->d_grid.regrid_xyzh = NULL;
    strcpy(d_ctx->d_grid.unst_grid_file,"hole_3d_cubit.inp");

    d_ctx->d_solution.checkpoint_interval = 0;
    d_ctx->d_vtk.visualization_interval = 0;

    d_ctx->d_external.ncell = 0;
    d_ctx->d_external.nrecv = 0;
    d_ctx->d_external.nsend = 0;
    d_ctx->d_external.nghost = 0;
    d_ctx->d_external.nfringe = 0;
    d_ctx->d_external.soln = NULL;
    d_ctx->d_external.recv = NULL;
    d_ctx->d_external.send = NULL;
    d_ctx->d_external.sizes = NULL;
    d_ctx->d_external.cell_geom = NULL;
    d_ctx->d_external.cell_info = NULL;
    d_ctx->d_external.face_info = NULL;

    /* read inputs from file */
    if(noinput==0){
        utilities_find_keyword_integer(filename,"restart_flag:",&d_ctx->d_initialize.restart_flag);
        utilities_find_keyword_string(filename,"restart_file:",d_ctx->d_initialize.restart_file);

        utilities_find_keyword_integer(filename,"unstructured_flag:",&d_ctx->d_grid.unstructured_flag);
        utilities_find_keyword_string(filename,"unstructured_file:",d_ctx->d_grid.unst_grid_file);
        utilities_find_keyword_three_integers(filename,"periodic_flag:",&d_ctx->d_grid.periodic[0],&d_ctx->d_grid.periodic[1],&d_ctx->d_grid.periodic[2]);
        utilities_find_keyword_three_doubles(filename,"domain_lo:",&d_ctx->d_grid.xlo[0],&d_ctx->d_grid.xlo[1],&d_ctx->d_grid.xlo[2]);
        utilities_find_keyword_three_doubles(filename,"domain_hi:",&d_ctx->d_grid.xhi[0],&d_ctx->d_grid.xhi[1],&d_ctx->d_grid.xhi[2]);
        utilities_find_keyword_double(filename,"min_dx:",&d_ctx->d_grid.min_dx);
        utilities_find_keyword_integer(filename,"max_amr_level:",&d_ctx->d_grid.max_level);
        utilities_find_keyword_integer(filename,"min_amr_level:",&d_ctx->d_grid.min_level);
        utilities_find_keyword_integer(filename,"refined_box:",&d_ctx->d_grid.refined_box);
        utilities_find_keyword_three_doubles(filename,"box_lo:",&d_ctx->d_grid.box_xlo[0],&d_ctx->d_grid.box_xlo[1],&d_ctx->d_grid.box_xlo[2]);
        utilities_find_keyword_three_doubles(filename,"box_hi:",&d_ctx->d_grid.box_xhi[0],&d_ctx->d_grid.box_xhi[1],&d_ctx->d_grid.box_xhi[2]);
        utilities_inputs_amr_message(d_ctx,d_ctx->d_mpi.rank);
    }
}

void initialize_save_inputs(char *filename,ctx_t *d_ctx,int noinput){
    char *backup_name;
    FILE *fp;

    /* save final inputs */
    if(noinput){
        backup_name = malloc(strlen("sample.input.amr")+1);
        strcpy(backup_name,"sample.input.amr");
    }else{
        backup_name = malloc(strlen("WRK/saved.")+strlen(filename)+1);
        strcpy(backup_name,"WRK/saved.");
        strcat(backup_name,filename);
    }

    if(d_ctx->d_mpi.rank==0){
        fp = fopen(backup_name,"w");

        fprintf(fp,"<INITIALIZATION>\n");
        fprintf(fp,"restart_flag: %d\n",d_ctx->d_initialize.restart_flag);
        fprintf(fp,"restart_file: %s\n",d_ctx->d_initialize.restart_file);
        fprintf(fp,"</INITIALIZATION>\n\n");

        fprintf(fp,"<GRID>\n");
        fprintf(fp,"unstructured_flag: %d\n",d_ctx->d_grid.unstructured_flag);
        fprintf(fp,"unstructured_file: %s\n",d_ctx->d_grid.unst_grid_file);
        fprintf(fp,"periodic_flag: %d %d %d\n",d_ctx->d_grid.periodic[0],d_ctx->d_grid.periodic[1],d_ctx->d_grid.periodic[2]);
        fprintf(fp,"domain_lo: %f %f %f\n",d_ctx->d_grid.xlo[0],d_ctx->d_grid.xlo[1],d_ctx->d_grid.xlo[2]);
        fprintf(fp,"domain_hi: %f %f %f\n",d_ctx->d_grid.xhi[0],d_ctx->d_grid.xhi[1],d_ctx->d_grid.xhi[2]);
        fprintf(fp,"min_dx: %f\n",d_ctx->d_grid.min_dx);
        fprintf(fp,"max_amr_level: %d\n",d_ctx->d_grid.max_level);
        fprintf(fp,"min_amr_level: %d\n",d_ctx->d_grid.min_level);
        fprintf(fp,"refined_box: %d\n",d_ctx->d_grid.refined_box);
        fprintf(fp,"box_lo: %f %f %f\n",d_ctx->d_grid.box_xlo[0],d_ctx->d_grid.box_xlo[1],d_ctx->d_grid.box_xlo[2]);
        fprintf(fp,"box_hi: %f %f %f\n",d_ctx->d_grid.box_xhi[0],d_ctx->d_grid.box_xhi[1],d_ctx->d_grid.box_xhi[2]);
        fprintf(fp,"</GRID>\n\n");

        fclose(fp);
        free(backup_name);
    }

    if(noinput){
        /* exit program when no input file is provided */
        exit(0);
    }
}

void initialize_solver(ctx_t *d_ctx,p4est_t **p4est,p4est_connectivity_t **conn){
    int load_data = 1;
    int autopartition = 1;
    int broadcasthead = 0;

    if(d_ctx->d_initialize.restart_flag){
        /* load the restart file */
        *p4est =
            p4est_load_ext(d_ctx->d_initialize.restart_file,
                           d_ctx->d_mpi.comm,sizeof(quad_data_t),
                           load_data,autopartition,broadcasthead,d_ctx,conn);

        /* get istep from filename */
        char *p = d_ctx->d_initialize.restart_file;
        while(*p){
            if (isdigit(*p)){
                long val = strtol(p,&p,10);
                d_ctx->istep = (int) val;
            } else {
                p++;
            }
        }
        utilities_restart_amr_message(d_ctx->d_initialize.restart_file,
                                      d_ctx->istep,
                                      d_ctx->d_mpi.rank);

        /* call external solver initialization callback function */
        external_initialize_solver(d_ctx);
    }else{
        /* build new grid */
        initialize_grid(d_ctx,conn);

        /* call external solver initialization callback function */
        external_initialize_solver(d_ctx);

        *p4est =
            p4est_new_ext(d_ctx->d_mpi.comm,                 /* mpi communicator */
                         *conn,                              /* p4est connectivity */
                          0,                                 /* minimum quadrants per MPI process */
                          d_ctx->d_grid.min_level,           /* minimum level of mesh refinement */
                          1,                                 /* fill the forest with a uniform mesh instead of the coarsest possible one */
                          sizeof(quad_data_t),               /* user quadrant data type */
                          initialize_quadrant_data,          /* p4est call back function from user for data initialization */
                          d_ctx);                            /* user data pointer */

        d_ctx->istep = 0;
    }
}

void initialize_grid(ctx_t *d_ctx,p4est_connectivity_t **conn){
    double x_length;
    double y_length;
    double z_length;
    double mesh_scale;

    int nelem_cLevel_x;
    int nelem_cLevel_y;
    int nelem_cLevel_z;

    int i;
    int nelem;

    /* determine the coarse level size */
    x_length = d_ctx->d_grid.xhi[0] - d_ctx->d_grid.xlo[0];
    y_length = d_ctx->d_grid.xhi[1] - d_ctx->d_grid.xlo[1];
    z_length = d_ctx->d_grid.xhi[2] - d_ctx->d_grid.xlo[2];

    mesh_scale = d_ctx->d_grid.min_dx * pow(2.0,d_ctx->d_grid.max_level);
    d_ctx->d_grid.scale = mesh_scale;
    nelem_cLevel_x = ceil(x_length/mesh_scale);
    nelem_cLevel_y = ceil(y_length/mesh_scale);
    nelem_cLevel_z = ceil(z_length/mesh_scale);

    d_ctx->d_grid.xhi[0] = d_ctx->d_grid.xlo[0] + mesh_scale*(double)nelem_cLevel_x;
    d_ctx->d_grid.xhi[1] = d_ctx->d_grid.xlo[1] + mesh_scale*(double)nelem_cLevel_y;
    d_ctx->d_grid.xhi[2] = d_ctx->d_grid.xlo[2] + mesh_scale*(double)nelem_cLevel_z;

    /* display grid information */
    if(d_ctx->d_mpi.rank==0){
        printf("[ amr ] AMR Levels (min,max): %d,%d\n",d_ctx->d_grid.min_level,d_ctx->d_grid.max_level);

        if(!d_ctx->d_grid.unstructured_flag){
            printf("[ amr ] Grid Sizes (min,max): %f,%f\n",d_ctx->d_grid.min_dx,mesh_scale);
#ifdef P4_TO_P8
            nelem = nelem_cLevel_x*nelem_cLevel_y*nelem_cLevel_z;
            printf("[ amr ] Coarsest Level: %d elements\n"
                   "        NX: %d\n"
                   "        NY: %d\n"
                   "        NZ: %d\n"
                   ,nelem,nelem_cLevel_x,nelem_cLevel_y,nelem_cLevel_z);
            printf("[ amr ] Resizing Domain: \n"
                   "        xlo: %f   xhi: %f \n"
                   "        ylo: %f   yhi: %f \n"
                   "        zlo: %f   zhi: %f \n"
                   ,d_ctx->d_grid.xlo[0],d_ctx->d_grid.xhi[0]
                   ,d_ctx->d_grid.xlo[1],d_ctx->d_grid.xhi[1]
                   ,d_ctx->d_grid.xlo[2],d_ctx->d_grid.xhi[2]);
#else
            nelem = nelem_cLevel_x*nelem_cLevel_y;
            printf("[ amr ] Coarsest Level: %d elements\n"
                   "        NX: %d\n"
                   "        NY: %d\n"
                   ,nelem,nelem_cLevel_x,nelem_cLevel_y);
            printf("[ amr ] Resizing Domain:\n"
                   "        xlo: %f   xhi: %f \n"
                   "        ylo: %f   yhi: %f \n"
                   ,d_ctx->d_grid.xlo[0],d_ctx->d_grid.xhi[0]
                   ,d_ctx->d_grid.xlo[1],d_ctx->d_grid.xhi[1]);
#endif
        }else{
            printf("[ amr ] Unstructured Grid file: %s\n",d_ctx->d_grid.unst_grid_file);
        }
    }

    /* create a new forest */
    if(!d_ctx->d_grid.unstructured_flag){

        /* Cartesian domain */
#ifndef P4_TO_P8
        *conn = p4est_connectivity_new_brick(
            nelem_cLevel_x,
            nelem_cLevel_y,
            d_ctx->d_grid.periodic[0],
            d_ctx->d_grid.periodic[1]);
#else
        *conn = p8est_connectivity_new_brick(
            nelem_cLevel_x,
            nelem_cLevel_y,
            nelem_cLevel_z,
            d_ctx->d_grid.periodic[0],
            d_ctx->d_grid.periodic[1],
            d_ctx->d_grid.periodic[2]);
#endif

        /* scale mesh and translate */
        for(i=0;i<(*conn)->num_vertices;++i){
            (*conn)->vertices[3*i+0] = (*conn)->vertices[3*i+0]*mesh_scale + d_ctx->d_grid.xlo[0];
            (*conn)->vertices[3*i+1] = (*conn)->vertices[3*i+1]*mesh_scale + d_ctx->d_grid.xlo[1];
            (*conn)->vertices[3*i+2] = (*conn)->vertices[3*i+2]*mesh_scale + d_ctx->d_grid.xlo[2];
        }

    }else{

        /* read unstructured grid */
        *conn = p4est_connectivity_read_inp(d_ctx->d_grid.unst_grid_file);

#ifdef P4EST_WITH_METIS
        p4est_connectivity_reorder(d_ctx->d_mpi.comm,0,*conn,P4EST_CONNECT_FULL);
#endif /* P4EST_WITH_METIS */
    }
}


void initialize_quadrant_data(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *q){
    double h;
    double xyz[3*P4EST_CHILDREN];
    int display_quad = 0;


    /* access our simulation context data */
    ctx_t *d_ctx = (ctx_t *) p4est->user_pointer;

    /* access our per quadrant data pointer */
    quad_data_t *data = (quad_data_t *) q->p.user_data;

    /* element size */
    h = d_ctx->d_grid.scale * P4EST_QUADRANT_LEN(q->level) / P4EST_ROOT_LEN;
    p4est_utilities_quad_coordinates(p4est,which_tree,q,xyz,display_quad);

    /* initialize the solution data in the quadrant via call-back function */
    external_initialize_quadrant(d_ctx->d_grid.dim,d_ctx->d_grid.unstructured_flag,NFIELDS,NXPATCH,&data->type,h,xyz,data->soln);
}