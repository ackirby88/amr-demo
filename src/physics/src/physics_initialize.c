/**
 * \file   physics_initialize.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Initialize Physics interface parameters.
 * 
 * Created on August 15, 2018, 1:59 PM
 */


/* header files */
#include "physics_initialize.h"


void physics_initialize(int rank,int dim,double mesh_scale,int nlevels,
                        int nfields,int nfringe,int nelem_subgrid,
                        int extern_cell_size,int extern_face_size,
                        int extern_edge_size,int extern_node_size,
                        int extern_mpi_comm_size){
    
    
    int noinput = 0;
    char input_file[] = "input.physics";
    
    utilities_write_physics_message(rank);
    
    // read the input file
    if(rank==0){
        printf("[physics] Input file: %s\n",input_file);
    }
    
    if(access(input_file,F_OK ) == -1){
        noinput = 1;
        if(rank==0){
            printf("[physics] The physics input file does not exist!\n");
        }
    }
        
    
    
    d_physics_interface.nlevels = nlevels;
    d_physics_interface.coarse_grid_h = mesh_scale;
    
    d_physics_interface.dim = dim;
    d_physics_interface.nfields = nfields;
    d_physics_interface.nfringe = nfringe;
    d_physics_interface.nelem_subgrid = nelem_subgrid;
    d_physics_interface.extern_cell_size = extern_cell_size;
    d_physics_interface.extern_face_size = extern_face_size;
    d_physics_interface.extern_edge_size = extern_edge_size;
    d_physics_interface.extern_node_size = extern_node_size;
    d_physics_interface.extern_mpi_comm_size = extern_mpi_comm_size;
    
    physics_initialize_inputs(input_file,noinput);
    physics_initialize_save_inputs(input_file,rank,noinput);
    
}


void physics_initialize_inputs(char *filename,int noinput){
    
    int notfound;
    
    /* default all inputs */
    d_physics_fluid.initial_condition = 1;
    d_physics_simulation.time_scheme = 1;
    d_physics_simulation.time_steps = 1;
    d_physics_simulation.dt = 0.001;
    d_physics_fluid.alpha = 0.0;
    d_physics_fluid.beta = 0.0;
    d_physics_fluid.bc[0] = 1;
    d_physics_fluid.bc[1] = 1;
    d_physics_fluid.bc[2] = 1;
    d_physics_fluid.bc[3] = 1;
    d_physics_fluid.bc[4] = 1;
    d_physics_fluid.bc[5] = 1;
    d_physics_fluid.cfl = 1.0;
    d_physics_fluid.gamma = 1.4;
    d_physics_fluid.density = 1.0;
    d_physics_fluid.pressure = 1.0/1.4;
    d_physics_fluid.mach = 0.1;
    d_physics_fluid.amr_ntag_methods = 0;
    d_physics_fluid.amr_tag_tolerances = NULL;
    d_physics_fluid.amr_tag_methods = NULL;
    
    
    /* read inputs from file */
    if(!noinput){
        
        utilities_find_keyword_integer(filename,"checkpoint_interval:",&d_physics_simulation.checkpoint_interval);
        utilities_find_keyword_integer(filename,"visualization_interval:",&d_physics_simulation.visualization_interval);
        utilities_find_keyword_integer(filename,"regrid_interval:",&d_physics_simulation.regrid_interval);
        
        
        notfound = utilities_find_keyword_integer(filename,"amr_ntag_methods:",&d_physics_fluid.amr_ntag_methods);
        if(notfound==0){
            /* allocate and read tagging methods */
            
            
            switch(d_physics_fluid.amr_ntag_methods){
                case(0):
                    d_physics_fluid.amr_tag_methods = (int *) malloc(sizeof(int));
                    d_physics_fluid.amr_tag_tolerances = (double *) malloc(sizeof(double));
                    
                    d_physics_fluid.amr_tag_methods[0] = 0;
                    d_physics_fluid.amr_tag_tolerances[0] = 1.0;
                    break;
                case(1):
                    d_physics_fluid.amr_tag_methods = (int *) malloc(sizeof(int));
                    d_physics_fluid.amr_tag_tolerances = (double *) malloc(sizeof(double));
                    
                    utilities_find_keyword_integer(filename,"amr_tag_methods:",&d_physics_fluid.amr_tag_methods[0]);
                    utilities_find_keyword_double(filename,"amr_tag_tolerances:",&d_physics_fluid.amr_tag_tolerances[0]);
                    break;
                case(2):
                    d_physics_fluid.amr_tag_methods = (int *) malloc(2*sizeof(int));
                    d_physics_fluid.amr_tag_tolerances = (double *) malloc(2*sizeof(double));
                    
                    utilities_find_keyword_two_integers(filename,"amr_tag_methods:",&d_physics_fluid.amr_tag_methods[0],
                                                                                    &d_physics_fluid.amr_tag_methods[1]);
                    utilities_find_keyword_two_doubles(filename,"amr_tag_tolerances:",&d_physics_fluid.amr_tag_tolerances[0],
                                                                                      &d_physics_fluid.amr_tag_tolerances[1]);
                    break;
                case(3):
                    d_physics_fluid.amr_tag_methods = (int *) malloc(3*sizeof(int));
                    d_physics_fluid.amr_tag_tolerances = (double *) malloc(3*sizeof(double));
                    
                    utilities_find_keyword_three_integers(filename,"amr_tag_methods:",&d_physics_fluid.amr_tag_methods[0],
                                                                                      &d_physics_fluid.amr_tag_methods[1],
                                                                                      &d_physics_fluid.amr_tag_methods[2]);
                    utilities_find_keyword_three_doubles(filename,"amr_tag_tolerances:",&d_physics_fluid.amr_tag_tolerances[0],
                                                                                        &d_physics_fluid.amr_tag_tolerances[1],
                                                                                        &d_physics_fluid.amr_tag_tolerances[2]);
                    break;
                case(4):
                    d_physics_fluid.amr_tag_methods = (int *) malloc(4*sizeof(int));
                    d_physics_fluid.amr_tag_tolerances = (double *) malloc(4*sizeof(double));
                    
                    utilities_find_keyword_four_integers(filename,"amr_tag_methods:",&d_physics_fluid.amr_tag_methods[0],
                                                                                     &d_physics_fluid.amr_tag_methods[1],
                                                                                     &d_physics_fluid.amr_tag_methods[2],
                                                                                     &d_physics_fluid.amr_tag_methods[3]);
                    utilities_find_keyword_four_doubles(filename,"amr_tag_tolerances:",&d_physics_fluid.amr_tag_tolerances[0],
                                                                                       &d_physics_fluid.amr_tag_tolerances[1],
                                                                                       &d_physics_fluid.amr_tag_tolerances[2],
                                                                                       &d_physics_fluid.amr_tag_tolerances[3]);
                    break;
                case(5):
                    d_physics_fluid.amr_tag_methods = (int *) malloc(5*sizeof(int));
                    d_physics_fluid.amr_tag_tolerances = (double *) malloc(5*sizeof(double));
                    
                    utilities_find_keyword_five_integers(filename,"amr_tag_methods:",&d_physics_fluid.amr_tag_methods[0],
                                                                                     &d_physics_fluid.amr_tag_methods[1],
                                                                                     &d_physics_fluid.amr_tag_methods[2],
                                                                                     &d_physics_fluid.amr_tag_methods[3],
                                                                                     &d_physics_fluid.amr_tag_methods[4]);
                    utilities_find_keyword_five_doubles(filename,"amr_tag_tolerances:",&d_physics_fluid.amr_tag_tolerances[0],
                                                                                       &d_physics_fluid.amr_tag_tolerances[1],
                                                                                       &d_physics_fluid.amr_tag_tolerances[2],
                                                                                       &d_physics_fluid.amr_tag_tolerances[3],
                                                                                       &d_physics_fluid.amr_tag_tolerances[4]);
                    break;
            }
                    
            
            
        }
        
        
        
        utilities_find_keyword_integer(filename,"initial_condition:",&d_physics_fluid.initial_condition);
        utilities_find_keyword_integer(filename,"time_scheme:",&d_physics_simulation.time_scheme);
        utilities_find_keyword_integer(filename,"time_steps:",&d_physics_simulation.time_steps);
        utilities_find_keyword_double(filename,"dt:",&d_physics_simulation.dt);

        utilities_find_keyword_integer(filename,"bc_xlo:",&d_physics_fluid.bc[0]);
        utilities_find_keyword_integer(filename,"bc_xhi:",&d_physics_fluid.bc[1]);
        utilities_find_keyword_integer(filename,"bc_ylo:",&d_physics_fluid.bc[2]);
        utilities_find_keyword_integer(filename,"bc_yhi:",&d_physics_fluid.bc[3]);
        utilities_find_keyword_integer(filename,"bc_zlo:",&d_physics_fluid.bc[4]);
        utilities_find_keyword_integer(filename,"bc_zhi:",&d_physics_fluid.bc[5]);

        utilities_find_keyword_double(filename,"alpha:",&d_physics_fluid.alpha);
        utilities_find_keyword_double(filename,"beta:",&d_physics_fluid.beta);
        utilities_find_keyword_double(filename,"cfl:",&d_physics_fluid.cfl);
        utilities_find_keyword_double(filename,"gamma:",&d_physics_fluid.gamma);
        utilities_find_keyword_double(filename,"density:",&d_physics_fluid.density);
        utilities_find_keyword_double(filename,"pressure:",&d_physics_fluid.pressure);
        utilities_find_keyword_double(filename,"mach:",&d_physics_fluid.mach);
        
        // scale cfl depending on time scheme
        d_physics_fluid.cfl = 0.25*d_physics_fluid.cfl*d_physics_simulation.time_scheme;

    }
    
}

void physics_initialize_save_inputs(char *filename, int rank, int noinput){
    
    FILE *fp;
    char *backup_name;
    
    /* save final inputs */
    
    if(noinput){
        backup_name = malloc(strlen("sample.input.physics"));
        strcpy(backup_name,"sample.input.physics");
    }else{
        backup_name = malloc(strlen("WRK/saved.") + strlen(filename)+1);
        strcpy(backup_name,"WRK/saved.");
        strcat(backup_name,filename);
    }
    
    
    if(rank==0){
        
        fp = fopen(backup_name,"w");
        
        fprintf(fp,"<SIMULATION PROPERTIES>\n");
        fprintf(fp,"checkpoint_interval: %d\n",d_physics_simulation.checkpoint_interval);
        fprintf(fp,"visualization_interval: %d\n",d_physics_simulation.visualization_interval);
        fprintf(fp,"regrid_interval: %d\n",d_physics_simulation.regrid_interval);
        
        fprintf(fp,"amr_ntag_methods: %d\n",d_physics_fluid.amr_ntag_methods);
        switch(d_physics_fluid.amr_ntag_methods){
            case(0):
                fprintf(fp,"amr_tag_methods: %d\n",1);
                fprintf(fp,"amr_tag_tolerances: %f\n",0.0);
                break;
            case(1):
                fprintf(fp,"amr_tag_methods: %d\n",
                            d_physics_fluid.amr_tag_methods[0]);
                fprintf(fp,"amr_tag_tolerances: %f\n",
                            d_physics_fluid.amr_tag_tolerances[0]);
                break;
            case(2):
                fprintf(fp,"amr_tag_methods: %d %d\n",
                            d_physics_fluid.amr_tag_methods[0],
                            d_physics_fluid.amr_tag_methods[1]);
                fprintf(fp,"amr_tag_tolerances: %f %f\n",
                            d_physics_fluid.amr_tag_tolerances[0],
                            d_physics_fluid.amr_tag_tolerances[1]);
                break;
            case(3):
                fprintf(fp,"amr_tag_methods: %d %d %d\n",
                            d_physics_fluid.amr_tag_methods[0],
                            d_physics_fluid.amr_tag_methods[1],
                            d_physics_fluid.amr_tag_methods[2]);
                fprintf(fp,"amr_tag_tolerances: %f %f %f\n",
                            d_physics_fluid.amr_tag_tolerances[0],
                            d_physics_fluid.amr_tag_tolerances[1],
                            d_physics_fluid.amr_tag_tolerances[2]);
                break;
            case(4):
                fprintf(fp,"amr_tag_methods: %d %d %d %d\n",
                            d_physics_fluid.amr_tag_methods[0],
                            d_physics_fluid.amr_tag_methods[1],
                            d_physics_fluid.amr_tag_methods[2],
                            d_physics_fluid.amr_tag_methods[3]);
                fprintf(fp,"amr_tag_tolerances: %f %f %f %f\n",
                            d_physics_fluid.amr_tag_tolerances[0],
                            d_physics_fluid.amr_tag_tolerances[1],
                            d_physics_fluid.amr_tag_tolerances[2],
                            d_physics_fluid.amr_tag_tolerances[3]);
                break;
            case(5):
                fprintf(fp,"amr_tag_methods: %d %d %d %d %d\n",
                            d_physics_fluid.amr_tag_methods[0],
                            d_physics_fluid.amr_tag_methods[1],
                            d_physics_fluid.amr_tag_methods[2],
                            d_physics_fluid.amr_tag_methods[3],
                            d_physics_fluid.amr_tag_methods[4]);
                fprintf(fp,"amr_tag_tolerances: %f %f %f %f %f\n",
                            d_physics_fluid.amr_tag_tolerances[0],
                            d_physics_fluid.amr_tag_tolerances[1],
                            d_physics_fluid.amr_tag_tolerances[2],
                            d_physics_fluid.amr_tag_tolerances[3],
                            d_physics_fluid.amr_tag_tolerances[4]);
                break;
        }
        
        
        fprintf(fp,"time_scheme: %d\n",d_physics_simulation.time_scheme);
        fprintf(fp,"time_steps: %d\n",d_physics_simulation.time_steps);
        fprintf(fp,"dt: %f\n",d_physics_simulation.dt);
        fprintf(fp,"</SIMULATION PROPERTIES>\n\n");
        
        fprintf(fp,"<FLUID PROPERTIES>\n");
        fprintf(fp,"initial_condition: %d\n",d_physics_fluid.initial_condition);
        fprintf(fp,"bc_xlo: %d\n",d_physics_fluid.bc[0]);
        fprintf(fp,"bc_xhi: %d\n",d_physics_fluid.bc[1]);
        fprintf(fp,"bc_ylo: %d\n",d_physics_fluid.bc[2]);
        fprintf(fp,"bc_yhi: %d\n",d_physics_fluid.bc[3]);
        fprintf(fp,"bc_zlo: %d\n",d_physics_fluid.bc[4]);
        fprintf(fp,"bc_zhi: %d\n",d_physics_fluid.bc[5]);
        fprintf(fp,"alpha: %f\n",d_physics_fluid.alpha);
        fprintf(fp,"beta: %f\n",d_physics_fluid.beta);
        fprintf(fp,"cfl: %f\n",d_physics_fluid.cfl);
        fprintf(fp,"gamma: %f\n",d_physics_fluid.gamma);
        fprintf(fp,"density: %f\n",d_physics_fluid.density);
        fprintf(fp,"pressure: %f\n",d_physics_fluid.pressure);
        fprintf(fp,"mach: %f\n",d_physics_fluid.mach);
        fprintf(fp,"</FLUID PROPERTIES>\n\n");
        
        fclose(fp);
        free(backup_name);
        
    }
    
    
    /* fix inputs */
    const int MAX_INT = 2147483647;
    if(d_physics_simulation.checkpoint_interval<=0)d_physics_simulation.checkpoint_interval = MAX_INT;
    if(d_physics_simulation.visualization_interval<=0)d_physics_simulation.visualization_interval = MAX_INT;
    if(d_physics_simulation.regrid_interval<=0)d_physics_simulation.regrid_interval = MAX_INT;
    
    
    if(noinput){
        /* exit program when no input file is provided */
        exit(0);
    }
    
}
