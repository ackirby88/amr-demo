/**
 * \file   amr_utilities.c
 * \ingroup amr_group
 * \author akirby
 * 
 * \brief Utility functions for the AMR code module.
 * 
 * Created on August 16, 2018, 10:10 AM
 */

/* header files */
#include "amr_utilities.h"


void utilities_write_amr_message(int sim_dimension, int mpi_rank){
    
    char *svn_version = AMR_SVN_REVISION;
    svn_version++;
    if(mpi_rank==0){
     
        printf("\n");
        printf("+==========================================+\n");
        printf("|              HPC AMR Driver              |\n");
        printf("|                                          |\n");
        printf("|              Build ID: %d  %14s|\n",     AMR_BUILD_NUMBER,"");
        printf("|           SVN %s  %11s|\n",              svn_version,"");
        printf("+==========================================+\n");
        printf("[ amr ] Simulation Dimension: %d\n",sim_dimension);
        
    }

}


void utilities_inputs_amr_message(ctx_t *d_ctx, int mpi_rank){
    
    if(mpi_rank==0){
     
        printf("\n");
        printf("+==========================================+\n");
        printf(" AMR Inputs:                                \n");
        printf(" ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾                        \n");
        
        printf(" restart_flag: %d\n",d_ctx->d_initialize.restart_flag);
        if(d_ctx->d_initialize.restart_flag){
            printf(" restart_file: %s\n",d_ctx->d_initialize.restart_file);
        }
        
        printf(" unstructured_flag: %d\n",d_ctx->d_grid.unstructured_flag);
        if(d_ctx->d_grid.unstructured_flag){
            printf(" unstructured_file: %s\n",d_ctx->d_grid.unst_grid_file);
        }
        
        printf(" periodic_flag: %d %d %d\n",d_ctx->d_grid.periodic[0],d_ctx->d_grid.periodic[1],d_ctx->d_grid.periodic[2]);
        printf(" domain_lo: %f %f %f\n",d_ctx->d_grid.xlo[0],d_ctx->d_grid.xlo[1],d_ctx->d_grid.xlo[2]);
        printf(" domain_hi: %f %f %f\n",d_ctx->d_grid.xhi[0],d_ctx->d_grid.xhi[1],d_ctx->d_grid.xhi[2]);
        printf(" min_dx: %f\n",d_ctx->d_grid.min_dx);
        printf(" max_amr_level: %d\n",d_ctx->d_grid.max_level);
        printf(" min_amr_level: %d\n",d_ctx->d_grid.min_level);
        printf("\n+==========================================+\n");
        
    }
    
    
}


void utilities_restart_amr_message(char *filename,int index,int mpi_rank){
    
    if(mpi_rank==0){
     
        printf("\n");
        printf("+==========================================+\n");
        printf(" Restarting Solution:                       \n");
        printf(" ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾                        \n");
        printf(" file: %s\n",filename);
        printf(" Index: %d\n",index);
        printf("+==========================================+\n\n");
        
    }
    
}


void utilities_write_amr_regrid(int mpi_rank,double max_time,
                                double balance_time,double partition_time,
                                double max_point_time,double max_box_time,
                                double max_feature_time){
    
    if(mpi_rank==0){
        
        printf("[regrid] Regrid Time (sec): %1.2e\n",max_time);
        printf("         ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n");
        printf("         •Balance:          %1.2e\n",balance_time);
        printf("         •Partition:        %1.2e\n",partition_time);
        printf("         ---------------------------\n");
        printf("         •Regrid Point:     %1.2e\n",max_point_time);
        printf("         •Regrid Box:       %1.2e\n",max_box_time);
        printf("         •Regrid Feature:   %1.2e\n",max_feature_time);
        printf("         ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n");
        
    }
    
}


void utilities_write_amr_final(int mpi_rank,double total_time,
                                double compute_time,double residual_time){
    
    if(mpi_rank==0){
     
        printf("\n");
        printf("+=============================================+\n");
        printf(" Simulation Complete:                       \n");
        printf(" ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾                        \n");
        printf(" Run Time      (Compute  + Commun): %f\n", total_time);
        printf(" Compute Time  (Residual + Update): %f\n", compute_time);
        printf(" Residual Time (Spatial  Residual): %f\n", residual_time);
        printf(" Communication Percent:             %f\n",(total_time-compute_time)/total_time * 100.0);
        printf("+=============================================+\n");
        
    }
    
}


void utilities_create_directories(int mpi_rank,sc_MPI_Comm mpicomm){
    
    if(mpi_rank==0){
        
        printf("[ amr ] creating WRK directories\n");
        
#ifndef _WINDOWS
        mkdir("WRK", S_IRWXU | S_IRWXG | S_IRWXO);
        mkdir("WRK/solution", S_IRWXU | S_IRWXG | S_IRWXO);
        mkdir("WRK/solution/volume", S_IRWXU | S_IRWXG | S_IRWXO);
        mkdir("WRK/checkpoint", S_IRWXU | S_IRWXG | S_IRWXO);
#else
        _mkdir("WRK");
        _mkdir("WRK/solution");
        _mkdir("WRK/solution/volume");
        _mkdir("WRK/checkpoint");
#endif
        
    }
    sc_MPI_Barrier(mpicomm);
  
}


int utilities_find_keyword_double(char* filename, char *keyword, double *dbl){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            *dbl = atof(&buff[length]);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_integer(char* filename, char *keyword, int *integer){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            *integer = atoi(&buff[length]);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_two_integers(char* filename, char *keyword, 
                                        int *integer1, 
                                        int *integer2){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
          sscanf(&buff[length],"%d %d", integer1,integer2);
          sscanf(&buff[length],"%d, %d",integer1,integer2);
          
          fclose(fp);
          return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_three_integers(char* filename, char *keyword, 
                                            int *integer1, 
                                            int *integer2, 
                                            int *integer3){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
            sscanf(&buff[length],"%d %d %d",  integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d, %d",integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d %d", integer1,integer2,integer3);
            sscanf(&buff[length],"%d %d, %d", integer1,integer2,integer3);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_three_doubles(char* filename, char *keyword, 
                                            double *dbl1, 
                                            double *dbl2, 
                                            double *dbl3){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
            sscanf(&buff[length],"%le %le %le",  dbl1,dbl2,dbl3);
            sscanf(&buff[length],"%le, %le, %le",dbl1,dbl2,dbl3);
            sscanf(&buff[length],"%le %le, %le", dbl1,dbl2,dbl3);
            sscanf(&buff[length],"%le, %le %le", dbl1,dbl2,dbl3);
            
            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


int utilities_find_keyword_string(char *filename, char *keyword, 
                                    char *string){
    FILE *fp;
    char buff[BUFF_SIZE];
    unsigned long i,length;

    fp = fopen(filename,"r");
    length = strlen(keyword);
  
    while(fp!=NULL && fgets(buff,sizeof(buff),fp) !=NULL){
        if(strstr(buff,keyword)){
      
            int ind=0;
            for(i=length;i<BUFF_SIZE;i++){
                if(buff[i] != ' '){
                    string[ind] = buff[i];
                    ind = ind+1;
                }
            }

            string[strcspn(string, "\n")] = 0;

            if(strlen(string) >= BUFF_SIZE){
                printf("[ AMR ] string length is greater than allowable");
                exit(0);
            }

            fclose(fp);
            return 0;
        }
    }
  
    fclose(fp);
    return 1;
}


void utilities_timer(double *time_out){
    *time_out = MPI_Wtime();
}


double utilities_mpireducemax_double(sc_MPI_Comm mpicomm,double double_in){
    
    double max_double;
    
    MPI_Reduce(&double_in,&max_double,1,MPI_DOUBLE,MPI_MAX,0,mpicomm);
    return max_double;
    
}
