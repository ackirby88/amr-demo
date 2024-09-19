/**
 * \file   physics_runtime_statistics.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Displays the run-time performance statistics for the numeric 
 *        discretization including information related to parallel process 
 *        wall-clock times.
 * 
 * Created on August 23, 2018, 10:36 AM
 */

/* header files */
#include "physics_runtime_statistics.h"


void runtime_statistics(int commrank,int mpisize,int steps,int nfields,int dof,
                        int nresiduals,int extern_cell_size,double t1,double t2,
                        double residual_time,double compute_time,
                        MPI_Comm mpicomm_ext){
    
    char percent[2];
    strcpy(percent,"%");
    
    double min_time_per_dof;
    double max_time_per_dof;
    double avg_time_per_dof;
    double min_time_per_dof_compute;
    double max_time_per_dof_compute;
    double avg_time_per_dof_compute;
    double min_time_per_dof_residual;
    double max_time_per_dof_residual;
    double avg_time_per_dof_residual;
    
    double time_per_dof = 0.0;
    double time_per_dof_compute = 0.0;
    double time_per_dof_residual = 0.0;
    
    double avg_time_residual;
    double avg_time_compute;
    double std_dev_compute;
    double diff_time_compute;
    
    data_rank_t residual_data;
    data_rank_t compute_data;
    data_rank_t min_time_residual;
    data_rank_t max_time_residual;
    data_rank_t min_time_compute;
    data_rank_t max_time_compute;
    

    
    residual_data.time = residual_time;
    residual_data.mpi_rank = commrank;
    
    compute_data.time = compute_time;
    compute_data.mpi_rank = commrank;
    
    // max,min,avg residual time
    MPI_Reduce(&residual_data, &min_time_residual, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, mpicomm_ext);
    MPI_Reduce(&residual_data, &max_time_residual, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, mpicomm_ext);
    MPI_Reduce(&residual_time, &avg_time_residual, 1, MPI_DOUBLE, MPI_SUM,        0, mpicomm_ext);
    avg_time_residual /= (double) mpisize;
    
    
    // max,min,avg compute time
    MPI_Allreduce(&compute_data, &max_time_compute, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpicomm_ext);
    MPI_Allreduce(&compute_time, &avg_time_compute, 1, MPI_DOUBLE, MPI_SUM, mpicomm_ext);
    
    MPI_Reduce(&compute_data, &min_time_compute, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, mpicomm_ext);
    MPI_Reduce(&compute_data, &max_time_compute, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, mpicomm_ext);
    MPI_Reduce(&compute_time, &avg_time_compute, 1, MPI_DOUBLE, MPI_SUM,        0, mpicomm_ext);
    avg_time_compute /= (double) mpisize;
    
    
    // standard deviation of compute times
    diff_time_compute = (compute_time - avg_time_compute)*(compute_time - avg_time_compute);
    MPI_Reduce(&diff_time_compute, &std_dev_compute, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm_ext);
    std_dev_compute = sqrt(std_dev_compute/(double) mpisize);
    
    
    if(dof){
        time_per_dof        = ((double) nfields)*(t2-t1)/( (double) dof)/((double) steps)/(double) nresiduals;
        time_per_dof_compute = ((double) nfields)*compute_time/( (double) dof)/((double) steps)/(double) nresiduals;
        time_per_dof_residual = ((double) nfields)*residual_time/( (double) dof)/((double) steps)/(double) nresiduals;
    }
  
    MPI_Reduce(&time_per_dof, &min_time_per_dof, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm_ext);
    MPI_Reduce(&time_per_dof, &max_time_per_dof, 1, MPI_DOUBLE, MPI_MAX, 0, mpicomm_ext);
    MPI_Reduce(&time_per_dof, &avg_time_per_dof, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm_ext);

    MPI_Reduce(&time_per_dof_compute, &min_time_per_dof_compute, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm_ext);
    MPI_Reduce(&time_per_dof_compute, &max_time_per_dof_compute, 1, MPI_DOUBLE, MPI_MAX, 0, mpicomm_ext);
    MPI_Reduce(&time_per_dof_compute, &avg_time_per_dof_compute, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm_ext);
    
    MPI_Reduce(&time_per_dof_residual, &min_time_per_dof_residual, 1, MPI_DOUBLE, MPI_MIN, 0, mpicomm_ext);
    MPI_Reduce(&time_per_dof_residual, &max_time_per_dof_residual, 1, MPI_DOUBLE, MPI_MAX, 0, mpicomm_ext);
    MPI_Reduce(&time_per_dof_residual, &avg_time_per_dof_residual, 1, MPI_DOUBLE, MPI_SUM, 0, mpicomm_ext);

    avg_time_per_dof /= (double) mpisize;
    avg_time_per_dof_compute /= (double) mpisize;
    avg_time_per_dof_residual /= (double) mpisize;
    
    
    // Output statistics
    if(commrank==0){
        
        printf("[physics] time/dof          (min,max,avg):   %e,   %e, %e\n",min_time_per_dof,max_time_per_dof,avg_time_per_dof);
        printf("[physics] time/dof compute  (min,max,avg):   %e,   %e, %e\n",
          min_time_per_dof_compute,max_time_per_dof_compute,avg_time_per_dof_compute);
        printf("[physics] time/dof residual (min,max,avg):   %e,   %e, %e\n",
          min_time_per_dof_residual,max_time_per_dof_residual,avg_time_per_dof_residual);
        printf("[physics] compute  times    (min,max,avg):[%d]%e,[%d]%e, %e,(std) %1.0e:%3.0f%s\n",
          min_time_compute.mpi_rank,min_time_compute.time,
          max_time_compute.mpi_rank,max_time_compute.time,
          avg_time_compute,std_dev_compute,std_dev_compute/avg_time_compute*100.0,percent);
        printf("[physics] residual times    (min,max,avg):[%d]%e,[%d]%e, %e\n",
          min_time_residual.mpi_rank,min_time_residual.time,
          max_time_residual.mpi_rank,max_time_residual.time,
          avg_time_residual);
        printf("\n");
        
    }
    
    
    //slowest processor details
    //if(commrank == max_time_compute.mpi_rank){
    //    for(p=0;p<ntypes;p++){
    //        if(cells_per_cell_type[p]){
    //            printf("[Slow#] [%d]: cell_type[%d], cells[%d]\n",
    //              max_time_compute.mpi_rank,p,cells_per_cell_type[p]);
    //        }
    //    }
    //}
    
}
