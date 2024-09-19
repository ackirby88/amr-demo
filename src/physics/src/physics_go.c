/**
 * \file   physics_go.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Physics CPU go function that sets up mpi communication 
 *        and calls solve.
 * 
 * Created on August 22, 2018, 2:00 PM
 */

/* header files */
#include "physics_go.h"


void physics_go(int* dim,
                int* unstructured,
                int* nfields,
                int* nelem_subgrid,
                int* ncell,
                int* nghost,
                int* nface,
                int* dof,
                int* cell_info,
                int* face_info,
                int* nsend,
                int* nrecv,
                int* send,
                int* recv,
                int* sizes,
                double* geom,
                double* soln,
                int* istep,
                int* regrid,
                int* visualize,
                int* checkpoint,
                int* evolve_solution,
                MPI_Comm mpicomm_ext,
                double* total_time,
                double* compute_time,
                double* residual_time){
    
    int i;
    int cell;
    
    int min_dof;
    int max_dof;
    int avg_dof;
    int total_dof;
    int dof_loc = *dof/d_physics_interface.nfields;
  
    double t1,t2;
    double compute_time_loc;
    double residual_time_loc;
    double residual_norm;

    
    int commrank;
    int mpisize;
    
    MPI_Comm_rank(mpicomm_ext,&commrank);
    MPI_Comm_size(mpicomm_ext,&mpisize);
    
    
    /* exit simulation conditional */
    if(*istep >= d_physics_simulation.time_steps){
        *evolve_solution = 0;
        return;
    }
    
    /* set up MPI biffers before solving */
    //========================================================================//
    int send_buf_size = 0;
    for(i=0;i<*nsend;++i){
        cell = send[d_physics_interface.extern_mpi_comm_size*i+0];
        send_buf_size += sizes[cell];
    }
    
    int* counts          = (int *)         malloc(sizeof(int)*(*ncell+*nghost));
    double* sendbuf      = (double *)      malloc(sizeof(double)*send_buf_size);
    MPI_Request* request = (MPI_Request *) malloc(sizeof(MPI_Request)*(*nrecv+*nsend));
    
    counts[0] = 0;
    for(i=1;i<*ncell+*nghost;++i){
        counts[i] = counts[i-1]+sizes[i-1];
    }
    //========================================================================//
    
    
    /* solve */
    physics_solve_interface(dim,unstructured,nfields,nelem_subgrid,ncell,nghost,nface,dof,
                            cell_info,face_info,nsend,nrecv,send,recv,
                            sizes,geom,soln,istep,regrid,visualize,checkpoint,
                            mpicomm_ext,counts,sendbuf,request,
                            &i,&t1,&t2,&compute_time_loc,&residual_time_loc,
                            &residual_norm);
    
    
    /* clean up MPI buffers after solving */
    physics_mpi_free_buffers(sendbuf,request,counts);
    
    
    /* display time step information */
    //========================================================================//
    MPI_Reduce(&dof_loc,   &min_dof, 1, MPI_INT, MPI_MIN, 0, mpicomm_ext);
    MPI_Reduce(&dof_loc,   &max_dof, 1, MPI_INT, MPI_MAX, 0, mpicomm_ext);
    MPI_Reduce(&dof_loc, &total_dof, 1, MPI_INT, MPI_SUM, 0, mpicomm_ext);
    avg_dof = (int) ((double) total_dof / (double) mpisize);
    
    residual_norm /= pow(total_dof,1.0/(double) *dim);
    
    if(commrank==0){
        printf("[ amr ] istep: %d\n",*istep);
        printf("[ amr ] processor dofs (min,max,avg,total): %d, %d, %d, %d \n",
                min_dof,max_dof,avg_dof,total_dof);
        printf("[physics] iterations: %d, cpu time: %e, residual: %1.15e\n",
                i,t2-t1,residual_norm);
    }
    
    *total_time += t2-t1;
    *compute_time += compute_time_loc;
    *residual_time += residual_time_loc;
    
    runtime_statistics( commrank,mpisize,i,d_physics_interface.nfields,*dof,
                        d_physics_simulation.time_scheme,
                        d_physics_interface.extern_cell_size,
                        t1,t2,residual_time_loc,compute_time_loc,mpicomm_ext);
    //========================================================================//
    
    
    
}
