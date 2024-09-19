/**
 * \file    physics_mpi.c
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief MPI solution data setup and exchange for the PHYSICS code module.
 *
 * Created on August 22, 2018, 3:12 PM
 */

/* header files */
#include "physics_mpi.h"


void physics_mpi_initialize_buffers(double* sendbuf,MPI_Request* request, int* counts,
                                    int ncell,int nghost,int nsend,int nrecv,int* send,int* sizes){
    
    int i;
    
    // add up total size needed for buffer
    int send_buf_size = 0;
    for(i=0;i<nsend;++i){
        const int cell = send[d_physics_interface.extern_mpi_comm_size*i+0];
        send_buf_size += sizes[cell];
    }
  
    // for speed could allocate this stuff once per time step instead of each ghost exchange
    sendbuf = (double *) malloc(sizeof(double)*send_buf_size);
    request = (MPI_Request *) malloc(sizeof(MPI_Request)*(nrecv+nsend));//over allocate for now
    counts = (int *) malloc(sizeof(int)*(ncell+nghost));
  
    // needed to access solution or could use cell_info already passed to external solver
    counts[0] = 0; //reset for accumulation
    for(i=1;i<ncell+nghost;++i){
        counts[i] = counts[i-1]+sizes[i-1];
    }
  
}

void physics_mpi_free_buffers(double* sendbuf,MPI_Request* request, int* counts){
    free(sendbuf);
    free(request);
    free(counts);
}


void physics_mpi_real_to_ghost(int nsend,int nrecv,int dof,int* send,int* recv,
                               int* sizes, double* vec,MPI_Comm mpicomm_ext,
                               double* sendbuf,MPI_Request* request, int* counts){
  
    int i,j;
  
    
    // ghosts are always at the end and stacked in order of processor so can directly insert them into vec
    int ind_prev    = dof;              // beginning of ghost storage (or end of real storage+1), this is passed to external solver
    int count       = sizes[recv[0]];   // size of data for first receiving cell
    int proc_prev   = recv[1];          // receiving processor id
    int nrequest    = 0;                // initialize number of requests
  
    // post non-blocking receive mpi calls
    for(i=1;i<nrecv;++i){
        const int cell = recv[d_physics_interface.extern_mpi_comm_size*i+0];
        const int proc_next = recv[d_physics_interface.extern_mpi_comm_size*i+1];
        
        if(proc_prev == proc_next){
            count += sizes[cell];
        } else {
            MPI_Irecv(&vec[ind_prev],count,MPI_DOUBLE,proc_prev,MPI_ANY_TAG,mpicomm_ext,&request[nrequest++]);
            ind_prev += count;
            proc_prev = proc_next;
            count = sizes[cell];
        }
    }
    //final receive
    if(nrecv) {
        MPI_Irecv(&vec[ind_prev],count,MPI_DOUBLE,proc_prev,MPI_ANY_TAG,mpicomm_ext,&request[nrequest++]);
    }
    
  
    // one cell can send to multiple ranks and send might not be contiguous fill in a send buffer
    int ind = 0;
    for(i=0;i<nsend;++i){
        const int cell = send[d_physics_interface.extern_mpi_comm_size*i+0];
        const int n = sizes[cell];
        //const int k = ctx->unst.cell_info[UNST_CELL_SIZE*cell+2];// or could use cell info or construct counts above
        const int k = counts[cell];     // need to index solution but all cells could be different sizes
        
        for(j=0;j<n;++j){
            sendbuf[ind++] = vec[k+j];  // copy over entire cell worth of data
        }
    }
    
    // send buffer is contiguous and stacked in order of processor
    const int tag = 35; // make some unique numbering
    ind_prev = 0; // using send buffer so start at 0
    proc_prev = send[1]; // processor to send to for first send cell
    count = sizes[send[0]]; // size of first send cell data
    
    // post non-blocking send mpi calls
    for(i=1;i<nsend;++i){
        const int cell = send[d_physics_interface.extern_mpi_comm_size*i+0];
        const int proc_next = send[d_physics_interface.extern_mpi_comm_size*i+1];
        if(proc_prev == proc_next){
            count += sizes[cell];
        } else {
            MPI_Isend(&sendbuf[ind_prev],count,MPI_DOUBLE,proc_prev,tag,mpicomm_ext,&request[nrequest++]);
            ind_prev += count;
            proc_prev = proc_next;
            count = sizes[cell];
        }
    }
    // final send
    if(nsend){
        MPI_Isend(&sendbuf[ind_prev],count,MPI_DOUBLE,proc_prev,tag,mpicomm_ext,&request[nrequest++]);
    }
    
    
    // wait for all send and receives to be placed
    MPI_Waitall(nrequest,request,MPI_STATUSES_IGNORE);
    
  
}

