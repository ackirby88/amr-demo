/** MPI wrapper functions for AMR code module.
 * 
 * \file   amr_mpi.c
 * \ingroup amr_group
 * \author akirby
 * 
 * \brief Wrapper functions around mpi and p4est initialization and finalize.
 *
 * Created on August 14, 2018, 4:55 PM
 */

/* header files */
#include "amr_mpi.h"


int mpi_init(int argc,char** argv,int *mpi_rank,sc_MPI_Comm *comm){
    
    int mpi_return;
    
    mpi_return = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpi_return);
    *comm = sc_MPI_COMM_WORLD;
    
    sc_init(*comm,1,1,NULL,SC_LP_ALWAYS);
    p4est_init(NULL,SC_LP_ALWAYS);
    
    sc_MPI_Comm_rank(*comm,mpi_rank);
    //printf("[ AMR ] My mpi rank: %d\n",*mpi_rank); 
    
    return(mpi_return);
}


int mpi_finalize(){
    
    int mpi_return;
    
    sc_finalize();
    mpi_return = sc_MPI_Finalize();
    SC_CHECK_MPI(mpi_return);
    
    return(mpi_return);
}

