/**
 * \file   physics_runtime_statistics.h
 * \author akirby
 *
 * Created on August 23, 2018, 10:34 AM
 */

#ifndef PHYSICS_RUNTIME_STATISTICS_H
#define PHYSICS_RUNTIME_STATISTICS_H

/* header files */
#ifndef DOXYGEN_IGNORE
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif


/** Diplays the run-time statistics related to simulation communication & computation.
 * 
 * @param [in] commrank         MPI rank
 * @param [in] mpisize          Total MPI ranks
 * @param [in] steps            Number of time steps
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] dof              Number of degrees of freedom
 * @param [in] nresiduals       Number of residuals calculated per step
 * @param [in] extern_cell_size Number of constants in the cell_info structure
 * @param [in] t1               Wall-clock start time (total)
 * @param [in] t2               Wall-clock stop time (total)
 * @param [in] residual_time    Wall-clock time for residual calculation
 * @param [in] compute_time     Wall-clock time for residual and data update
 * @param [in] mpicomm_ext      MPI communicator
 */
void runtime_statistics(int commrank,int mpisize,int steps,int nfields,int dof,
                        int nresiduals,int extern_cell_size,double t1,double t2,
                        double residual_time,double compute_time,
                        MPI_Comm mpicomm_ext);


/**
 * data_rank_t is a struct containing time and mpi_rank that used to determine 
 * which rank had what wall-clock time when using mpi_allreduce functions
 */
typedef struct {
    double time;    /**< Wall-clock time for this MPI rank */
    int mpi_rank;   /**< MPI rank*/
}
data_rank_t;


#endif /* PHYSICS_RUNTIME_STATISTICS_H */

