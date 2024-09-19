/**
 * \file   amr_utilities.h 
 * \author akirby
 *
 * Created on August 14, 2018, 4:08 PM
 */

#ifndef AMR_UTILITIES_H
#define AMR_UTILITIES_H


/* header files */
#include "var_defines.h"
#include "var_ctx.h"
#include "amr_version_tag.h"

#include <mpi.h>

#ifndef DOXYGEN_IGNORE
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef WINDOWS
#include<windows.h>
#else
#include <sys/stat.h>
#endif
#endif

/** Display AMR module message to terminal
 * 
 * @param [in] sim_dimension     simulation dimension (2D/3D)
 * @param [in] mpi_rank          mpi rank
 */
void utilities_write_amr_message(int sim_dimension,int mpi_rank);

/** Display AMR input options from file to terminal
 * 
 * @param [in] d_ctx        context data related to simulation
 * @param [in] mpi_rank     mpi rank
 */
void utilities_inputs_amr_message(ctx_t *d_ctx, int mpi_rank);

/** Display AMR restart simulation message to terminal
 * 
 * @param [in] filename     restart file name read from inputs
 * @param [in] index        file index number that was extracted from file name
 * @param [in] mpi_rank     mpi rank
 */
void utilities_restart_amr_message(char *filename,int index,int mpi_rank);

/** Display AMR regrid timings
 * 
 * @param [in] mpi_rank         mpi rank
 * @param [in] max_time         total regrid time
 * @param [in] balance_time     total amr 2:1 balancing time
 * @param [in] partition_time   total mpi partitioning time
 * @param [in] max_point_time   total point regridding time
 * @param [in] max_box_time     total box regridding time
 * @param [in] max_feature_time total feature-based regridding time
 */
void utilities_write_amr_regrid(int mpi_rank,double max_time,double balance_time,double partition_time,double max_point_time,double max_box_time,double max_feature_time);

/** Display final AMR module message to terminal
 * 
 * @param [in] mpi_rank         mpi rank
 * @param [in] total_time       total simulation time
 * @param [in] compute_time     total computation time
 * @param [in] residual_time    total residual time
 */
void utilities_write_amr_final(int mpi_rank,double total_time,double compute_time,double residual_time);

/** Build the output directories for the simulation
 * 
 * @param [in] mpi_rank     mpi rank
 * @param [in] mpicomm      mpi communicator
 */
void utilities_create_directories(int mpi_rank,sc_MPI_Comm mpicomm);

/** Input file read line helper function: one integer
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer           integer variable to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_integer(char* filename,char *keyword,int *integer);

/** Input file read line helper function: two integers
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer1          integer variable #1 to read
 * @param [out] integer2          integer variable #2 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_two_integers(char* filename,char *keyword,int *integer1,int *integer2);

/** Input file read line helper function: three integers
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] integer1          integer variable #1 to read
 * @param [out] integer2          integer variable #2 to read
 * @param [out] integer3          integer variable #3 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_three_integers(char* filename,char *keyword,int *integer1,int *integer2,int *integer3);

/** Input file read line helper function: one double
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] dbl         double precision variable to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_double(char* filename,char *keyword,double *dbl);

/** Input file read line helper function: three doubles
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] dbl1        double precision variable #1 to read
 * @param [out] dbl2        double precision variable #2 to read
 * @param [out] dbl3        double precision variable #3 to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_three_doubles(char* filename,char *keyword,double *dbl1,double *dbl2,double *dbl3);

/** Input file read line helper function: string
 * 
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] string      string to read
 * @return                  returns 0 if successfully read variable
 */
int utilities_find_keyword_string(char *filename,char *keyword,char *string);

/** MPI wallclock timer wrapper function
 * 
 * @param [out] time_out  mpi wall time
 */
void utilities_timer(double *time_out);

/** MPI reduce double max function wrapper
 * 
 * @param [in] mpicomm      mpi communicator     
 * @param [in] double_in    double precision number to find global max value
 * @return  global max double value 
 */
double utilities_mpireducemax_double(sc_MPI_Comm mpicomm,double double_in);


#endif /* AMR_UTILITIES_H */
