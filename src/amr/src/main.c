/** Top-level main function of the AMR framework.
 * \file main.c
 * \ingroup amr_group
 * \author akirby
 *
 * \brief Main function for the AMR code module.
 *
 * Created on August 14, 2018, 4:00 PM
 *
 */

/* header files */
#include "main.h"
#include "var_mpi.h"
#include "var_ctx.h"

#include "amr_mpi.h"
#include "amr_utilities.h"
#include "amr_initialize.h"
#include "amr_evolve.h"
#include "amr_regrid.h"

/** Main AMR driving function
 *
 * \param [in] argc     number of command line arguments
 * \param [in] argv     command line arguments
 *
 */
int main(int argc, char **argv){
    ctx_t d_ctx;
    p4est_t *p4est;
    p4est_connectivity_t *conn;
    int mpiret;

    d_ctx.d_grid.dim = P4EST_DIM;
    mpiret = mpi_init(argc,argv,&d_ctx.d_mpi.rank,&d_ctx.d_mpi.comm);
    utilities_write_amr_message(P4EST_DIM,d_ctx.d_mpi.rank);
    initialize(argc,argv,&d_ctx,&p4est,&conn);

    evolve(&d_ctx,&p4est,&conn);

    utilities_write_amr_final(d_ctx.d_mpi.rank,d_ctx.total_time,
                              d_ctx.compute_time,d_ctx.residual_time);

    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn);
    mpiret = mpi_finalize();
    return(mpiret);
}