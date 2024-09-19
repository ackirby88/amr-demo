/**
 * \file   amr_evolve.c
 * \ingroup amr_group
 * \author akirby
 *
 * \brief Primary evolve function for the AMR code module:
 *        Evolve/Visualize/Checkpoint/Regrid solution.
 *
 * Created on August 20, 2018, 10:42 AM
 */

/* header files */
#include "amr_evolve.h"

void evolve(ctx_t *d_ctx,p4est_t **p4est,p4est_connectivity_t **conn){
    int regrid;
    int initial;
    int visualize;
    int save_data;
    int checkpoint;
    int point_regrid;
    int feature_regrid;
    int evolve_solution;
    int allow_coarsening;
    char filename[37];

    d_ctx->total_time = 0.0;
    d_ctx->compute_time = 0.0;
    d_ctx->residual_time = 0.0;

    initial = 1;
    save_data = 1;
    point_regrid = 1;
    feature_regrid = 1;
    evolve_solution = 1;
    allow_coarsening = 1;

    /* adapt and visualize initial solution */
    if(d_ctx->istep==0){
        regrid_solution(d_ctx,*p4est,initial,point_regrid,feature_regrid,allow_coarsening);
        vtk_write_solution(*p4est,d_ctx->istep);
    }

    initial = 0;
    while(evolve_solution){
        /* step solution */
        external_evolve(d_ctx,p4est,&d_ctx->istep,&regrid,&visualize,&checkpoint,&evolve_solution);

        /* visualize and checkpoint at end of simulation */
        if(evolve_solution==0){
            visualize = 1;
            checkpoint = 1;
        }

        /* visualization solution */
        if(visualize){
            vtk_write_solution(*p4est,d_ctx->istep);
        }

        /* checkpoint solution */
        if(checkpoint){
            snprintf(filename,37,"WRK/checkpoint/amr_soln_%07d.bin",d_ctx->istep);
            p4est_save(filename,*p4est,save_data);
        }

        /* regrid solution */
        if(regrid){
            regrid_solution(d_ctx,*p4est,initial,point_regrid,feature_regrid,allow_coarsening);
        }
    }

    external_deallocate_soln(*p4est);
}