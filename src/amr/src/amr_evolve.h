/** 
 * \file   amr_evolve.h
 * \author akirby
 *
 * Created on August 20, 2018, 10:40 AM
 */

#ifndef AMR_EVOLVE_H
#define AMR_EVOLVE_H


/* header files */
#include "amr_external.h"
#include "amr_regrid.h"
#include "amr_vtk.h"
#include "var_ctx.h"


/** AMR stepping kernel function:
 * 
 * Primary evolve function for AMR code module:
 *        Evolve/Visualize/Checkpoint/Regrid solution.
 * 
 * @param [in] d_ctx        context data related to each quadrant
 * @param [in,out] p4est    p4est tree data structure
 * @param [in,out] conn     p4est tree connectivity
 */
void evolve(ctx_t *d_ctx,p4est_t **p4est,p4est_connectivity_t **conn);


#endif /* AMR_EVOLVE_H */

