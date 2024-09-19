/**
 * \file   amr_initialize.h
 * \author akirby
 * 
 * Created on August 15, 2018, 1:55 PM
 */

#ifndef AMR_INITIALIZE_H
#define AMR_INITIALIZE_H


/* header files */
#include "var_ctx.h"
#include "var_quad.h"
#include "main.h"
#include "amr_regrid.h"
#include "amr_external.h"
#include "amr_utilities.h"
#include "amr_p4est_utilities.h"

#ifndef DOXYGEN_IGNORE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#endif

#ifndef P4_TO_P8
#include <p4est_extended.h>
#else
#include <p8est_extended.h>
#endif

/** Initializes the AMR code module
 * 
 * Builds output directories, reads input files, builds the mesh, and 
 * initializes the solution data.
 * 
 * @param [in]  argc    number of command line arguments
 * @param [in]  argv    command line arguments
 * @param [out] d_ctx   simulation context data
 * @param [out] p4est   p4est tree structure
 * @param [out] conn    p4est tree connectivity
 */
void initialize(int argc,char** argv,ctx_t *d_ctx,p4est_t **p4est,p4est_connectivity_t **conn);

/** Reads the input file from the command line.
 * 
 * @param [in]  filename    input file name
 * @param [out] d_ctx       simulation context data
 * @param [in]  noinput     flag if input file was found
 */
void initialize_inputs(char *filename,ctx_t *d_ctx,int noinput);

/** Saves the inputs read from the input file to WRK/saved.(filename)
 * 
 * @param [in]  filename    input file name
 * @param [in]  d_ctx       simulation context data
 * @param [in]  noinput     flag if input file was found
 */
void initialize_save_inputs(char *filename,ctx_t *d_ctx,int noinput);

/** Builds new grid or reads in restart file. Initializes external solver data.
 * 
 * @param [in]  d_ctx   simulation context data
 * @param [out] p4est   p4est tree structure
 * @param [out] conn    p4est tree connectivity
 */
void initialize_solver(ctx_t *d_ctx,p4est_t **p4est,p4est_connectivity_t **conn);

/** Builds new Cartesian grid given user dimensions.
 * 
 * @param [in,out] d_ctx   simulation context data
 * @param [out] conn    p4est tree connectivity
 */
void initialize_grid(ctx_t *d_ctx,p4est_connectivity_t **conn);

/** Callback function to initialize solution data in each amr quadrant/octant.
 * 
 * @param [in]  p4est       p4est tree structure
 * @param [in]  which_tree  tree in the forest of octrees
 * @param [out] q           quadrant solution data
 */
void initialize_quadrant_data(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *q);


#endif /* AMR_INITIALIZE_H */

