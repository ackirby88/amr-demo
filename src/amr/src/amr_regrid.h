/**
 * \file   amr_regrid.h
 * \author akirby
 *
 * Created on September 5, 2018, 10:18 AM
 */

#ifndef AMR_REGRID_H
#define AMR_REGRID_H


/* header files */
#include "var_ctx.h"
#include "var_tag.h"
#include "amr_tag.h"
#include "amr_external.h"
#include "amr_utilities.h"
#include "amr_initialize.h"
#include "amr_p4est_utilities.h"

#ifndef P4_TO_P8
#include "p4est_search.h"
#else
#include "p8est_search.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif


#ifdef __cplusplus
}
#endif


/** Top-level regridding coordination function
 * 
 * @param [in] d_ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in] initial          flag if initial time step
 * @param [in] point_regrid     flag if regridding to user-specified list of points
 * @param [in] feature_regrid   flag if regridding to user-specified physics features
 * @param [in] allow_coarsening flag if regrid coarsening
 */
void regrid_solution(ctx_t *d_ctx,p4est_t *p4est,int initial,int point_regrid,int feature_regrid,int allow_coarsening);

/** Regrid coarsening function
 * 
 * @param [in] d_ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in]  initial         flag if initial time step
 * @param [out] coarsen_time
 */
void regrid_coarsen(ctx_t *d_ctx,p4est_t *p4est,int initial,double *coarsen_time);

/** Mesh 2:1 quadrant balancing and mpi partitioning function
 * 
 * @param [in] d_ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in]  initial         flag if initial time step
 * @param [out] bal_time        wall-clock time for balancing function
 * @param [out] part_time       wall-clock time for partitioning function
 */
void regrid_balance_partition(ctx_t *d_ctx,p4est_t *p4est,int initial,double *bal_time,double *part_time);

/** Regridding function for adapting to user-specified box region
 * 
 * @param [in] d_ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in] initial          flag if initial time step
 */
void regrid_box(ctx_t *d_ctx,p4est_t *p4est,int initial);

/** Regridding function for adapting to user-specified set of points (external)
 * 
 * @param [in] d_ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in] initial          flag if initial time step
 */
void regrid_points(ctx_t *d_ctx,p4est_t *p4est,int initial);

/** Regridding function for adapting to user-specified features (external)
 * 
 * @param [in] d_ctx            simulation context data
 * @param [in,out] p4est        p4est forest structure
 * @param [in] initial          flag if initial time step
 */
void regrid_feature(ctx_t *d_ctx,p4est_t *p4est,int initial);

/** Solution replacement function when coarsening or refining (external)
 * 
 * @param [in] p4est            p4est forest structure
 * @param [in] which_tree       tree id in p4est forest
 * @param [in] num_outgoing     number of quadrants being replaced
 * @param [out] outgoing        quadrants being replaced
 * @param [in] num_incoming     number of quadrants being introduced
 * @param [in] incoming         quadrant being introduced
 */
void regrid_replace_quads(p4est_t *p4est,p4est_topidx_t which_tree,int num_outgoing,p4est_quadrant_t *outgoing[],int num_incoming,p4est_quadrant_t *incoming[]);

/** Mpi partition weighting function per quadrant (external)
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     p4est quadrant that needs to be weighted
 * @return partition weight for this quadrant
 */
int  regrid_partition_weight(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant);


#endif /* AMR_REGRID_H */

