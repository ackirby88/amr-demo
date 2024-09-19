/**
 * \file   amr_tag.h
 * \author akirby
 *
 * Created on September 7, 2018, 9:18 AM
 */

#ifndef AMR_TAG_H
#define AMR_TAG_H


/* header files */
#include "var_ctx.h"
#include "var_tag.h"
#include "var_quad.h"
#include "amr_external.h"
#include "amr_p4est_utilities.h"


#ifdef __cplusplus
extern "C" {
#endif



#ifdef __cplusplus
}
#endif


/** Resets all refinement tags to no_tag
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] tag_default  tag number to set all quadrant tags
 * @return flag for refinement of quadrant (always 0)
 */
int tag_reset(p4est_t *p4est,int tag_default);

/** Tags all quadrants for refinement
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_all(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant);

/** Tags quadrant for refinement if found in or overlaps user-specified box region
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_box(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant);

/** Tags quadrant for refinement if user-specified point is found in quadrant
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_point(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant);

/** Tags quadrant for refinement if user-specified feature found in quadrant
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @return flag for refinement of quadrant
 */
int tag_feature(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant);

/** Tags quadrant for coarsening if no refinement tags found
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] children     quadrant children tagged for coarsening
 * @return flag for coarsening quadrant
 */
int tag_coarsen(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *children[]);

/** Point search in p4est forest for list of points
 * 
 * @param [in] p4est        p4est forest structure
 * @param [in] which_tree   tree id in p4est forest
 * @param [in] quadrant     quadrant in p4est tree
 * @param local_num         leaf id number 
 * @param points_in         user-specified point list
 * @return flag for point found in quadrant
 */
int tag_point_search(p4est_t *p4est,p4est_topidx_t which_tree,p4est_quadrant_t *quadrant,p4est_locidx_t local_num,void *points_in);

/** Overlap check function for two line segments
 * 
 * @param x1    lower coordinate of line #1
 * @param x2    upper coordinate of line #1
 * @param xx1   lower coordinate of line #2
 * @param xx2   upper coordinate of line #2
 * @return flag for overlap found
 */
int check_overlap(double x1,double x2,double xx1,double xx2);


#endif /* AMR_TAG_H */

