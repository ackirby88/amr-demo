/**
 * \file   physics_tag_feature.h
 * \author akirby
 *
 * Created on September 7, 2018, 5:09 PM
 */

#ifndef PHYSICS_TAG_FEATURE_H
#define PHYSICS_TAG_FEATURE_H

/* header files */
#include "physics_var_physics.h"

#ifdef __cplusplus
extern "C" {
#endif

/** External solver AMR tag feature function
 * 
 * @param [in] ntag             Number of feature tag methods
 * @param [in] tag_methods      Feature tagging methods flags
 * @param [in] tag_tolerances   Feature tagging tolerances for each method
 * @param [in] dim              Simulation spatial dimension
 * @param [in] unstructured     Unstructured grid flag
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nsub_elem        Number of sub-grid elements per AMR quad
 * @param [in] cell_type        Cell type
 * @param [in] elem_h           Length of AMR quad (structured)
 * @param [in] elem_xyz         AMR cell geometry coordinates
 * @param [in] soln             Cell Solution pointer
 * @param [out] tag             Cell AMR tag flag
 */
void solver_tag_feature(int *ntag,int *tag_methods, double *tag_tolerances,
                         int *dim,int *unstructured,int *nfields,int *nsub_elem,
                         int *cell_type,double *elem_h,double *elem_xyz,
                         double *soln,int *tag);

#ifdef __cplusplus
}
#endif


/** Physics tag feature interface to external solver tag feature. 
 * 
 * @param [in] ntag             Number of feature tag methods
 * @param [in] tag_methods      Feature tagging methods flags
 * @param [in] tag_tolerances   Feature tagging tolerances for each method
 * @param [in] dim              Simulation spatial dimension
 * @param [in] unstructured     Unstructured grid flag
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nsub_elem        Number of sub-grid elements per AMR quad
 * @param [in] cell_type        Cell type
 * @param [in] elem_h           Length of AMR quad (structured)
 * @param [in] elem_xyz         AMR cell geometry coordinates
 * @param [in] soln             Cell Solution pointer
 * @param [out] tag             Cell AMR tag flag
 */
void physics_tag_feature(int *ntag,int *tag_methods, double *tag_tolerances,
                         int *dim,int *unstructured,int *nfields,int *nsub_elem,
                         int *cell_type,double *elem_h,double *elem_xyz,
                         double *soln,int *tag);


#endif /* PHYSICS_TAG_FEATURE_H */

