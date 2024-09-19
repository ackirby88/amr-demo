/** 
 * \file   physics_tag_feature.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Physics tag feature interface to external solver tag feature.
 * 
 * Created on September 7, 2018, 5:07 PM
 */


/* header files */
#include "physics_tag_feature.h"


void physics_tag_feature(int *ntag,int *tag_methods, double *tag_tolerances,
                         int *dim,int *unstructured,int *nfields,int *nsub_elem,
                         int *cell_type,double *elem_h,double *elem_xyz,
                         double *soln,int *tag){
    
    solver_tag_feature(ntag,tag_methods,tag_tolerances,
                       dim,unstructured,nfields,nsub_elem,
                       cell_type,elem_h,elem_xyz,
                       soln,tag);
    
}
