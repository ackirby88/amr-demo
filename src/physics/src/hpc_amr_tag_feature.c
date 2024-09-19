/**
 * \file   hpc_amr_tag_feature.c
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief AMR-Physics interface function for user-specied 
 *       feature-based AMR tagging functions.
 *
 * Created on September 7, 2018, 5:06 PM
 */


/* header files */
#include "physics_tag_feature.h"


void hpc_amr_tag_feature(int *dim,int *unstructured,int *numfields,int *nsub_elem,int *cell_type,double *elem_h,double *elem_xyz,double *elem_soln,int *tag){

    physics_tag_feature(&d_physics_fluid.amr_ntag_methods,
                         d_physics_fluid.amr_tag_methods,
                         d_physics_fluid.amr_tag_tolerances,
                         dim,unstructured,numfields,nsub_elem,
                         cell_type,elem_h,elem_xyz,
                         elem_soln,tag);
    
}
