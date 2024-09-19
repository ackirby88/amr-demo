/**
 * \file   hpc_amr_initialize_quadrant.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief AMR-Physics interface function for initializing a quadrant with
 *        solution data.
 *
 * Created on August 21, 2018, 5:05 PM
 */

/* header files */
#include "physics_initialize_quadrant.h"

void hpc_amr_initialize_quadrant(int *dim,int *unstructured,int *numfields,
                                 int *nsub_elem,int *cell_type,double *elem_h,
                                 double *elem_xyz,double *elem_soln){

    physics_initialize_quadrant(dim,unstructured,numfields,nsub_elem,
                                cell_type,elem_h,elem_xyz,elem_soln);
}