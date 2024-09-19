/**
 * \file   physics_coarsen_operator_2d.c
 * \ingroup physics_group
 * \author akirby
 * 
 * \brief External solver function wrapper for 2D coarsen operator
 *
 * Created on September 6, 2018, 10:26 AM
 */


/* header files */
#include "physics_coarsen_operator.h"


void physics_coarsen_operator_2d(int *unstructured,int *numfields,int *nsub_elem,
                                 int *cell_type,double *elem_h,double *geom,
                                 double *uc_1,double *uc_2,double *uc_3,double *uc_4,
                                 double *soln){
    
    solver_coarsen_operator_2d(unstructured,numfields,nsub_elem,
                               cell_type,elem_h,geom,
                               uc_1,uc_2,uc_3,uc_4,
                               soln);
    
}
