/**
 * \file   hpc_amr_refine_operator_2d.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief AMR-Physics interface function for 2D refine operator.
 * 
 * Created on September 6, 2018, 11:19 AM
 */


/* header files */
#include "physics_refine_operator.h"


void hpc_amr_refine_operator_2d(int *unstructured,int *numfields,
                                int *nsub_elem,int *cell_type,double *elem_h,double *geom,
                                double *soln,
                                double *uc_1,double *uc_2,double *uc_3,double *uc_4){

    physics_refine_operator_2d(unstructured,numfields,
                               nsub_elem,cell_type,elem_h,geom,
                               soln,
                               uc_1,uc_2,uc_3,uc_4);
    
}
