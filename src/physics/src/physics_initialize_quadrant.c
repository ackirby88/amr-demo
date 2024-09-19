/**
 * \file   physics_initialize_quadrant.c
 * \ingroup physics_group
 * \author akirby
 *
 * \brief Initialize quadrant solution data external solver interface.
 *
 * Created on August 28, 2018, 4:17 PM
 */

/* header files */
#include "physics_initialize_quadrant.h"

void physics_initialize_quadrant(int *dim,int *unstructured,int *nfields,
                                 int *nsub_elem,int *cell_type,double *elem_h,
                                 double *elem_xyz,double *elem_soln){

    solver_initialize_quadrant(&d_physics_fluid.initial_condition,dim,unstructured,
                                nfields,nsub_elem,cell_type,elem_h,elem_xyz,elem_soln,
                                &d_physics_fluid.mach,&d_physics_fluid.alpha,
                                &d_physics_fluid.beta,&d_physics_fluid.gamma,
                                &d_physics_fluid.density,&d_physics_fluid.pressure);
}