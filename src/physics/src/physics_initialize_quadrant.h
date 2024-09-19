/**
 * \file   physics_initialize_quadrant.h
 * \author akirby
 *
 * Created on August 28, 2018, 4:19 PM
 */

#ifndef PHYSICS_INITIALIZE_QUADRANT_H
#define PHYSICS_INITIALIZE_QUADRANT_H

/* header files */
#include "physics_var_physics.h"

#ifdef __cplusplus
extern "C" {
#endif

/** External solver initialize quadrant with solution data.
 *
 * \ingroup solver_group
 *
 * @param [in] initial_condition    Initial condition flag
 * @param [in] dim                  Simulation spatial dimension
 * @param [in] unstructured         Unstructured grid flag
 * @param [in] nfields              Number of fields (unknowns) per grid point
 * @param [in] nsub_elem            Number of 1D subgrid elements per AMR quad
 * @param [in] cell_type            Cell type
 * @param [in] elem_h               AMR cell length (structured)
 * @param [in] elem_xyz             AMR cell geometry coordinates (corners)
 * @param [out] elem_soln           Solution
 * @param [in] mach                 Fluid Property: Mach number
 * @param [in] alpha                Fluid Property: angle of attack
 * @param [in] beta                 Fluid Property: angle of yaw
 * @param [in] gamma                Fluid Property: ratio of specific heats
 * @param [in] density              Fluid Property: freestream density
 * @param [in] pressure             Fluid Property: freestream pressure
 */
void solver_initialize_quadrant(int *initial_condition,int *dim,int *unstructured,
                             int *nfields,int *nsub_elem,int *cell_type,
                             double *elem_h,double *elem_xyz,double *elem_soln,
                             double *mach,double *alpha,double *beta,
                             double *gamma,double *density,double *pressure);

#ifdef __cplusplus
}
#endif

/** Function wrapper for external initialize quadrant with solution data.
 *
 * @param [in] dim              Simulation spatial dimension
 * @param [in] unstructured     Unstructured grid flag
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nsub_elem        Number of 1D subgrid elements per AMR quad
 * @param [in] cell_type        Cell type
 * @param [in] elem_h           AMR cell length (structured)
 * @param [in] elem_xyz         AMR cell geometry coordinates (corners)
 * @param [out] elem_soln       Solution
 */
void physics_initialize_quadrant(int *dim,int *unstructured,int *nfields,
                                 int *nsub_elem,int *cell_type,
                                 double *elem_h,double *elem_xyz,
                                 double *elem_soln);
#endif /* PHYSICS_INITIALIZE_QUADRANT_H */