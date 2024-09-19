/* 
 * \file   physics_time_step.h
 * \author akirby
 *
 * Created on September 19, 2018, 11:21 AM
 */

#ifndef PHYSICS_TIME_STEP_H
#define PHYSICS_TIME_STEP_H

/* header files */
#include "physics_var_physics.h"
#include "physics_var_interface.h"


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */
    
/** External solver function for calculating time step in each AMR element
 * 
 * \ingroup solver_group
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) at each grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in each AMR quad
 * @param [in] type             Cell type
 * @param [in] level            AMR level
 * @param [in] cell_geom        AMR cell geometry coordinates
 * @param [in] soln             Solution pointer for this AMR quad
 * @param [in] cfl              CFL user input
 * @param [in,out] dt           Max allowable time step
 */
void solver_time_step(int *dim,int *nfields,int *nelem_subgrid,int *type,int *level,double *cell_geom,double *soln,double *cfl,double *dt);


#ifdef __cplusplus
}
#endif


/** Function wrapper for external time step calculation.
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements per AMR quad
 * @param [in] ncell            Number of real cells
 * @param [in] cell_info        Cell info data structure
 * @param [in] geom             Cell geometries
 * @param [in] soln             Solution storage pointer
 * @param [in,out] dt           Maximum stable time step
 */
void physics_time_step(int *dim,int *nfields,int *nelem_subgrid,
                       int *ncell,int *cell_info,double *geom,double *soln, 
                       double *dt);


#endif /* PHYSICS_TIME_STEP_H */

