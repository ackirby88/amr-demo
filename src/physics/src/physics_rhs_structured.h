/**
 * \file   physics_rhs_structured.h
 * \author akirby
 *
 * Created on August 29, 2018, 4:45 PM
 */

#ifndef PHYSICS_RHS_STRUCTURED_H
#define PHYSICS_RHS_STRUCTURED_H

/* header files */
#ifndef DOXYGEN_IGNORE
#include <math.h>
#endif

#include "physics_var_physics.h"
#include "physics_var_interface.h"
#include "physics_residual_subgrid.h"
#include "physics_residual_bc_structured.h"
#include "physics_residual_face_full_structured.h"
#include "physics_residual_face_full_hang_structured.h"
#include "physics_residual_face_hang_full_structured.h"


/** Structured right-hand-side operator function
 * 
 * @param [in] dim              Simulation dimension
 * @param [in] nfields          Number of unknows at a solution point
 * @param [in] nelem_subgrid    Number of 1D sub-elements in each amr quad
 * @param [in] ncell            Number of real cells
 * @param [in] cell_info        Cell info data structure
 * @param [in] geom             Cell geometries
 * @param [in] nface            Number of quad faces
 * @param [in] face_info        Face info data structure
 * @param [in] dof              Total number of real degrees of freedom
 * @param [in] soln             Solution storage pointer
 * @param [out] rhs             Right-hand-side storage pointer
 */
void physics_rhs_structured(int *dim,int *nfields,int *nelem_subgrid,int *ncell,
                            int *cell_info,double *geom,
                            int *nface,int *face_info,int *dof,
                            double *soln,double *rhs);


#endif /* PHYSICS_RHS_STRUCTURED_H */

