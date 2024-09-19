/**
 * \file   physics_rhs.h
 * \author akirby
 *
 * Created on August 24, 2018, 3:49 PM
 */

#ifndef PHYSICS_RHS_H
#define PHYSICS_RHS_H

/* header files */
#include "physics_rhs_structured.h"
#include "physics_rhs_unstructured.h"


/** Physics right-hand-side (rhs) interface operator
 * 
 * @param [in] dim              Simulation dimension
 * @param [in] unstructured     Unstructured grid flag
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
void physics_rhs(int *dim,int *unstructured,int *nfields,int *nelem_subgrid,
                 int *ncell,int *cell_info,double *geom,int *nface,
                 int *face_info,int *dof,double *soln,double *rhs);


#endif /* PHYSICS_RHS_H */

