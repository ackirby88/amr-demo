/**
 * \file   physics_residual_bc_unstructured.h
 * \author akirby
 * 
 * Created on August 31, 2018, 11:57 AM
 */

#ifndef PHYSICS_RESIDUAL_BC_UNSTRUCTURED_H
#define PHYSICS_RESIDUAL_BC_UNSTRUCTURED_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */
    
/** External solver function for calculating boundary condition residual 
 *  on unstructured grids
 * 
 * \ingroup solver_group
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face number) of AMR quad
 * @param [in] type             Cell type 
 * @param [in] bc               Boundary condition type
 * @param [in] geom             Cell geometry coordinates
 * @param [in] soln             Solution pointer for this quad
 * @param [in,out] res          Residual pointer for this quad
 */
void solver_residual_bc_unstructured(int *dim,int *nfields,int *nelem_subgrid,int *side,
                                    int *type,int *bc,double *geom,double *soln,double *res);


#ifdef __cplusplus
}
#endif


/** Solver function wrapper for calculating boundary condition residual 
 *  on unstructured grids
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face number) of AMR quad
 * @param [in] type             Cell type 
 * @param [in] bc               Boundary condition type
 * @param [in] geom             Cell geometry coordinates
 * @param [in] soln             Solution pointer for this quad
 * @param [in,out] res          Residual pointer for this quad
 */
void physics_residual_bc_unstructured(int *dim,int *nfields,int *nelem_subgrid,int *side,
                                      int *type,int *bc,double *geom,double *soln,double *res);


#endif /* PHYSICS_RESIDUAL_BC_UNSTRUCTURED_H */
