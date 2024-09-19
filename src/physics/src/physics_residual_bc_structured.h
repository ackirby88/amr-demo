/**
 * \file   physics_residual_bc_structured.h
 * \author akirby
 *
 * Created on September 4, 2018, 9:54 AM
 */

#ifndef PHYSICS_RESIDUAL_BC_STRUCTURED_H
#define PHYSICS_RESIDUAL_BC_STRUCTURED_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */
    
/** External solver function for calculating boundary condition residual 
 *  on structured grids
 * 
 * \ingroup solver_group
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face number) of AMR quad
 * @param [in] type             Cell type 
 * @param [in] level            AMR level
 * @param [in] bc               Boundary condition type
 * @param [in] h_base           Coarsest level element width (structured)
 * @param [in] xlo              x-coordinate of first node in AMR quad
 * @param [in] ylo              y-coordinate of first node in AMR quad
 * @param [in] zlo              z-coordinate of first node in AMR quad
 * @param [in] soln             Solution pointer for this quad
 * @param [in,out] res          Residual pointer for this quad
 */
void solver_residual_bc_structured(int *dim,int *nfields,int *nelem_subgrid,
                                   int *side,int *type,int *level,int *bc,
                                   double *h_base,double *xlo,double *ylo,double *zlo,
                                   double *soln,double *res);


#ifdef __cplusplus
}
#endif


/** Solver function wrapper for calculating boundary condition residual 
 *  on structured grids
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face number) of AMR quad
 * @param [in] type             Cell type 
 * @param [in] level            AMR level
 * @param [in] bc               Boundary condition type
 * @param [in] h_base           Coarsest level element width (structured)
 * @param [in] xlo              x-coordinate of first node in AMR quad
 * @param [in] ylo              y-coordinate of first node in AMR quad
 * @param [in] zlo              z-coordinate of first node in AMR quad
 * @param [in] soln             Solution pointer for this quad
 * @param [in,out] res          Residual pointer for this quad
 */
void physics_residual_bc_structured(int *dim,int *nfields,int *nelem_subgrid,
                                    int *side,int *type,int *level,int *bc,
                                    double *h_base,double *xlo,double *ylo,double *zlo,
                                    double *soln,double *res);


#endif /* PHYSICS_RESIDUAL_BC_STRUCTURED_H */

