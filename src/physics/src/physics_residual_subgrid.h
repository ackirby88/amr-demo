/**
 * \file   physics_residual_subgrid.h
 * \author akirby
 *
 * Created on August 27, 2018, 11:10 AM
 */

#ifndef PHYSICS_RESIDUAL_SUBGRID_H
#define PHYSICS_RESIDUAL_SUBGRID_H


#ifdef __cplusplus  
extern "C" {
#endif

/* external solver subroutines */
    
/** External solver function for calculating residual of subgrid in each AMR quad
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
 * @param [in,out] res          Residual pointer for this AMR quad
 */
void solver_residual_subgrid(int *dim,int *nfields,int *nelem_subgrid,int *type,int *level,double *cell_geom,double *soln,double *res);


#ifdef __cplusplus
}
#endif


/** External solver function wrapper for calculating residual of subgrid in each AMR quad
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) at each grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in each AMR quad
 * @param [in] type             Cell type
 * @param [in] level            AMR level
 * @param [in] cell_geom        AMR cell geometry coordinates
 * @param [in] soln             Solution pointer for this AMR quad
 * @param [in,out] res          Residual pointer for this AMR quad
 */
void physics_residual_subgrid(int *dim,int* nfields,int *nelem_subgrid,int *type,int *level,double *cell_geom,double *soln,double *res);


#endif /* PHYSICS_RESIDUAL_SUBGRID_H */

