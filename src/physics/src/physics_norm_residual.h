/**
 * \file   physics_norm_residual.h
 * \author akirby
 *
 * Created on August 29, 2018, 9:26 AM
 */

#ifndef PHYSICS_NORM_RESIDUAL_H
#define PHYSICS_NORM_RESIDUAL_H

/* header files */
#include "physics_var_interface.h"

/** L2 discrete vector norm
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements per AMR quad
 * @param [in] ncell            Number of real cells
 * @param [in] cell_info        Cell info data structure
 * @param [in] residual         Residual vector used to calculate discrete norm
 * @return L2 discrete vector norm 
 */
double physics_norm_L2_residual(int dim,int nfields,int nelem_subgrid,int ncell,int* cell_info,double *residual);


#endif /* PHYSICS_NORM_RESIDUAL_H */

