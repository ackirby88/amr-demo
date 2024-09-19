/**
 * \file   physics_residual_face_full_unstructured.h
 * \author akirby
 *
 * Created on August 31, 2018, 9:49 AM
 */

#ifndef PHYSICS_RESIDUAL_FACE_FULL_UNSTRUCTURED_H
#define PHYSICS_RESIDUAL_FACE_FULL_UNSTRUCTURED_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */

/** External solver function for calculating residual at 
 *        AMR quad face with no hanging elements
 * 
 *  \ingroup solver_group
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side_l           Side id (face) of AMR left quad 
 * @param [in] side_r           Side id (face) of AMR right quad 
 * @param [in] type_l           Cell type of left element
 * @param [in] type_r           Cell type of right element
 * @param [in] indl             Solution index of left element
 * @param [in] indr             Solution index of right element
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] geom_l           Left AMR cell geometry cell coordinates
 * @param [in] geom_r           Right AMR cell geometry cell coordinates
 * @param [in] ql               Solution pointer of left element
 * @param [in] qr               Solution pointer of right element
 * @param [in,out] res          Residual pointer (all elements)
 */
void solver_residual_face_full_unstructured(int *dim,int *nfields,int *nelem_subgrid,int *side_l,int *side_r,
                                            int *type_l,int *type_r,int *indl,int *indr,int *dof,
                                            double *geom_l,double *geom_r,
                                            double *ql,double *qr,double *res);
 

#ifdef __cplusplus
}
#endif

/** Solver function wrapper for calculating residual at 
 *        AMR quad face with no hanging elements
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side_l           Side id (face) of AMR left quad 
 * @param [in] side_r           Side id (face) of AMR right quad 
 * @param [in] type_l           Cell type of left element
 * @param [in] type_r           Cell type of right element
 * @param [in] indl             Solution index of left element
 * @param [in] indr             Solution index of right element
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] geom_l           Left AMR cell geometry cell coordinates
 * @param [in] geom_r           Right AMR cell geometry cell coordinates
 * @param [in] ql               Solution pointer of left element
 * @param [in] qr               Solution pointer of right element
 * @param [in,out] res          Residual pointer (all elements)
 */
void physics_residual_face_full_unstructured(int *dim,int *nfields,int *nelem_subgrid,int *side_l,int *side_r,
                                             int *type_l,int *type_r,int *indl,int *indr,int *dof,
                                             double *geom_l,double *geom_r,
                                             double *ql,double *qr,double *res);


#endif /* PHYSICS_RESIDUAL_FACE_FULL_UNSTRUCTURED_H */

