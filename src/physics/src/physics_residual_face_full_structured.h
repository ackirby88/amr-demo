/**
 * \file   physics_residual_face_full_structured.h
 * \author akirby
 *
 * Created on August 28, 2018, 12:44 PM
 */

#ifndef PHYSICS_RESIDUAL_FACE_FULL_STRUCTURED_H
#define PHYSICS_RESIDUAL_FACE_FULL_STRUCTURED_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */
    
/** External solver function for calculating residual at 
 *        AMR quad face with no hanging elements
 * 
 * \ingroup solver_group
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face) of AMR left quad 
 * @param [in] level            AMR level
 * @param [in] type_l           Cell type of left element
 * @param [in] type_r           Cell type of right element
 * @param [in] indl             Solution index of left element
 * @param [in] indr             Solution index of right element
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              lower corner x-coordinate of AMR quad
 * @param [in] ylo              lower corner y-coordinate of AMR quad
 * @param [in] zlo              lower corner z-coordinate of AMR quad
 * @param [in] ql               Solution pointer of left element
 * @param [in] qr               Solution pointer of right element
 * @param [in,out] res          Residual pointer (all elements)
 */
void solver_residual_face_full_structured(int *dim,int *nfields,int *nelem_subgrid,int *side,int *level,
                                          int *type_l,int *type_r,int *indl,int *indr,int *dof,
                                          double *h_base,double *xlo,double *ylo,double *zlo,
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
 * @param [in] side             Side id (face) of AMR left quad 
 * @param [in] level            AMR level
 * @param [in] type_l           Cell type of left element
 * @param [in] type_r           Cell type of right element
 * @param [in] indl             Solution index of left element
 * @param [in] indr             Solution index of right element
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              lower corner x-coordinate of AMR quad
 * @param [in] ylo              lower corner y-coordinate of AMR quad
 * @param [in] zlo              lower corner z-coordinate of AMR quad
 * @param [in] ql               Solution pointer of left element
 * @param [in] qr               Solution pointer of right element
 * @param [in,out] res          Residual pointer (all elements)
 */
void physics_residual_face_full_structured(int *dim,int *nfields,int *nelem_subgrid,int *side,int *level,
                                           int *type_l,int *type_r,int *indl,int *indr,int *dof,
                                           double *h_base,double *xlo,double *ylo,double *zlo,
                                           double *ql,double *qr,double *res);


#endif /* PHYSICS_RESIDUAL_FACE_FULL_STRUCTURED_H */

