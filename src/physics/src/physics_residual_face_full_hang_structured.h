/**
 * \file   physics_residual_face_full_hang_structured.h
 * \author akirby
 *
 * Created on September 5, 2018, 2:03 PM
 */

#ifndef PHYSICS_RESIDUAL_FACE_FULL_HANG_STRUCTURED_H
#define PHYSICS_RESIDUAL_FACE_FULL_HANG_STRUCTURED_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */

/** External solver function for calculating residual on 
 *        full-hang face on 2D structured grid
 * 
 * \ingroup solver_group
 * 
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face) of AMR left quad 
 * @param [in] level_l          Level of AMR left element
 * @param [in] level_r          Level of AMR right element
 * @param [in] type_l           Cell type of AMR left element
 * @param [in] type_r1          Cell type of AMR right element #1
 * @param [in] type_r2          Cell type of AMR right element #2
 * @param [in] ind_l            Solution index of AMR left element
 * @param [in] ind_r1           Solution index of AMR right element #1
 * @param [in] ind_r2           Solution index of AMR right element #2
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              Lower corner x-coordinate of AMR quad face
 * @param [in] ylo              Lower corner y-coordinate of AMR quad face
 * @param [in] zlo              Lower corner z-coordinate of AMR quad face
 * @param [in] ql               Solution pointer of left element
 * @param [in] qr_1             Solution pointer of right element #1
 * @param [in] qr_2             Solution pointer of right element #2
 * @param [in,out] rhs          Residual pointer (all elements)
 */
void solver_residual_face_full_hang_structured_2d(
                int *nfields,int *nelem_subgrid,int *side,
                int *level_l,int *level_r,
                int *type_l,int *type_r1,int *type_r2,
                int *ind_l,int *ind_r1,int *ind_r2,
                int *dof,double *h_base,double *xlo,double *ylo,double *zlo,
                double *ql,double *qr_1,double *qr_2,
                double *rhs);

/** External solver function for calculating residual on 
 *        full-hang face on 3D structured grid
 * 
 * \ingroup solver_group
 * 
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face) of AMR left quad 
 * @param [in] level_l          Level of AMR left element
 * @param [in] level_r          Level of AMR right element
 * @param [in] type_l           Cell type of AMR left element
 * @param [in] type_r1          Cell type of AMR right element #1
 * @param [in] type_r2          Cell type of AMR right element #2
 * @param [in] type_r3          Cell type of AMR right element #3
 * @param [in] type_r4          Cell type of AMR right element #4
 * @param [in] ind_l            Solution index of AMR left element
 * @param [in] ind_r1           Solution index of AMR right element #1
 * @param [in] ind_r2           Solution index of AMR right element #2
 * @param [in] ind_r3           Solution index of AMR right element #3
 * @param [in] ind_r4           Solution index of AMR right element #4
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              Lower corner x-coordinate of AMR element face
 * @param [in] ylo              Lower corner y-coordinate of AMR element face
 * @param [in] zlo              Lower corner z-coordinate of AMR element face
 * @param [in] ql               Solution pointer of AMR left element
 * @param [in] qr_1             Solution pointer of AMR right element #1
 * @param [in] qr_2             Solution pointer of AMR right element #2
 * @param [in] qr_3             Solution pointer of AMR right element #3
 * @param [in] qr_4             Solution pointer of AMR right element #4
 * @param [in,out] rhs          Residual pointer (all elements)
 */
void solver_residual_face_full_hang_structured_3d(
                int *nfields,int *nelem_subgrid,int *side,
                int *level_l,int *level_r,
                int *type_l,int *type_r1,int *type_r2,int *type_r3,int *type_r4,
                int *ind_l,int *ind_r1,int *ind_r2,int *ind_r3,int *ind_r4,
                int *dof,double *h_base,double *xlo,double *ylo,double *zlo,
                double *ql,double *qr_1,double *qr_2,double *qr_3,double *qr_4,
                double *rhs);

#ifdef __cplusplus
}
#endif


/** External solver function wrapper for calculating residual on 
 *        full-hang face on structured grid
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR quad
 * @param [in] side             Side id (face) of AMR left quad 
 * @param [in] level_l          Level of AMR left element
 * @param [in] level_r          Level of AMR right element
 * @param [in] type_l           Cell type of AMR left element
 * @param [in] type_r           Cell types of AMR right elements
 * @param [in] ind_l            Solution index of AMR left element
 * @param [in] ind_r            Solution indicies of AMR right elements
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              Lower corner x-coordinate of AMR quad face
 * @param [in] ylo              Lower corner y-coordinate of AMR quad face
 * @param [in] zlo              Lower corner z-coordinate of AMR quad face
 * @param [in] soln             Solution pointer (all elements)
 * @param [in,out] rhs          Residual pointer (all elements)
 */
void physics_residual_face_full_hang_structured(int *dim,int *nfields,int *nelem_subgrid,
                            int *side,int *level_l,int *level_r,int *type_l,int *type_r,
                            int *ind_l,int *ind_r,int *dof,double *h_base,
                            double *xlo,double *ylo,double *zlo,
                            double *soln,double *rhs);


#endif /* PHYSICS_RESIDUAL_FACE_FULL_HANG_STRUCTURED_H */
