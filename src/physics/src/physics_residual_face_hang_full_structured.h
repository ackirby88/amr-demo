/**
 * \file   physics_residual_face_hang_full_structured.h
 * \author akirby
 *
 * Created on September 5, 2018, 2:37 PM
 */

#ifndef PHYSICS_RESIDUAL_FACE_HANG_FULL_STRUCTURED_H
#define PHYSICS_RESIDUAL_FACE_HANG_FULL_STRUCTURED_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */
    
/** External solver function for calculating residual on 
 *        hang-full face on 2D structured grid
 * 
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR element
 * @param [in] side             Side id (face) of AMR left element 
 * @param [in] level_l          Level of AMR left element
 * @param [in] level_r          Level of AMR right element
 * @param [in] type_l1          Cell type of AMR left element #1
 * @param [in] type_l2          Cell type of AMR left element #2
 * @param [in] type_r           Cell type of AMR right element
 * @param [in] ind_l1           Solution index of AMR left element #1
 * @param [in] ind_l2           Solution index of AMR left element #2
 * @param [in] ind_r            Solution index of AMR right element
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              Lower corner x-coordinate of AMR element face
 * @param [in] ylo              Lower corner y-coordinate of AMR element face
 * @param [in] zlo              Lower corner z-coordinate of AMR element face
 * @param [in] ql_1             Solution pointer of left element #1
 * @param [in] ql_2             Solution pointer of left element #2
 * @param [in] qr               Solution pointer of right element
 * @param [in,out] rhs          Residual pointer (all elements)
 */
void solver_residual_face_hang_full_structured_2d(
                int *nfields,int *nelem_subgrid,int *side,
                int *level_l,int *level_r,
                int *type_l1,int *type_l2,int *type_r,
                int *ind_l1,int *ind_l2,int *ind_r,
                int *dof,double *h_base,double *xlo,double *ylo,double *zlo,
                double *ql_1,double *ql_2,double *qr,
                double *rhs);

/** External solver function for calculating residual on 
 *        hang-full face on 2D structured grid
 * 
 * \ingroup solver_group
 * 
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR element
 * @param [in] side             Side id (face) of AMR left element 
 * @param [in] level_l          Level of AMR left element
 * @param [in] level_r          Level of AMR right element
 * @param [in] type_l1          Cell type of AMR left element #1
 * @param [in] type_l2          Cell type of AMR left element #2
 * @param [in] type_l3          Cell type of AMR left element #3
 * @param [in] type_l4          Cell type of AMR left element #4
 * @param [in] type_r           Cell type of AMR right element
 * @param [in] ind_l1           Solution index of AMR left element #1
 * @param [in] ind_l2           Solution index of AMR left element #2
 * @param [in] ind_l3           Solution index of AMR left element #3
 * @param [in] ind_l4           Solution index of AMR left element #4
 * @param [in] ind_r            Solution index of AMR right element
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              Lower corner x-coordinate of AMR element face
 * @param [in] ylo              Lower corner y-coordinate of AMR element face
 * @param [in] zlo              Lower corner z-coordinate of AMR element face
 * @param [in] ql_1             Solution pointer of left element #1
 * @param [in] ql_2             Solution pointer of left element #2
 * @param [in] ql_3             Solution pointer of left element #3
 * @param [in] ql_4             Solution pointer of left element #4
 * @param [in] qr               Solution pointer of right element
 * @param [in,out] rhs          Residual pointer (all elements)
 */
void solver_residual_face_hang_full_structured_3d(
                int *nfields,int *nelem_subgrid,int *side,
                int *level_l,int *level_r,
                int *type_l1,int *type_l2,int *type_l3,int *type_l4,int *type_r,
                int *ind_l1,int *ind_l2,int *ind_l3,int *ind_l4,int *ind_r,
                int *dof,double *h_base,double *xlo,double *ylo,double *zlo,
                double *ql_1,double *ql_2,double *ql_3,double *ql_4,double *qr,
                double *rhs);

#ifdef __cplusplus
}
#endif


/** External solver function wrapper for calculating residual on 
 *        full-hang face on structured grid
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR element
 * @param [in] side             Side id (face) of AMR left element 
 * @param [in] level_l          Level of AMR left element
 * @param [in] level_r          Level of AMR right element
 * @param [in] type_l           Cell types of AMR left elements
 * @param [in] type_r           Cell type of AMR right elements
 * @param [in] ind_l            Solution indicies of AMR left elements
 * @param [in] ind_r            Solution index of AMR right elements
 * @param [in] dof              Number of real degrees of freedom
 * @param [in] h_base           Coarse level grid element length
 * @param [in] xlo              Lower corner x-coordinate of AMR element face
 * @param [in] ylo              Lower corner y-coordinate of AMR element face
 * @param [in] zlo              Lower corner z-coordinate of AMR element face
 * @param [in] soln             Solution pointer (all elements)
 * @param [in,out] rhs          Residual pointer (all elements)
 */
void physics_residual_face_hang_full_structured(int *dim,int *nfields,int *nelem_subgrid,
                            int *side,int *level_l,int *level_r,int *type_l,int *type_r,
                            int *ind_l,int *ind_r,int *dof,double *h_base,
                            double *xlo,double *ylo,double *zlo,
                            double *soln,double *rhs);



#endif /* PHYSICS_RESIDUAL_FACE_HANG_FULL_STRUCTURED_H */

