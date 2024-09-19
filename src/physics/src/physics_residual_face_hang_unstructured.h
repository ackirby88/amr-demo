/**
 * \file   physics_residual_face_hang_unstructured.h
 * \author akirby
 *
 * Created on September 5, 2018, 4:01 PM
 */

#ifndef PHYSICS_RESIDUAL_FACE_HANG_UNSTRUCTURED_H
#define PHYSICS_RESIDUAL_FACE_HANG_UNSTRUCTURED_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */

/** External solver function for calculating residual on hanging face on 2D unstructured
 * 
 * \ingroup solver_group
 * 
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR element
 * @param [in] side_l           Side of AMR left element
 * @param [in] side_r1          Side of AMR right element #1
 * @param [in] side_r2          Side of AMR right element #2
 * @param [in] type_l           Cell type of AMR left element
 * @param [in] type_r1          Cell type of AMR right element #1
 * @param [in] type_r2          Cell type of AMR right element #2
 * @param [in] ind_l            Solution index of AMR left element
 * @param [in] ind_r1           Solution index of AMR right element #1
 * @param [in] ind_r2           Solution index of AMR right element #2
 * @param [in] dof              Total number of real degrees of freedom
 * @param [in] geom_l           Geometry coordinates of AMR left element
 * @param [in] geom_r1          Geometry coordinates of AMR right element #1
 * @param [in] geom_r2          Geometry coordinates of AMR right element #2
 * @param [in] ql               Solution pointer of AMR left element
 * @param [in] qr_1             Solution pointer of AMR right element #1
 * @param [in] qr_2             Solution pointer of AMR right element #2
 * @param [in,out] rhs          Residual pointer
 */
void solver_residual_face_hang_unstructured_2d(
            int *nfields,int *nelem_subgrid,
            int *side_l,int *side_r1,int *side_r2,
            int *type_l,int *type_r1,int *type_r2,
            int *ind_l,int *ind_r1,int *ind_r2,
            int *dof,
            double *geom_l,
            double *geom_r1,
            double *geom_r2,
            double *ql,
            double *qr_1,
            double *qr_2,
            double *rhs);

/** External solver function for calculating residual on hanging face on 3D unstructured
 * 
 * \ingroup solver_group
 * 
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR element
 * @param [in] side_l           Side of AMR left element
 * @param [in] side_r1          Side of AMR right element #1
 * @param [in] side_r2          Side of AMR right element #2
 * @param [in] side_r3          Side of AMR right element #3
 * @param [in] side_r4          Side of AMR right element #4
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
 * @param [in] dof              Total number of real degrees of freedom
 * @param [in] geom_l           Geometry coordinates of AMR left element
 * @param [in] geom_r1          Geometry coordinates of AMR right element #1
 * @param [in] geom_r2          Geometry coordinates of AMR right element #2
 * @param [in] geom_r3          Geometry coordinates of AMR right element #3
 * @param [in] geom_r4          Geometry coordinates of AMR right element #4
 * @param [in] ql               Solution pointer of AMR left element
 * @param [in] qr_1             Solution pointer of AMR right element #1
 * @param [in] qr_2             Solution pointer of AMR right element #2
 * @param [in] qr_3             Solution pointer of AMR right element #3
 * @param [in] qr_4             Solution pointer of AMR right element #4
 * @param [in,out] rhs          Residual pointer
 */
void solver_residual_face_hang_unstructured_3d(
            int *nfields,int *nelem_subgrid,
            int *side_l,int *side_r1,int *side_r2,int *side_r3,int *side_r4,
            int *type_l,int *type_r1,int *type_r2,int *type_r3,int *type_r4,
            int *ind_l,int *ind_r1,int *ind_r2,int *ind_r3,int *ind_r4,
            int *dof,
            double *geom_l,
            double *geom_r1,
            double *geom_r2,
            double *geom_r3,
            double *geom_r4,
            double *ql,
            double *qr_1,
            double *qr_2,
            double *qr_3,
            double *qr_4,
            double *rhs);



#ifdef __cplusplus
}
#endif

/** External solver function wrapper for calculating residual on hanging face on unstructured
 * 
 * @param [in] dim              Simulation spatial dimension
 * @param [in] nfields          Number of fields (unknowns) per grid point
 * @param [in] nelem_subgrid    Number of 1D subgrid elements in AMR element
 * @param [in] side_l           Side of AMR left element
 * @param [in] side_r           Side of AMR right element
 * @param [in] type_l           Cell type of AMR left element
 * @param [in] type_r           Cell type of AMR right element
 * @param [in] cell_l           Cell index into geom of AMR left element
 * @param [in] cell_r           Cell index into geom of AMR right element
 * @param [in] ind_l            Cell index into soln of AMR left element
 * @param [in] ind_r            Cell index into soln of AMR right element
 * @param [in] dof              Total number of real degrees of freedom
 * @param [in] geom             Cell geometry coordinates pointer
 * @param [in] soln             Solution pointer
 * @param [in,out] rhs          Residual pointer
 */
void physics_residual_face_hang_unstructured(int *dim,int *nfields,int *nelem_subgrid,
                                             int *side_l,int *side_r,
                                             int *type_l,int *type_r,
                                             int *cell_l,int *cell_r,
                                             int *ind_l,int *ind_r,
                                             int *dof,double *geom,
                                             double *soln,double *rhs);


#endif /* PHYSICS_RESIDUAL_FACE_HANG_UNSTRUCTURED_H */

