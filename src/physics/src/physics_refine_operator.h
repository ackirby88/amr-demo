/**
 * \file   physics_refine_operator.h
 * \author akirby
 *
 * Created on September 6, 2018, 11:22 AM
 */

#ifndef PHYSICS_REFINE_OPERATOR_H
#define PHYSICS_REFINE_OPERATOR_H


#ifdef __cplusplus
extern "C" {
#endif

/* external solver subroutines */
    
/** External solver 2D refine operator
 * 
 * \ingroup solver_group
 * 
 * @param [in] unstructured     Flag for unstructured grid
 * @param [in] numfields        Number of fields (unknowns) per grid point
 * @param [in] nsub_elem        Number of 1D subgrid elements in each AMR quad
 * @param [in] cell_type        Cell type
 * @param [in] elem_h           AMR quadrant side length (structured)
 * @param [in] geom             AMR quadrant corner xyz geometry coordinates
 * @param [out] soln            coarse cell solution pointer
 * @param [in] uc1              fine cell #1 solution pointer
 * @param [in] uc2              fine cell #2 solution pointer
 * @param [in] uc3              fine cell #3 solution pointer
 * @param [in] uc4              fine cell #4 solution pointer
 */
void solver_refine_operator_2d(int *unstructured,int *numfields,int *nsub_elem,
                               int *cell_type,double *elem_h,double *geom,double *soln,
                               double *uc1,double *uc2,double *uc3,double *uc4);

/** External solver 3D refine operator
 * 
 * \ingroup solver_group
 * 
 * @param [in] unstructured     Flag for unstructured grid
 * @param [in] numfields        Number of fields (unknowns) per grid point
 * @param [in] nsub_elem        Number of 1D subgrid elements in each AMR quad
 * @param [in] cell_type        Cell type
 * @param [in] elem_h           AMR quadrant side length (structured)
 * @param [in] geom             AMR quadrant corner xyz geometry coordinates
 * @param [out] soln            coarse cell solution pointer
 * @param [in] uc1              fine cell #1 solution pointer
 * @param [in] uc2              fine cell #2 solution pointer
 * @param [in] uc3              fine cell #3 solution pointer
 * @param [in] uc4              fine cell #4 solution pointer
 * @param [in] uc5              fine cell #5 solution pointer
 * @param [in] uc6              fine cell #6 solution pointer
 * @param [in] uc7              fine cell #7 solution pointer
 * @param [in] uc8              fine cell #8 solution pointer
 */
void solver_refine_operator_3d(int *unstructured,int *numfields,int *nsub_elem,
                               int *cell_type,double *elem_h,double *geom,double *soln,
                               double *uc1,double *uc2,double *uc3,double *uc4,
                               double *uc5,double *uc6,double *uc7,double *uc8);


#ifdef __cplusplus
}
#endif

/** Physics solver 2D refine operator interface
 * 
 * @param [in] unstructured     Flag for unstructured grid
 * @param [in] numfields        Number of fields (unknowns) per grid point
 * @param [in] nsub_elem        Number of 1D subgrid elements in each AMR quad
 * @param [in] cell_type        Cell type
 * @param [in] elem_h           AMR quadrant side length (structured)
 * @param [in] geom             AMR quadrant corner xyz geometry coordinates
 * @param [out] soln            coarse cell solution pointer
 * @param [in] uc1              fine cell #1 solution pointer
 * @param [in] uc2              fine cell #2 solution pointer
 * @param [in] uc3              fine cell #3 solution pointer
 * @param [in] uc4              fine cell #4 solution pointer
 */
void physics_refine_operator_2d(int *unstructured,int *numfields,int *nsub_elem,
                                int *cell_type,double *elem_h,double *geom,double *soln,
                                double *uc1,double *uc2,double *uc3,double *uc4);

/** Physics solver 3D refine operator interface
 * 
 * @param [in] unstructured     Flag for unstructured grid
 * @param [in] numfields        Number of fields (unknowns) per grid point
 * @param [in] nsub_elem        Number of 1D subgrid elements in each AMR quad
 * @param [in] cell_type        Cell type
 * @param [in] elem_h           AMR quadrant side length (structured)
 * @param [in] geom             AMR quadrant corner xyz geometry coordinates
 * @param [out] soln            coarse cell solution pointer
 * @param [in] uc1              fine cell #1 solution pointer
 * @param [in] uc2              fine cell #2 solution pointer
 * @param [in] uc3              fine cell #3 solution pointer
 * @param [in] uc4              fine cell #4 solution pointer
 * @param [in] uc5              fine cell #5 solution pointer
 * @param [in] uc6              fine cell #6 solution pointer
 * @param [in] uc7              fine cell #7 solution pointer
 * @param [in] uc8              fine cell #8 solution pointer
 */
void physics_refine_operator_3d(int *unstructured,int *numfields,int *nsub_elem,
                                int *cell_type,double *elem_h,double *geom,double *soln,
                                double *uc1,double *uc2,double *uc3,double *uc4,
                                double *uc5,double *uc6,double *uc7,double *uc8);


#endif /* PHYSICS_REFINE_OPERATOR_H */

