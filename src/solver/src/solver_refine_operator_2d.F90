!>     
!> \file:   solver_refine_operator_2d.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D refine operators
!>
!> Created on September 6, 2018, 11:27 AM
!>


!> 2D refine operator: Unweighted Interpolation
!>
!> @param unstructured     Flag for unstructured grid
!> @param nfields          Number of fields (unknowns) per grid point
!> @param nelem_subgrid    Number of 1D subgrid elements in each AMR quad
!> @param cell_type        Cell type
!> @param h                AMR quadrant side length (structured)
!> @param geom             AMR quadrant corner xyz geometry coordinates
!> @param soln             coarse cell solution pointer
!> @param uc1              fine cell #1 solution pointer
!> @param uc2              fine cell #2 solution pointer
!> @param uc3              fine cell #3 solution pointer
!> @param uc4              fine cell #4 solution pointer
!> 
!>      This refine operator performs a simple unweighted interpolation of the cell solutions. 
!>      For proper refinement on unstructured grids, the cell volumes should be taken into account.
!>
subroutine solver_refine_operator_2d(unstructured,nfields,nelem_subgrid,cell_type,h,geom,&
                                     soln,uc1,uc2,uc3,uc4) bind(C)
    use iso_c_binding
    implicit none
    

    integer(c_int),intent(in)   :: unstructured
    integer(c_int),intent(in)   :: nfields
    integer(c_int),intent(in)   :: nelem_subgrid
    integer(c_int),intent(in)   :: cell_type
    real(c_double),intent(in)   :: h
    real(c_double),intent(in)   :: geom(3,4)
    real(c_double),intent(in)   :: soln(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(out)  :: uc1(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(out)  :: uc2(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(out)  :: uc3(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(out)  :: uc4(nfields*nelem_subgrid*nelem_subgrid)
    


    !set children equal to parent
    uc1 = soln
    uc2 = soln
    uc3 = soln
    uc4 = soln
    

end subroutine
