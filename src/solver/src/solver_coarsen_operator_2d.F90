
!>     
!> \file   solver_coarsen_operator_2d.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D coarsen operators
!>
!> Created on September 6, 2018, 10:37 AM
!>


!> 2D coarsen operator: Unweighted Average
!>
!> @param unstructured     Flag for unstructured grid
!> @param nfields          Number of fields (unknowns) per grid point
!> @param nelem_subgrid    Number of 1D subgrid elements in each AMR quad
!> @param cell_type        Cell type
!> @param h                AMR quadrant side length (structured)
!> @param geom             AMR quadrant corner xyz geometry coordinates
!> @param uc1              fine cell #1 solution pointer
!> @param uc2              fine cell #2 solution pointer
!> @param uc3              fine cell #3 solution pointer
!> @param uc4              fine cell #4 solution pointer
!> @param soln             coarse cell solution pointer
!> 
!>      This coarsen operator performs a simple unweighted averaging of the cell solutions. 
!>      For proper coarsening on unstructured grids, the cell volumes should be taken into account.
!>
subroutine solver_coarsen_operator_2d(unstructured,nfields,nelem_subgrid,cell_type,h,geom,&
                                      uc1,uc2,uc3,uc4,soln) bind(C)
    use iso_c_binding
    implicit none
    

    integer(c_int),intent(in)   :: unstructured
    integer(c_int),intent(in)   :: nfields
    integer(c_int),intent(in)   :: nelem_subgrid
    integer(c_int),intent(in)   :: cell_type
    real(c_double),intent(in)   :: h
    real(c_double),intent(in)   :: geom(3,4)
    real(c_double),intent(in)   :: uc1(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(in)   :: uc2(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(in)   :: uc3(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(in)   :: uc4(nfields*nelem_subgrid*nelem_subgrid)
    real(c_double),intent(out)  :: soln(nfields*nelem_subgrid*nelem_subgrid)


    !average all the values together
    soln = uc1 + uc2 + uc3 + uc4
    soln = 0.25*soln


end subroutine
