!>     
!> \file   solver_residual_face_full_unstructured.F90
!> \author akirby
!> \ingroup solver_group
!> \brief Solver full face residual unstructured grid interface
!>
!> Created on August 31, 2018, 9:53 AM
!>


!> Solver full face residual unstructured grid interface
!>
!> @param dim              Simulation spatial dimension
!> @param nfields          Number of fields (unknowns) per grid point
!> @param nelem_subgrid    Number of 1D subgrid elements in AMR quad
!> @param side_l           Side id (face) of AMR left quad 
!> @param side_r           Side id (face) of AMR right quad 
!> @param cell_type_l      Cell type of left element
!> @param cell_type_r      Cell type of right element
!> @param ind_l            Solution index of left element
!> @param ind_r            Solution index of right element
!> @param dof              Number of real degrees of freedom
!> @param geom_l           Left AMR cell geometry cell coordinates
!> @param geom_r           Right AMR cell geometry cell coordinates
!> @param ql               Solution pointer of left element
!> @param qr               Solution pointer of right element
!> @param res              Residual pointer (all elements)
!>
subroutine solver_residual_face_full_unstructured(dim,nfields,nelem_subgrid,side_l,side_r,&
                                  cell_type_l,cell_type_r,ind_l,ind_r,dof,geom_l,geom_r,ql,qr,res) bind(C)
    use iso_c_binding
    use solver_residual_face_full_unstructured_2d_module
    use solver_residual_face_full_unstructured_3d_module
    implicit none

    integer(c_int),intent(in)      :: dim
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side_l
    integer(c_int),intent(in)      :: side_r
    integer(c_int),intent(in)      :: cell_type_l
    integer(c_int),intent(in)      :: cell_type_r
    integer(c_int),intent(in)      :: ind_l
    integer(c_int),intent(in)      :: ind_r
    integer(c_int),intent(in)      :: dof
    real(c_double),intent(in)      :: geom_l(3,2**dim)
    real(c_double),intent(in)      :: geom_r(3,2**dim)
    real(c_double),intent(in)      :: ql(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(in)      :: qr(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(inout)   :: res(dof)



    if(dim==2)then
        call solver_residual_face_full_unstructured_2d(nfields,nelem_subgrid,side_l,side_r,&
                                  cell_type_l,cell_type_r,ind_l,ind_r,dof,geom_l,geom_r,ql,qr,res)
    else
        call solver_residual_face_full_unstructured_3d(nfields,nelem_subgrid,side_l,side_r,&
                                  cell_type_l,cell_type_r,ind_l,ind_r,dof,geom_l,geom_r,ql,qr,res)
    end if


end subroutine
