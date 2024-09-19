
!>     
!> \file   solver_residual_face_full_structured.F90
!> \author akirby
!> \ingroup solver_group
!> \brief Solver full face residual structured grid interface
!>
!> Created on August 28, 2018, 12:58 PM
!>


!> Solver full face residual structured grid interface
!>
!> @param dim               Simulation spatial dimension
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
!> @param side              Side id (face) of AMR left quad 
!> @param level             AMR level
!> @param cell_type_l       Cell type of left element
!> @param cell_type_r       Cell type of right element
!> @param ind_l             Solution index of left element
!> @param ind_r             Solution index of right element
!> @param dof               Number of real degrees of freedom
!> @param h_base            Coarse level grid element length
!> @param xlo               lower corner x-coordinate of AMR quad
!> @param ylo               lower corner y-coordinate of AMR quad
!> @param zlo               lower corner z-coordinate of AMR quad
!> @param ql                Solution pointer of left element
!> @param qr                Solution pointer of right element
!> @param res               Residual pointer (all elements)
!>
subroutine solver_residual_face_full_structured(dim,nfields,nelem_subgrid,side,level,&
                                  cell_type_l,cell_type_r,ind_l,ind_r,dof,h_base,xlo,ylo,zlo,ql,qr,res) bind(C)
    use iso_c_binding
    use solver_residual_face_full_structured_2d_module
    use solver_residual_face_full_structured_3d_module
    implicit none
    
    integer(c_int),intent(in)      :: dim
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side
    integer(c_int),intent(in)      :: level
    integer(c_int),intent(in)      :: cell_type_l
    integer(c_int),intent(in)      :: cell_type_r
    integer(c_int),intent(in)      :: ind_l
    integer(c_int),intent(in)      :: ind_r
    integer(c_int),intent(in)      :: dof
    real(c_double),intent(in)      :: h_base
    real(c_double),intent(in)      :: xlo(2)
    real(c_double),intent(in)      :: ylo(2)
    real(c_double),intent(in)      :: zlo(2)
    real(c_double),intent(in)      :: ql(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(in)      :: qr(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(inout)   :: res(dof)
    
    
    
    if(dim==2)then
        call solver_residual_face_full_structured_2d(nfields,nelem_subgrid,side,level,&
                                  cell_type_l,cell_type_r,ind_l,ind_r,dof,h_base,xlo,ylo,zlo,ql,qr,res)
    else
        call solver_residual_face_full_structured_3d(nfields,nelem_subgrid,side,level,&
                                  cell_type_l,cell_type_r,ind_l,ind_r,dof,h_base,xlo,ylo,zlo,ql,qr,res)
    end if
        
    
end subroutine
