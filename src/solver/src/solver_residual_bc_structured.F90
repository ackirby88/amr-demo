
!>     
!> \file   solver_residual_bc_structured.F90
!> \author akirby
!> \ingroup solver_group
!> \brief Solver boundary condition residual structured grid interface
!> 
!> Created on September 4, 2018, 9:57 AM
!>


!> Solver boundary condition residual structured grid interface
!> 
!> @param dim               Simulation spatial dimension
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
!> @param side              Side id (face number) of AMR quad
!> @param cell_type         Cell type 
!> @param level             AMR level
!> @param bc                Boundary condition type
!> @param h_base            Coarsest level element width (structured)
!> @param xlo               x-coordinate of first node in AMR quad
!> @param ylo               y-coordinate of first node in AMR quad
!> @param zlo               z-coordinate of first node in AMR quad
!> @param soln              Solution pointer for this quad
!> @param res               Residual pointer for this quad
!>
subroutine solver_residual_bc_structured(dim,nfields,nelem_subgrid,side,cell_type,level,bc,h_base,xlo,ylo,zlo,soln,res) bind(C)
    use iso_c_binding
    use solver_residual_bc_structured_2d_module
    use solver_residual_bc_structured_3d_module
    implicit none

    integer(c_int),intent(in)      :: dim
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side
    integer(c_int),intent(in)      :: cell_type
    integer(c_int),intent(in)      :: level
    integer(c_int),intent(in)      :: bc
    real(c_double),intent(in)      :: h_base
    real(c_double),intent(in)      :: xlo
    real(c_double),intent(in)      :: ylo
    real(c_double),intent(in)      :: zlo
    real(c_double),intent(in)      :: soln(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(inout)   :: res(nfields,(nelem_subgrid)**dim)


    if(dim==2)then
        call solver_residual_bc_structured_2d(nfields,nelem_subgrid,side,cell_type,level,bc,h_base,xlo,ylo,zlo,soln,res)
    else
        call solver_residual_bc_structured_3d(nfields,nelem_subgrid,side,cell_type,level,bc,h_base,xlo,ylo,zlo,soln,res)
    end if
    
    
end subroutine
