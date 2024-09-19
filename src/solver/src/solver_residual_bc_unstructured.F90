
!>     
!> \file   solver_residual_bc_unstructured.F90
!> \author akirby
!> \ingroup solver_group
!> \brief Solver boundary condition residual unstructured grid interface
!>
!> Created on August 31, 2018, 12:02 PM
!>


!> Solver boundary condition residual unstructured grid interface
!>
!> @param dim               Simulation spatial dimension
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
!> @param side              Side id (face number) of AMR quad
!> @param cell_type         Cell type 
!> @param bc                Boundary condition type
!> @param geom              Cell geometry coordinates
!> @param Q                 Solution pointer for this quad
!> @param res               Residual pointer for this quad
!>
subroutine solver_residual_bc_unstructured(dim,nfields,nelem_subgrid,side,cell_type,bc,geom,Q,res) bind(C)
    use iso_c_binding
    use solver_residual_bc_unstructured_2d_module
    use solver_residual_bc_unstructured_3d_module
    implicit none
    
    integer(c_int),intent(in)      :: dim
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side
    integer(c_int),intent(in)      :: cell_type
    integer(c_int),intent(in)      :: bc
    real(c_double),intent(in)      :: geom(3,2**dim)
    real(c_double),intent(in)      :: Q(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(inout)   :: res(nfields,(nelem_subgrid)**dim)
    
    
    if(dim==2)then
        call solver_residual_bc_unstructured_2d(nfields,nelem_subgrid,side,cell_type,bc,geom,Q,res)
    else
        call solver_residual_bc_unstructured_3d(nfields,nelem_subgrid,side,cell_type,bc,geom,Q,res)
    end if
        
    
end subroutine
