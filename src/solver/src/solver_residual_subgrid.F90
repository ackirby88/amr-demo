!>     
!> \file   solver_residual_subgrid.F90
!> \author akirby
!> \ingroup solver_group
!> \brief Solver subgrid residual interface
!>
!> Created on August 27, 2018, 11:47 AM
!>


!> Solver subgrid residual interface
!> 
!> @param dim               Simulation spatial dimension
!> @param nfields           Number of fields (unknowns) at each grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in each AMR quad
!> @param cell_type         Cell type
!> @param level             AMR level
!> @param cell_geom         AMR cell geometry coordinates
!> @param Q                 Solution pointer for this AMR quad
!> @param res               Residual pointer for this AMR quad
!> 
subroutine solver_residual_subgrid(dim,nfields,nelem_subgrid,cell_type,level,cell_geom,Q,res) bind(C)
    use iso_c_binding
    use solver_residual_subgrid_2d_module
    use solver_residual_subgrid_3d_module
    implicit none
    
    integer(c_int),intent(in)      :: dim
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: cell_type
    integer(c_int),intent(in)      :: level
    real(c_double),intent(in)      :: cell_geom(3,2**dim)
    real(c_double),intent(in)      :: Q(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(inout)   :: res(nfields,(nelem_subgrid)**dim)
    
    
    if(dim==2)then
        call solver_residual_subgrid_2d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,res)
    else
        call solver_residual_subgrid_3d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,res)
    end if
    
    
end subroutine
