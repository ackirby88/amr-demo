!>     
!> \file   solver_time_step.F90
!> \author akirby
!> \ingroup solver_group
!> \brief Solver time step function interface
!>
!> Created on September 19, 2018, 11:41 AM
!>

!> Solver time step function interface
!> 
!> @param dim               Simulation spatial dimension
!> @param nfields           Number of fields (unknowns) at each grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in each AMR quad
!> @param cell_type         Cell type
!> @param level             AMR level
!> @param cell_geom         AMR cell geometry coordinates
!> @param Q                 Solution pointer for this AMR quad
!> @param cfl               CFL user input
!> @param dt                Maximum stable time step
!> 
subroutine solver_time_step(dim,nfields,nelem_subgrid,cell_type,level,cell_geom,Q,cfl,dt) bind(C)
    use iso_c_binding
    use solver_time_step_2d_module
    use solver_time_step_3d_module
    implicit none
    
    integer(c_int),intent(in)      :: dim
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: cell_type
    integer(c_int),intent(in)      :: level
    real(c_double),intent(in)      :: cell_geom(3,2**dim)
    real(c_double),intent(in)      :: Q(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(in)      :: cfl
    real(c_double),intent(inout)   :: dt
    
    
    if(dim==2)then
        call solver_time_step_2d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,cfl,dt)
    else
        call solver_time_step_3d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,cfl,dt)
    end if
    
    
end subroutine
