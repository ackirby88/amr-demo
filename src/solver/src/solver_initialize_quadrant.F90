
!>
!> \file   solver_initialize_quadrant.F90
!> \author akirby
!> \ingroup solver_group
!> \brief Solver initialize quadrant interface
!>
!> Created on August 28, 2018, 4:24 PM
!>


!> Solver initialize quadrant interface
!>
!> @param initial_condition     Initial condition flag
!> @param dim                   Simulation spatial dimension
!> @param unstructured          Unstructured grid flag
!> @param nfields               Number of fields (unknowns) per grid point
!> @param nelem_subgrid         Number of 1D subgrid elements per AMR quad
!> @param cell_type             Cell type
!> @param h                     AMR cell length (structured)
!> @param cell_geom             AMR cell geometry coordinates (corners)
!> @param Q                     Solution
!> @param mach                  Fluid Property: Mach number
!> @param alpha                 Fluid Property: angle of attack
!> @param beta                  Fluid Property: angle of yaw
!> @param gamma                 Fluid Property: ratio of specific heats
!> @param density               Fluid Property: freestream density
!> @param pressure              Fluid Property: freestream pressure
!>
subroutine solver_initialize_quadrant(initial_condition,dim,unstructured,nfields,nelem_subgrid,cell_type,h,cell_geom,Q,&
                                   mach,alpha,beta,gamma,density,pressure) bind(C)
    use iso_c_binding
    use solver_initialize_quadrant_2d_module
    use solver_initialize_quadrant_3d_module
    implicit none

    integer(c_int),intent(in)      :: initial_condition
    integer(c_int),intent(in)      :: dim
    integer(c_int),intent(in)      :: unstructured
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(inout)   :: cell_type
    real(c_double),intent(in)      :: h
    real(c_double),intent(in)      :: cell_geom(3,2**dim)
    real(c_double),intent(out)     :: Q(nfields,(nelem_subgrid)**dim)
    real(c_double),intent(in)      :: mach
    real(c_double),intent(in)      :: alpha
    real(c_double),intent(in)      :: beta
    real(c_double),intent(in)      :: gamma
    real(c_double),intent(in)      :: density
    real(c_double),intent(in)      :: pressure

    if(dim==2)then
        call solver_initialize_quadrant_2d(initial_condition,unstructured,nfields,nelem_subgrid,cell_type,h,cell_geom,Q,&
                                        mach,alpha,beta,gamma,density,pressure)
    else
        call solver_initialize_quadrant_3d(initial_condition,unstructured,nfields,nelem_subgrid,cell_type,h,cell_geom,Q,&
                                        mach,alpha,beta,gamma,density,pressure)
    end if
end subroutine