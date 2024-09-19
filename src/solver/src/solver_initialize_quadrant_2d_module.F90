
!>
!> \file   solver_initialize_quadrant_2d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D initialize quadrant
!>
!> Created on August 30, 2018, 9:35 AM
!>

module solver_initialize_quadrant_2d_module
    implicit none
    contains

    !> 2D initialize quadrant and boundary conditions
    !>
    !> @param initial_condition     Initial condition flag
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
    subroutine solver_initialize_quadrant_2d(initial_condition,unstructured,nfields,nelem_subgrid,cell_type,h,cell_geom,Q,&
                                      mach,alpha,beta,gamma,density,pressure)
        use my_kinddefs
        use initial_boundary_2d_module
        use initial_condition_2d_module
        implicit none

        integer(i4),intent(in)      :: initial_condition
        integer(i4),intent(in)      :: unstructured
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(inout)   :: cell_type
        real(dp),   intent(in)      :: h
        real(dp),   intent(in)      :: cell_geom(3,4) !always 3 dimensional
        real(dp),   intent(out)     :: Q(nfields,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(in)      :: mach
        real(dp),   intent(in)      :: alpha
        real(dp),   intent(in)      :: beta
        real(dp),   intent(in)      :: gamma
        real(dp),   intent(in)      :: density
        real(dp),   intent(in)      :: pressure


        call initial_boundary_2d(unstructured,cell_geom,cell_type)
        call initial_condition_2d(initial_condition,unstructured,nfields,nelem_subgrid,h,cell_geom,Q,&
                                  mach,alpha,beta,gamma,density,pressure)
    end subroutine
end module