
!     
! \file   integrate_3d_module.F90
! \author akirby
!
! Created on August 30, 2018, 8:59 AM
!


module integrate_3d_module
    implicit none
    contains

    
    subroutine integrate_3d_pos(i,j,k,nfields,nelem_subgrid,flux,res)
        use my_kinddefs
        implicit none

        integer(i4),intent(in)      :: i
        integer(i4),intent(in)      :: j
        integer(i4),intent(in)      :: k
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        real(dp),   intent(in)      :: flux(5)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)

        res(:,i,j,k) = res(:,i,j,k) + flux(:)

    end subroutine


    subroutine integrate_3d_neg(i,j,k,nfields,nelem_subgrid,flux,res)
        use my_kinddefs
        implicit none

        integer(i4),intent(in)      :: i
        integer(i4),intent(in)      :: j
        integer(i4),intent(in)      :: k
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        real(dp),   intent(in)      :: flux(5)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)

        res(:,i,j,k) = res(:,i,j,k) - flux(:)

    end subroutine


end module
