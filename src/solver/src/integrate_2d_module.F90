
!     
! \file   integrate_2d_module.F90
! \author akirby
!
! Created on August 30, 2018, 8:58 AM
!


module integrate_2d_module
    implicit none
    contains
    
    
    subroutine integrate_2d_pos(i,j,nfields,nelem_subgrid,flux,res)
        use my_kinddefs
        implicit none

        integer(i4),intent(in)      :: i
        integer(i4),intent(in)      :: j
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        real(dp),   intent(in)      :: flux(4)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid)

        res(:,i,j) = res(:,i,j) + flux(:)

    end subroutine


    subroutine integrate_2d_neg(i,j,nfields,nelem_subgrid,flux,res)
        use my_kinddefs
        implicit none

        integer(i4),intent(in)      :: i
        integer(i4),intent(in)      :: j
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        real(dp),   intent(in)      :: flux(4)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid)

        res(:,i,j) = res(:,i,j) - flux(:)

    end subroutine


end module
