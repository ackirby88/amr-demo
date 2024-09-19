
!     
! \file   cell_volume_module.F90
! \author akirby
!
! Created on August 30, 2018, 9:47 AM
!


module cell_volume_module
    implicit none
    contains


    subroutine cell_volume(xc,vol)
        use my_kinddefs
        implicit none

        real(dp),   intent(in) :: xc(3,8)
        real(dp),   intent(out):: vol


        integer(i4) :: iface
        real(dp)    :: scalarProduct

        integer(i4) :: fconn(4,6)

        
        fconn(:,1) = (/ 5,6,2,1 /)
        fconn(:,2) = (/ 5,7,8,6 /)
        fconn(:,3) = (/ 6,8,4,2 /)
        fconn(:,4) = (/ 2,4,3,1 /)
        fconn(:,5) = (/ 5,1,3,7 /)
        fconn(:,6) = (/ 7,3,4,8 /)


        vol = 0.0_dp
        do iface = 1,6

            call scalar_product(xc(:,fconn(1,iface)),xc(:,fconn(2,iface)),xc(:,fconn(3,iface)),scalarProduct)         
            vol = vol - 0.25_dp*scalarProduct
            
            call scalar_product(xc(:,fconn(1,iface)),xc(:,fconn(3,iface)),xc(:,fconn(4,iface)),scalarProduct)
            vol = vol - 0.25_dp*scalarProduct
            
            call scalar_product(xc(:,fconn(1,iface)),xc(:,fconn(2,iface)),xc(:,fconn(4,iface)),scalarProduct)
            vol = vol - 0.25_dp*scalarProduct

            call scalar_product(xc(:,fconn(2,iface)),xc(:,fconn(3,iface)),xc(:,fconn(4,iface)),scalarProduct)
            vol = vol - 0.25_dp*scalarProduct
            
        end do
        vol = vol/3.0_dp


    end subroutine


    subroutine scalar_product(a,b,c,scalarProduct)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: a(3),b(3),c(3)
        real(dp),intent(out):: scalarProduct

        scalarProduct = a(1)*b(2)*c(3) - a(1)*b(3)*c(2) &
                       +a(2)*b(3)*c(1) - a(2)*b(1)*c(3) & 
                       +a(3)*b(1)*c(2) - a(3)*b(2)*c(1)

        
    end subroutine


end module
