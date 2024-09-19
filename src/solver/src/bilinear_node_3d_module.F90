
!     
! \file   bilinear_node_3d_module.F90
! \author akirby
!
! Created on August 30, 2018, 8:53 AM
!


module bilinear_node_3d_module
    implicit none
    contains


    subroutine bilinear_node_3d(cell_geom,node_natural,node_physical)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: cell_geom(3,8)
        real(dp),intent(in) :: node_natural(3)
        real(dp),intent(out):: node_physical(3)

        real(dp) :: oneMxi,oneMeta,oneMsigma
        real(dp) :: onePxi,onePeta,onePsigma

        oneMxi      = 0.5_dp*(1.0_dp - node_natural(1))
        onePxi      = 0.5_dp*(1.0_dp + node_natural(1))
        oneMeta     = 0.5_dp*(1.0_dp - node_natural(2))
        onePeta     = 0.5_dp*(1.0_dp + node_natural(2))
        oneMsigma   = 0.5_dp*(1.0_dp - node_natural(3))
        onePsigma   = 0.5_dp*(1.0_dp + node_natural(3))

        node_physical(:) =  cell_geom(:,1)*(oneMxi)*(oneMeta)*(oneMsigma) + &
                            cell_geom(:,2)*(onePxi)*(oneMeta)*(oneMsigma) + &
                            cell_geom(:,3)*(oneMxi)*(onePeta)*(oneMsigma) + &
                            cell_geom(:,4)*(onePxi)*(onePeta)*(oneMsigma) + &
                            cell_geom(:,5)*(oneMxi)*(oneMeta)*(onePsigma) + &
                            cell_geom(:,6)*(onePxi)*(oneMeta)*(onePsigma) + &
                            cell_geom(:,7)*(oneMxi)*(onePeta)*(onePsigma) + &
                            cell_geom(:,8)*(onePxi)*(onePeta)*(onePsigma)


    end subroutine
    
    
end module

