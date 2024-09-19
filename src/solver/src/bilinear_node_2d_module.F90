
!     
! \file   bilinear_node_2d_module.F90
! \author akirby
!
! Created on August 28, 2018, 3:08 PM
!

module bilinear_node_2d_module
    implicit none
    contains
    
    
    subroutine bilinear_node_2d(cell_geom,node_natural,node_physical)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: cell_geom(3,4)
        real(dp),intent(in) :: node_natural(2)
        real(dp),intent(out):: node_physical(2)

        real(dp) :: oneMxi,oneMeta
        real(dp) :: onePxi,onePeta

        oneMxi      = 0.5_dp*(1.0_dp - node_natural(1))
        onePxi      = 0.5_dp*(1.0_dp + node_natural(1))
        oneMeta     = 0.5_dp*(1.0_dp - node_natural(2))
        onePeta     = 0.5_dp*(1.0_dp + node_natural(2))

        node_physical(:) =  cell_geom(1:2,1)*(oneMxi)*(oneMeta) + &
                            cell_geom(1:2,2)*(onePxi)*(oneMeta) + &
                            cell_geom(1:2,3)*(oneMxi)*(onePeta) + &
                            cell_geom(1:2,4)*(onePxi)*(onePeta)


    end subroutine
    
    
end module
