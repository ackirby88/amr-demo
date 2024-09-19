
!     
! \file   normals_3d_module.F90
! \author akirby
!
! Created on August 30, 2018, 8:48 AM
!


module normals_3d_module
    implicit none
    contains
    
    
    subroutine normals_3d(nelem_subgrid,side,cell_geom,normvec)
        use my_kinddefs
        use bilinear_node_3d_module
        implicit none
        
        integer(i4),intent(in) :: nelem_subgrid
        integer(i4),intent(in) :: side
        real(dp),   intent(in) :: cell_geom(3,8)
        real(dp),   intent(out):: normvec(3,nelem_subgrid,nelem_subgrid)

        
        integer(i4) :: i,j,k
        
        real(dp)    :: vec1(3)
        real(dp)    :: vec2(3)

        real(dp)    :: node1(3)
        real(dp)    :: node2(3)
        real(dp)    :: node3(3)
        real(dp)    :: node4(3)
        real(dp)    :: node1_natural(3)
        real(dp)    :: node2_natural(3)
        real(dp)    :: node3_natural(3)
        real(dp)    :: node4_natural(3)
        real(dp)    :: scale_factor
        real(dp)    :: val1,val2
        
        
        scale_factor = 2.0_dp / real(nelem_subgrid,dp)
        
        
        ! bilinear mapping
        !
        !
        !     n3 o------------o n4
        !       /|           /|     n1 = (-1,-1,-1),   N1 = 0.125_dp*(1-xi)*(1-eta)*(1-sigma)
        !      / |          / |     n2 = ( 1,-1,-1),   N2 = 0.125_dp*(1+xi)*(1-eta)*(1-sigma)
        !     /  |         /  |     n3 = (-1, 1,-1),   N3 = 0.125_dp*(1-xi)*(1+eta)*(1-sigma)
        ! n7 o------------o n8|     n4 = ( 1, 1,-1),   N4 = 0.125_dp*(1+xi)*(1+eta)*(1-sigma)
        !    |n1 o--------|---o n2  n5 = (-1,-1, 1),   N5 = 0.125_dp*(1-xi)*(1-eta)*(1+sigma)
        !    |  /         |  /      n6 = ( 1,-1, 1),   N6 = 0.125_dp*(1+xi)*(1-eta)*(1+sigma)
        !    | /          | /       n7 = (-1, 1, 1),   N7 = 0.125_dp*(1-xi)*(1+eta)*(1+sigma)
        !    |/           |/        n8 = ( 1, 1, 1),   N4 = 0.125_dp*(1+xi)*(1+eta)*(1+sigma)
        ! n5 o------------o n6
        
        
        
        select case(side)
            
            case(0) !xlo
                
                node1_natural(1) = -1.0_dp
                node2_natural(1) = -1.0_dp
                node3_natural(1) = -1.0_dp
                node4_natural(1) = -1.0_dp
                    
                do k = 1,nelem_subgrid
                    
                    val1 = scale_factor*real(k-1,dp) - 1.0_dp
                    val2 = scale_factor*real(k,dp) - 1.0_dp
                    
                    node1_natural(3) = val1
                    node2_natural(3) = val1
                    node3_natural(3) = val2
                    node4_natural(3) = val2
                    
                    do j = 1,nelem_subgrid
                        
                        val1 = scale_factor*real(j-1,dp) - 1.0_dp
                        val2 = scale_factor*real(j,dp) - 1.0_dp
                        
                        node1_natural(2) = val1
                        node2_natural(2) = val2
                        node3_natural(2) = val1
                        node4_natural(2) = val2
                        
                        call bilinear_node_3d(cell_geom,node1_natural,node1)
                        call bilinear_node_3d(cell_geom,node2_natural,node2)
                        call bilinear_node_3d(cell_geom,node3_natural,node3)
                        call bilinear_node_3d(cell_geom,node4_natural,node4)
                        
                        ! calculate normal vector
                        vec1 = node4 - node1
                        vec2 = node3 - node2

                        ! normal vector
                        normvec(1,j,k) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                        normvec(2,j,k) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                        normvec(3,j,k) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))

                    end do
                    
                end do
                
            case(1) !xhi
                
                node1_natural(1) = 1.0_dp
                node2_natural(1) = 1.0_dp
                node3_natural(1) = 1.0_dp
                node4_natural(1) = 1.0_dp
                    
                do k = 1,nelem_subgrid
                    
                    val1 = scale_factor*real(k-1,dp) - 1.0_dp
                    val2 = scale_factor*real(k,dp) - 1.0_dp
                    
                    node1_natural(3) = val1
                    node2_natural(3) = val1
                    node3_natural(3) = val2
                    node4_natural(3) = val2
                    
                    do j = 1,nelem_subgrid
                        
                        val1 = scale_factor*real(j-1,dp) - 1.0_dp
                        val2 = scale_factor*real(j,dp) - 1.0_dp
                        
                        node1_natural(2) = val1
                        node2_natural(2) = val2
                        node3_natural(2) = val1
                        node4_natural(2) = val2
                        
                        call bilinear_node_3d(cell_geom,node1_natural,node1)
                        call bilinear_node_3d(cell_geom,node2_natural,node2)
                        call bilinear_node_3d(cell_geom,node3_natural,node3)
                        call bilinear_node_3d(cell_geom,node4_natural,node4)
                        
                        ! calculate normal vector
                        vec1 = node4 - node1
                        vec2 = node3 - node2

                        ! normal vector
                        normvec(1,j,k) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                        normvec(2,j,k) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                        normvec(3,j,k) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))

                    end do
                    
                end do
                
            case(2) !ylo
                
                node1_natural(2) = -1.0_dp
                node2_natural(2) = -1.0_dp
                node3_natural(2) = -1.0_dp
                node4_natural(2) = -1.0_dp
                
                do k = 1,nelem_subgrid
                    
                    val1 = scale_factor*real(k-1,dp) - 1.0_dp
                    val2 = scale_factor*real(k,dp) - 1.0_dp
                    
                    node1_natural(3) = val1
                    node2_natural(3) = val1
                    node3_natural(3) = val2
                    node4_natural(3) = val2
                    
                    do i = 1,nelem_subgrid
                        
                        val1 = scale_factor*real(i-1,dp) - 1.0_dp
                        val2 = scale_factor*real(i,dp) - 1.0_dp

                        node1_natural(1) = val1
                        node2_natural(1) = val2
                        node3_natural(1) = val1
                        node4_natural(1) = val2
                        
                        call bilinear_node_3d(cell_geom,node1_natural,node1)
                        call bilinear_node_3d(cell_geom,node2_natural,node2)
                        call bilinear_node_3d(cell_geom,node3_natural,node3)
                        call bilinear_node_3d(cell_geom,node4_natural,node4)
                        
                        ! calculate normal vector
                        vec1 = node4 - node1
                        vec2 = node3 - node2

                        ! normal vector
                        normvec(1,i,k) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                        normvec(2,i,k) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                        normvec(3,i,k) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))
                    
                    end do
                    
                end do
                
            case(3) !yhi
                
                node1_natural(2) = 1.0_dp
                node2_natural(2) = 1.0_dp
                node3_natural(2) = 1.0_dp
                node4_natural(2) = 1.0_dp
                
                do k = 1,nelem_subgrid
                    
                    val1 = scale_factor*real(k-1,dp) - 1.0_dp
                    val2 = scale_factor*real(k,dp) - 1.0_dp
                    
                    node1_natural(3) = val1
                    node2_natural(3) = val1
                    node3_natural(3) = val2
                    node4_natural(3) = val2
                    
                    do i = 1,nelem_subgrid
                        
                        val1 = scale_factor*real(i-1,dp) - 1.0_dp
                        val2 = scale_factor*real(i,dp) - 1.0_dp

                        node1_natural(1) = val1
                        node2_natural(1) = val2
                        node3_natural(1) = val1
                        node4_natural(1) = val2
                        
                        call bilinear_node_3d(cell_geom,node1_natural,node1)
                        call bilinear_node_3d(cell_geom,node2_natural,node2)
                        call bilinear_node_3d(cell_geom,node3_natural,node3)
                        call bilinear_node_3d(cell_geom,node4_natural,node4)
                        
                        ! calculate normal vector
                        vec1 = node4 - node1
                        vec2 = node3 - node2

                        ! normal vector
                        normvec(1,i,k) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                        normvec(2,i,k) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                        normvec(3,i,k) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))
                    
                    end do
                    
                end do
                
            case(4) !zlo
                
                node1_natural(3) = -1.0_dp
                node2_natural(3) = -1.0_dp
                node3_natural(3) = -1.0_dp
                node4_natural(3) = -1.0_dp
                
                do j = 1,nelem_subgrid
                    
                    val1 = scale_factor*real(j-1,dp) - 1.0_dp
                    val2 = scale_factor*real(j,dp) - 1.0_dp
                    
                    node1_natural(2) = val1
                    node2_natural(2) = val1
                    node3_natural(2) = val2
                    node4_natural(2) = val2
                    
                    do i = 1,nelem_subgrid
                        
                        val1 = scale_factor*real(i-1,dp) - 1.0_dp
                        val2 = scale_factor*real(i,dp) - 1.0_dp

                        node1_natural(1) = val1
                        node2_natural(1) = val2
                        node3_natural(1) = val1
                        node4_natural(1) = val2
                        
                        call bilinear_node_3d(cell_geom,node1_natural,node1)
                        call bilinear_node_3d(cell_geom,node2_natural,node2)
                        call bilinear_node_3d(cell_geom,node3_natural,node3)
                        call bilinear_node_3d(cell_geom,node4_natural,node4)
                        
                        ! calculate normal vector
                        vec1 = node4 - node1
                        vec2 = node3 - node2

                        ! normal vector
                        normvec(1,i,j) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                        normvec(2,i,j) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                        normvec(3,i,j) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))
                        
                    end do
                end do
                
            case(5) !zhi
                
                node1_natural(3) = 1.0_dp
                node2_natural(3) = 1.0_dp
                node3_natural(3) = 1.0_dp
                node4_natural(3) = 1.0_dp
                
                do j = 1,nelem_subgrid
                    
                    val1 = scale_factor*real(j-1,dp) - 1.0_dp
                    val2 = scale_factor*real(j,dp) - 1.0_dp
                    
                    node1_natural(2) = val1
                    node2_natural(2) = val1
                    node3_natural(2) = val2
                    node4_natural(2) = val2
                    
                    do i = 1,nelem_subgrid

                        val1 = scale_factor*real(i-1,dp) - 1.0_dp
                        val2 = scale_factor*real(i,dp) - 1.0_dp

                        node1_natural(1) = val1
                        node2_natural(1) = val2
                        node3_natural(1) = val1
                        node4_natural(1) = val2
                        
                        call bilinear_node_3d(cell_geom,node1_natural,node1)
                        call bilinear_node_3d(cell_geom,node2_natural,node2)
                        call bilinear_node_3d(cell_geom,node3_natural,node3)
                        call bilinear_node_3d(cell_geom,node4_natural,node4)
                        
                        ! calculate normal vector
                        vec1 = node4 - node1
                        vec2 = node3 - node2

                        ! normal vector
                        normvec(1,i,j) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                        normvec(2,i,j) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                        normvec(3,i,j) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))
                        
                    end do
                    
                end do
                
        end select
        
        
    end subroutine
    
    
end module
