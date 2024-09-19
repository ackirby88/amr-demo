
!     
! \file   normals_2d_module.F90
! \author akirby
!
! Created on August 28, 2018, 3:14 PM
!


module normals_2d_module
    implicit none
    contains
    
    
    subroutine normals_2d(nelem_subgrid,side,cell_geom,normvec)
        use my_kinddefs
        use bilinear_node_2d_module
        implicit none
        
        integer(i4),intent(in) :: nelem_subgrid
        integer(i4),intent(in) :: side
        real(dp),   intent(in) :: cell_geom(3,4)
        real(dp),   intent(out):: normvec(2,nelem_subgrid)

        
        integer(i4) :: i,j
        
        real(dp)    :: scale_factor
        real(dp)    :: node1(2)
        real(dp)    :: node2(2)    
        real(dp)    :: node1_natural(2)
        real(dp)    :: node2_natural(2)
        
        scale_factor = 2.0_dp / real(nelem_subgrid,dp)
        
        
        
        ! bilinear mapping
        ! n3 o-----------o n4     
        !    |           |       n1 = (-1,-1),   N1 = 0.25_dp*(1-xi)*(1-eta)
        !    |           |       n2 = ( 1,-1),   N2 = 0.25_dp*(1+xi)*(1-eta)
        !    |           |       n3 = (-1, 1),   N3 = 0.25_dp*(1-xi)*(1+eta)
        !    |           |       n4 = ( 1, 1),   N4 = 0.25_dp*(1+xi)*(1+eta)
        ! n1 o-----------o n2
        
        
        
        select case(side)
            
            case(0) !xlo
                
                node1_natural(1) = -1.0_dp
                node2_natural(1) = -1.0_dp
                    
                do j = 1,nelem_subgrid
                    
                    node1_natural(2) = scale_factor*real(j,dp) - 1.0_dp
                    node2_natural(2) = scale_factor*real(j-1,dp) - 1.0_dp
                    
                    call bilinear_node_2d(cell_geom,node1_natural,node1)
                    call bilinear_node_2d(cell_geom,node2_natural,node2)
                    
                    normvec(1,j) = node2(2) - node1(2)  ! dy
                    normvec(2,j) = node1(1) - node2(1)  !-dx

                end do
            
            case(1) !xhi
                
                node1_natural(1) = 1.0_dp
                node2_natural(1) = 1.0_dp
                    
                do j = 1,nelem_subgrid
                    
                    node1_natural(2) = scale_factor*real(j-1,dp) - 1.0_dp
                    node2_natural(2) = scale_factor*real(j,dp) - 1.0_dp
                    
                    call bilinear_node_2d(cell_geom,node1_natural,node1)
                    call bilinear_node_2d(cell_geom,node2_natural,node2)
                    
                    normvec(1,j) = node2(2) - node1(2)  ! dy
                    normvec(2,j) = node1(1) - node2(1)  !-dx

                end do
                
            case(2) !ylo
                
                node1_natural(2) = -1.0_dp 
                node2_natural(2) = -1.0_dp 
                    
                do i = 1,nelem_subgrid
                    
                    node1_natural(1) = scale_factor*real(i-1,dp) - 1.0_dp
                    node2_natural(1) = scale_factor*real(i,dp) - 1.0_dp    
                    
                    call bilinear_node_2d(cell_geom,node1_natural,node1)
                    call bilinear_node_2d(cell_geom,node2_natural,node2)
                    
                    normvec(1,i) = node2(2) - node1(2)  ! dy
                    normvec(2,i) = node1(1) - node2(1)  !-dx
                    
                end do
                
            case(3) !yhi
                
                node1_natural(2) = 1.0_dp 
                node2_natural(2) = 1.0_dp
                    
                do i = 1,nelem_subgrid
                    
                    node1_natural(1) = scale_factor*real(i,dp) - 1.0_dp
                    node2_natural(1) = scale_factor*real(i-1,dp) - 1.0_dp
                    
                    call bilinear_node_2d(cell_geom,node1_natural,node1)
                    call bilinear_node_2d(cell_geom,node2_natural,node2)
                    
                    normvec(1,i) = node2(2) - node1(2)  ! dy
                    normvec(2,i) = node1(1) - node2(1)  !-dx
                    
                end do
                
        end select
        

    end subroutine
    
    
end module
