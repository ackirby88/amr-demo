!>     
!> \file   solver_time_step_3d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 3D subgrid residual
!>
!> Created on September 19, 2018, 12:50 PM
!>

module solver_time_step_3d_module
    implicit none
    contains

    
    !> 3D time step calculation
    !>
    !> @param nfields           Number of fields (unknowns) at each grid point
    !> @param nelem_subgrid     Number of 1D subgrid elements in each AMR quad
    !> @param cell_type         Cell type
    !> @param level             AMR level
    !> @param cell_geom         AMR cell geometry coordinates
    !> @param Q                 Solution pointer for this AMR quad
    !> @param cfl               CFL user input
    !> @param dt                Maximum stable time step for this mpi rank
    !>
    subroutine solver_time_step_3d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,cfl,dt)
        use my_kinddefs
        use cell_volume_module
        use bilinear_node_3d_module
        implicit none

        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: cell_type
        integer(i4),intent(in)      :: level
        real(dp),   intent(in)      :: cell_geom(3,8) !always 3 dimensional
        real(dp),   intent(in)      :: Q(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(in)      :: cfl
        real(dp),   intent(inout)   :: dt


        real(dp)    :: gamma = 1.4_dp
        real(dp)    :: gm1 = 0.4_dp
        
        integer(i4) :: i,j,k

        real(dp)    :: u,v,w,c
        real(dp)    :: V2
        real(dp)    :: oneOrho
        real(dp)    :: pOverRho
        real(dp)    :: eig
        
        real(dp)    :: node1_natural(3)
        real(dp)    :: node2_natural(3)
        real(dp)    :: node3_natural(3)
        real(dp)    :: node4_natural(3)
        real(dp)    :: node5_natural(3)
        real(dp)    :: node6_natural(3)
        real(dp)    :: node7_natural(3)
        real(dp)    :: node8_natural(3)
        real(dp)    :: xc(3,8)

        real(dp)    :: scale_factor,val1,val2
        real(dp)    :: volume

        
        scale_factor = 2.0_dp / real(nelem_subgrid,dp)

        
        do k = 1,nelem_subgrid
            val1 = scale_factor*real(k-1,dp) - 1.0_dp
            val2 = scale_factor*real(k,dp) - 1.0_dp

            node1_natural(3) = val1
            node2_natural(3) = val1
            node3_natural(3) = val1
            node4_natural(3) = val1
            node5_natural(3) = val2
            node6_natural(3) = val2
            node7_natural(3) = val2
            node8_natural(3) = val2

            do j = 1,nelem_subgrid

                val1 = scale_factor*real(j-1,dp) - 1.0_dp
                val2 = scale_factor*real(j,dp) - 1.0_dp

                node1_natural(2) = val1
                node2_natural(2) = val1
                node3_natural(2) = val2
                node4_natural(2) = val2
                node5_natural(2) = val1
                node6_natural(2) = val1
                node7_natural(2) = val2
                node8_natural(2) = val2

                do i = 1,nelem_subgrid

                    val1 = scale_factor*real(i-1,dp) - 1.0_dp
                    val2 = scale_factor*real(i,dp) - 1.0_dp

                    node1_natural(1) = val1
                    node2_natural(1) = val2
                    node3_natural(1) = val1
                    node4_natural(1) = val2
                    node5_natural(1) = val1
                    node6_natural(1) = val2
                    node7_natural(1) = val1
                    node8_natural(1) = val2

                    call bilinear_node_3d(cell_geom,node1_natural,xc(:,1))
                    call bilinear_node_3d(cell_geom,node2_natural,xc(:,2))
                    call bilinear_node_3d(cell_geom,node3_natural,xc(:,3))
                    call bilinear_node_3d(cell_geom,node4_natural,xc(:,4))
                    call bilinear_node_3d(cell_geom,node5_natural,xc(:,5))
                    call bilinear_node_3d(cell_geom,node6_natural,xc(:,6))
                    call bilinear_node_3d(cell_geom,node7_natural,xc(:,7))
                    call bilinear_node_3d(cell_geom,node8_natural,xc(:,8))

                    call cell_volume(xc,volume)
                    
                    oneOrho = 1._dp/Q(1,i,j,k)

                    u = Q(2,i,j,k)*oneOrho
                    v = Q(3,i,j,k)*oneOrho
                    w = Q(4,i,j,k)*oneOrho

                    V2 = u*u + v*v + w*w

                    pOverRho = gm1*(Q(5,i,j,k)*oneOrho - half*V2)
                    c = sqrt(gamma * pOverRho)

                    eig = sqrt(V2) + c
                    dt = min(dt,cfl*(volume**(1.0_dp/3.0_dp))/eig)
                    

                end do

            end do

        end do
        
        
    end subroutine
    
    
end module
