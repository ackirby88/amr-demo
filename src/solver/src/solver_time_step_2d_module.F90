!>     
!> \file   solver_time_step_2d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D subgrid residual
!>
!> Created on September 19, 2018, 11:46 AM
!>

module solver_time_step_2d_module
    implicit none
    contains

    
    !> 2D time step calculation
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
    subroutine solver_time_step_2d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,cfl,dt)
        use my_kinddefs
        use bilinear_node_2d_module
        implicit none

        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: cell_type
        integer(i4),intent(in)      :: level
        real(dp),   intent(in)      :: cell_geom(3,4) !always 3 dimensional
        real(dp),   intent(in)      :: Q(nfields,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(in)      :: cfl
        real(dp),   intent(inout)   :: dt


        real(dp)    :: gamma = 1.4_dp
        real(dp)    :: gm1 = 0.4_dp
        
        integer(i4) :: i,j

        real(dp)    :: u,v,c
        real(dp)    :: V2
        real(dp)    :: oneOrho
        real(dp)    :: pOverRho
        real(dp)    :: eig
        
        real(dp)    :: node1(2)
        real(dp)    :: node2(2)
        real(dp)    :: node3(2)
        real(dp)    :: node4(2)
        real(dp)    :: node1_natural(2)
        real(dp)    :: node2_natural(2)
        real(dp)    :: node3_natural(2)
        real(dp)    :: node4_natural(2)
        real(dp)    :: vec1(2)
        real(dp)    :: vec2(2)

        real(dp)    :: scale_factor
        real(dp)    :: val1,val2
        real(dp)    :: volume


        scale_factor = 2.0_dp / real(nelem_subgrid,dp)

        
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

                call bilinear_node_2d(cell_geom,node1_natural,node1)
                call bilinear_node_2d(cell_geom,node2_natural,node2)
                call bilinear_node_2d(cell_geom,node3_natural,node3)
                call bilinear_node_2d(cell_geom,node4_natural,node4)

                ! calculate volume
                ! cross vector 1: n2->n8
                vec1 = node4 - node1

                ! cross vector 2: n4->n6
                vec2 = node3 - node2

                !volume
                volume = 0.5_dp*sqrt((vec1(1)*vec2(2) - vec1(2)*vec2(1))**2)
                
                !fluid properties
                oneOrho = 1._dp/Q(1,i,j)
                u = Q(2,i,j)*oneOrho
                v = Q(3,i,j)*oneOrho
                V2 = u*u + v*v

                pOverRho = gm1*(Q(4,i,j)*oneOrho - half*V2)
                c = sqrt(gamma * pOverRho)
                eig = sqrt(V2) + c
                
                dt = min(dt,cfl*sqrt(volume)/eig)

            end do

        end do
        
    end subroutine
    
    
end module
