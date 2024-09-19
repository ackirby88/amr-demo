!>
!> \file   solver_residual_subgrid_2d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D subgrid residual
!>
!> Created on August 30, 2018, 9:43 AM
!>

module solver_residual_subgrid_2d_module
    implicit none
    contains

    !> 2D subgrid residual
    !>
    !> @param nfields           Number of fields (unknowns) at each grid point
    !> @param nelem_subgrid     Number of 1D subgrid elements in each AMR quad
    !> @param cell_type         Cell type
    !> @param level             AMR level
    !> @param cell_geom         AMR cell geometry coordinates
    !> @param Q                 Solution pointer for this AMR quad
    !> @param res               Residual pointer for this AMR quad
    !>
    !>
    !>      bilinear mapping
    !>    n3 o-----------o n4
    !>       |           |       n1 = (-1,-1),   N1 = 0.25_dp*(1-xi)*(1-eta)
    !>       |           |       n2 = ( 1,-1),   N2 = 0.25_dp*(1+xi)*(1-eta)
    !>       |           |       n3 = (-1, 1),   N3 = 0.25_dp*(1-xi)*(1+eta)
    !>       |           |       n4 = ( 1, 1),   N4 = 0.25_dp*(1+xi)*(1+eta)
    !>    n1 o-----------o n2
    subroutine solver_residual_subgrid_2d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,res)
        use my_kinddefs
        use flux_2d_module
        use bilinear_node_2d_module
        implicit none

        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: cell_type
        integer(i4),intent(in)      :: level
        real(dp),   intent(in)      :: cell_geom(3,4) !always 3 dimensional
        real(dp),   intent(in)      :: Q(nfields,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid)


        integer(i4) :: i,j

        real(dp)    :: flux(4)
        real(dp)    :: normvec(2)
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

        ! bilinear mapping
        ! n3 o-----------o n4
        !    |           |       n1 = (-1,-1),   N1 = 0.25_dp*(1-xi)*(1-eta)
        !    |           |       n2 = ( 1,-1),   N2 = 0.25_dp*(1+xi)*(1-eta)
        !    |           |       n3 = (-1, 1),   N3 = 0.25_dp*(1-xi)*(1+eta)
        !    |           |       n4 = ( 1, 1),   N4 = 0.25_dp*(1+xi)*(1+eta)
        ! n1 o-----------o n2

        !print*,"node1: ",cell_geom(1,1),cell_geom(2,1)
        !print*,"node2: ",cell_geom(1,2),cell_geom(2,2)
        !print*,"node3: ",cell_geom(1,3),cell_geom(2,3)
        !print*,"node4: ",cell_geom(1,4),cell_geom(2,4)

        !x faces --> n2-n4
        do j = 1,nelem_subgrid
            node1_natural(2) = scale_factor*real(j-1,dp) - 1.0_dp
            node2_natural(2) = scale_factor*real(j,dp) - 1.0_dp

            do i = 1,nelem_subgrid-1
                node1_natural(1) = scale_factor*real(i,dp) - 1.0_dp
                node2_natural(1) = scale_factor*real(i,dp) - 1.0_dp

                call bilinear_node_2d(cell_geom,node1_natural,node1)
                call bilinear_node_2d(cell_geom,node2_natural,node2)

                ! calculate normal vector
                !n = (dy,-dx)
                normvec(1) =  node2(2) - node1(2)
                normvec(2) =  node1(1) - node2(1)

                call flux_roe_2d(Q(:,i,j),Q(:,i+1,j),normvec,flux)
                res(:,i  ,j) = res(:,i  ,j) + flux
                res(:,i+1,j) = res(:,i+1,j) - flux
            end do
        end do

        !y faces --> n4-n3
        do j = 1,nelem_subgrid-1
            node1_natural(2) = scale_factor*real(j,dp) - 1.0_dp
            node2_natural(2) = scale_factor*real(j,dp) - 1.0_dp

            do i = 1,nelem_subgrid
                node1_natural(1) = scale_factor*real(i,dp) - 1.0_dp
                node2_natural(1) = scale_factor*real(i-1,dp) - 1.0_dp

                call bilinear_node_2d(cell_geom,node1_natural,node1)
                call bilinear_node_2d(cell_geom,node2_natural,node2)

                ! calculate normal vector
                !n = (dy,-dx)
                normvec(1) =  node2(2) - node1(2)
                normvec(2) =  node1(1) - node2(1)

                call flux_roe_2d(Q(:,i,j),Q(:,i,j+1),normvec,flux)
                res(:,i,j  ) = res(:,i,j  ) + flux
                res(:,i,j+1) = res(:,i,j+1) - flux

            end do
        end do

        ! divide residual by the volume
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
                res(:,i,j) = res(:,i,j)/volume
            end do
        end do
    end subroutine
end module