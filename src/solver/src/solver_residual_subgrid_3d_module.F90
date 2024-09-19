!>     
!> \file   solver_residual_subgrid_3d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 3D subgrid residual
!>
!> Created on August 30, 2018, 9:44 AM
!>

module solver_residual_subgrid_3d_module
    implicit none
    contains


    !> 3D subgrid residual
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
    !>          bilinear mapping
    !>
    !>        n3 o------------o n4
    !>          /|           /|     n1 = (-1,-1,-1),   N1 = 0.125_dp*(1-xi)*(1-eta)*(1-sigma)
    !>         / |          / |     n2 = ( 1,-1,-1),   N2 = 0.125_dp*(1+xi)*(1-eta)*(1-sigma)
    !>        /  |         /  |     n3 = (-1, 1,-1),   N3 = 0.125_dp*(1-xi)*(1+eta)*(1-sigma)
    !>    n7 o------------o n8|     n4 = ( 1, 1,-1),   N4 = 0.125_dp*(1+xi)*(1+eta)*(1-sigma)
    !>       |n1 o--------|---o n2  n5 = (-1,-1, 1),   N5 = 0.125_dp*(1-xi)*(1-eta)*(1+sigma)
    !>       |  /         |  /      n6 = ( 1,-1, 1),   N6 = 0.125_dp*(1+xi)*(1-eta)*(1+sigma)
    !>       | /          | /       n7 = (-1, 1, 1),   N7 = 0.125_dp*(1-xi)*(1+eta)*(1+sigma)
    !>       |/           |/        n8 = ( 1, 1, 1),   N8 = 0.125_dp*(1+xi)*(1+eta)*(1+sigma)
    !>    n5 o------------o n6
    subroutine solver_residual_subgrid_3d(nfields,nelem_subgrid,cell_type,level,cell_geom,Q,res)
        use my_kinddefs
        use flux_3d_module
        use bilinear_node_3d_module
        use cell_volume_module
        implicit none

        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: cell_type
        integer(i4),intent(in)      :: level
        real(dp),   intent(in)      :: cell_geom(3,8)
        real(dp),   intent(in)      :: Q(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)


        integer(i4) :: i,j,k

        real(dp)    :: flux(5)
        real(dp)    :: normvec(3)
        real(dp)    :: node1(3)
        real(dp)    :: node2(3)
        real(dp)    :: node3(3)
        real(dp)    :: node4(3)
        real(dp)    :: node1_natural(3)
        real(dp)    :: node2_natural(3)
        real(dp)    :: node3_natural(3)
        real(dp)    :: node4_natural(3)
        real(dp)    :: node5_natural(3)
        real(dp)    :: node6_natural(3)
        real(dp)    :: node7_natural(3)
        real(dp)    :: node8_natural(3)
        real(dp)    :: vec1(3)
        real(dp)    :: vec2(3)
        real(dp)    :: xc(3,8)

        real(dp)    :: scale_factor,val1,val2
        real(dp)    :: volume



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
        !    |/           |/        n8 = ( 1, 1, 1),   N8 = 0.125_dp*(1+xi)*(1+eta)*(1+sigma)
        ! n5 o------------o n6



        !x faces --> n2-n4-n6-n8
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

                do i = 1,nelem_subgrid-1

                    val2 = scale_factor*real(i,dp) - 1.0_dp

                    node1_natural(1) = val2
                    node2_natural(1) = val2
                    node3_natural(1) = val2
                    node4_natural(1) = val2

                    call bilinear_node_3d(cell_geom,node1_natural,node1)
                    call bilinear_node_3d(cell_geom,node2_natural,node2)
                    call bilinear_node_3d(cell_geom,node3_natural,node3)
                    call bilinear_node_3d(cell_geom,node4_natural,node4)

                    ! calculate normal vector
                    ! cross vector 1: n2->n8
                    vec1 = node4 - node1

                    ! cross vector 2: n4->n6
                    vec2 = node3 - node2

                    ! normal vector
                    normvec(1) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                    normvec(2) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                    normvec(3) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))


                    !print*,"node[1]: ",node1
                    !print*,"node[2]: ",node2
                    !print*,"node[3]: ",node3
                    !print*,"node[4]: ",node4
                    !print*,"normal:  ",normvec,sqrt(normvec(1)*normvec(1) + normvec(2)*normvec(2) + normvec(3)*normvec(3))


                    !print*,"elem[",i,",",j,"]:",normvec(1),normvec(2),normvec(3)
                    call flux_roe_3d(Q(:,i,j,k),Q(:,i+1,j,k),normvec,flux)
                    res(:,i  ,j,k) = res(:,i  ,j,k) + flux
                    res(:,i+1,j,k) = res(:,i+1,j,k) - flux

                end do

            end do

        end do


        !y faces --> n3-n4-n7-n8
        do k = 1,nelem_subgrid

            val1 = scale_factor*real(k-1,dp) - 1.0_dp
            val2 = scale_factor*real(k,dp) - 1.0_dp

            node1_natural(3) = val1
            node2_natural(3) = val1
            node3_natural(3) = val2
            node4_natural(3) = val2

            do j = 1,nelem_subgrid-1

                val2 = scale_factor*real(j,dp) - 1.0_dp

                node1_natural(2) = val2
                node2_natural(2) = val2
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
                    ! cross vector 1: n4->n7
                    vec1 = node3 - node2
                    
                    ! cross vector 2: n3->n8
                    vec2 = node4 - node1

                    ! normal vector
                    normvec(1) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                    normvec(2) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                    normvec(3) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))


                    !print*,"node[1]: ",node1
                    !print*,"node[2]: ",node2
                    !print*,"node[3]: ",node3
                    !print*,"node[4]: ",node4
                    !print*,"normal:  ",normvec,sqrt(normvec(1)*normvec(1) + normvec(2)*normvec(2) + normvec(3)*normvec(3))


                    !print*,"elem[",i,",",j,"]:",normvec(1),normvec(2),normvec(3)
                    call flux_roe_3d(Q(:,i,j,k),Q(:,i,j+1,k),normvec,flux)
                    res(:,i,j  ,k) = res(:,i,j  ,k) + flux
                    res(:,i,j+1,k) = res(:,i,j+1,k) - flux

                end do

            end do

        end do


        !z faces --> n5-n6-n7-n8
        do k = 1,nelem_subgrid-1

            val2 = scale_factor*real(k,dp) - 1.0_dp

            node1_natural(3) = val2
            node2_natural(3) = val2
            node3_natural(3) = val2
            node4_natural(3) = val2


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
                    ! cross vector 1: n5->n8
                    vec1 = node4 - node1

                    ! cross vector 2: n6->n7
                    vec2 = node3 - node2

                    ! normal vector
                    normvec(1) = 0.5_dp*(vec1(2)*vec2(3) - vec1(3)*vec2(2))
                    normvec(2) = 0.5_dp*(vec1(3)*vec2(1) - vec1(1)*vec2(3))
                    normvec(3) = 0.5_dp*(vec1(1)*vec2(2) - vec1(2)*vec2(1))


                    !print*,"node[1]: ",node1
                    !print*,"node[2]: ",node2
                    !print*,"node[3]: ",node3
                    !print*,"node[4]: ",node4
                    !print*,"normal:  ",normvec,sqrt(normvec(1)*normvec(1) + normvec(2)*normvec(2) + normvec(3)*normvec(3))


                    !print*,"elem[",i,",",j,"]:",normvec(1),normvec(2),normvec(3)
                    call flux_roe_3d(Q(:,i,j,k),Q(:,i,j,k+1),normvec,flux)
                    res(:,i,j,k  ) = res(:,i,j,k  ) + flux
                    res(:,i,j,k+1) = res(:,i,j,k+1) - flux

                end do

            end do

        end do



        ! divide residual by the volume
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
                    res(:,i,j,k) = res(:,i,j,k)/volume

                end do

            end do

        end do


    end subroutine


end module
