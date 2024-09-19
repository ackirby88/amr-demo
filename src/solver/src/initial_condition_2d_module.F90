
!
! \file   initial_condition_2d_module.F90
! \author akirby
!
! Created on August 31, 2018, 11:09 AM
!
module initial_condition_2d_module
    implicit none
    contains

    subroutine initial_condition_2d(initial_condition,unstructured,nfields,nelem_subgrid,h,cell_geom,Q,&
                                      mach,alpha,beta,gamma,density,pressure)
        use my_kinddefs
        use bilinear_node_2d_module
        implicit none

        integer(i4),intent(in)      :: initial_condition
        integer(i4),intent(in)      :: unstructured
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        real(dp),   intent(in)      :: h
        real(dp),   intent(in)      :: cell_geom(3,4) !always 3 dimensional
        real(dp),   intent(out)     :: Q(nfields,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(in)      :: mach
        real(dp),   intent(in)      :: alpha
        real(dp),   intent(in)      :: beta
        real(dp),   intent(in)      :: gamma
        real(dp),   intent(in)      :: density
        real(dp),   intent(in)      :: pressure

        real(dp)    :: xo,yo
        real(dp)    :: rho,u,v,p
        real(dp)    :: ff,a,two_a,beta_L
        real(dp)    :: radius2,dT1,dT

        real(dp)    :: gm1,oneOgm1,oneOrho0
        real(dp)    :: RT0,oneORT0
        real(dp)    :: RC,RT0pRC
        real(dp)    :: C_0,V_0
        real(dp)    :: u_in,v_in
        real(dp)    :: rhou,rhov,rhoe
        real(dp)    :: time

        real(dp)    :: node_natural(2),node(2)
        real(dp)    :: scale_factor
        real(dp)    :: val1,val2

        integer(i4) :: i,j

        ! flow physics
        gm1         = gamma - 1.0_dp
        oneOgm1     = 1.0_dp / gm1
        RT0         = pressure / density
        oneORT0     = 1.0_dp / RT0
        RC          = 120.0_dp/291.15_dp * RT0
        RT0pRC      = RT0 + RC

        C_0         = sqrt(gamma * RT0)
        V_0         = mach * C_0

        oneOrho0    = 1.0_dp / density
        u_in        = V_0*cos(alpha*pi/180.0_dp)*cos(beta*pi/180.0_dp)
        v_in        = V_0*sin(alpha*pi/180.0_dp)*cos(beta*pi/180.0_dp)

        rhou = density*u_in
        rhov = density*v_in
        rhoe = (pressure * oneOgm1)  + half*density*(u_in*u_in + v_in*v_in)

        time = 0.0_dp

        select case(initial_condition)
            case(0) !free stream
                do j = 1,nelem_subgrid
                do i = 1,nelem_subgrid
                    Q(1,i,j) = density
                    Q(2,i,j) = rhou
                    Q(3,i,j) = rhov
                    Q(4,i,j) = rhoe
                end do
                end do

            case(1) !isentropic vortex
                xo = 0.0_dp
                yo = 0.0_dp

                u  = u_in
                v  = v_in

                xo = xo + u*time
                yo = yo + v*time

                a       = 0.5_dp
                beta_L  = 3.0_dp
                two_a   = 2.0_dp*a
                dT1     = -gm1*beta_L*beta_L/(16.0_dp*a*gamma*pi*pi)

                scale_factor = 2.0_dp / real(nelem_subgrid,dp)

                do j = 1,nelem_subgrid
                    val1 = scale_factor*real(j-1,dp) - 1.0_dp
                    val2 = scale_factor*real(j,dp) - 1.0_dp
                    node_natural(2) = 0.5_dp*(val1+val2)

                    do i = 1,nelem_subgrid
                        val1 = scale_factor*real(i-1,dp) - 1.0_dp
                        val2 = scale_factor*real(i,dp) - 1.0_dp
                        node_natural(1) = 0.5_dp*(val1+val2)

                        call bilinear_node_2d(cell_geom,node_natural,node)

                        u  = u_in
                        v  = v_in

                        radius2 = (node(1)- xo)*(node(1) - xo) + (node(2) - yo)*(node(2) - yo)
                        dT = dT1*exp(two_a*(1.0_dp - radius2))
                        rho = (1._dp + dT)**(oneOgm1)
                        ff = beta_L/(twopi)*exp(a*(1.0_dp - radius2))

                        u = u - (node(2) - yo)*ff
                        v = v + (node(1) - xo)*ff
                        p = pressure + rho**(gamma) - 1.0_dp

                        Q(1,i,j) = rho
                        Q(2,i,j) = rho*u
                        Q(3,i,j) = rho*v
                        Q(4,i,j) = p*oneOgm1 + half*rho*(u*u+v*v)
                    end do
                end do
        end select
    end subroutine
end module