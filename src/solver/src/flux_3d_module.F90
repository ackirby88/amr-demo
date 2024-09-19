!     
! \file   flux_3d_module.F90
! \author akirby
!
! Created on August 28, 2018, 10:41 AM
!

module flux_3d_module
    implicit none
    contains

    
    subroutine flux_bc_3d(q,normal,flux)
        use my_kinddefs
        implicit none
        
        real(dp),intent(in) :: q(5)
        real(dp),intent(in) :: normal(3)
        real(dp),intent(out):: flux(5)
        
        real(dp) :: fluxLoc(5,3)

        
        call flux_native_3d(q,fluxLoc)
        
        flux = fluxLoc(:,1)*normal(1) + fluxLoc(:,2)*normal(2) + fluxLoc(:,3)*normal(3)

    end subroutine


    
    subroutine flux_native_3d(q,flux)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: q(5) 
        real(dp),intent(out):: flux(5,3)
        

        real(dp) :: rho,oneOrho,u,v,w,p
        real(dp) :: gamma = 1.4_dp

        
        rho = q(1)
        oneOrho = 1.0_dp/rho
        u   = q(2)*oneOrho
        v   = q(3)*oneOrho
        w   = q(4)*oneOrho
        p   = (gamma-1.0_dp)*(q(5) - 0.5_dp*rho*(u*u + v*v + w*w))

        flux(1,1) = q(2)
        flux(2,1) = q(2)*u + p
        flux(3,1) = q(2)*v
        flux(4,1) = q(2)*w
        flux(5,1) = u*(q(5) + p)

        flux(1,2) = q(3)
        flux(2,2) = q(3)*u
        flux(3,2) = q(3)*v + p
        flux(4,2) = q(3)*w
        flux(5,2) = v*(q(5) + p)
        
        flux(1,3) = q(4)
        flux(2,3) = q(4)*u
        flux(3,3) = q(4)*v
        flux(4,3) = q(4)*w + p
        flux(5,3) = w*(q(5) + p)

    end subroutine


    
    subroutine flux_roe_3d(ql,qr,norm,flux)
        use my_kinddefs
        implicit none

        real(dp),intent(in)  :: ql(5)
        real(dp),intent(in)  :: qr(5)
        real(dp),intent(in)  :: norm(3)
        real(dp),intent(out) :: flux(5)
        
        
        real(dp) :: njk(3)                      ! Normalized normal
        real(dp) :: mag                         ! normal vector magnitude
        real(dp) :: nx,ny,nz                    ! Normal vector
        real(dp) :: ul,ur,vl,vr,wl,wr           ! Velocity component.
        real(dp) :: rhol,rhor,pl,pr             ! Primitive variables
        real(dp) :: oneOrhol,oneOrhor           ! one / density
        real(dp) :: qnl,qnr                     ! Normal velocities
        real(dp) :: al,ar,Hl,Hr                 ! Speed of sound, Total enthalpy
        real(dp) :: RT,rho,u,v,w,H,a,qn         ! Roe-averages
        real(dp) :: drho,dqn,dpress,LdU(4)          ! Wave strengths
        real(dp) :: du,dv,dw                    ! Velocity differences
        real(dp) :: ws(4),R(5,4)                ! Wave speeds and right-eigenvectors
        real(dp) :: dws(4)                      ! Width of a parabolic fit for entropy fix
        real(dp) :: fluxl(5),fluxr(5),diss(5)   ! Fluxes ad dissipation term
        real(dp) :: gamma = 1.4_dp              ! Ratio of specific heats
        real(dp) :: gm1 = 0.4_dp                ! 1 - gamma
        real(dp) :: zero = 0.0_dp,fifth = 0.2_dp


        
        !---> Geometry
        mag = sqrt(norm(1)*norm(1) + norm(2)*norm(2) + norm(3)*norm(3))
        njk = norm/mag
        

        ! Face normal vector (unit vector)
        nx = njk(1)
        ny = njk(2)
        nz = njk(3)


        !Primitive and other variables.
        !  Left state
        rhol        = ql(1)
        oneOrhol    = 1.0_dp/rhol
        
        ul = ql(2)*oneOrhol
        vl = ql(3)*oneOrhol
        wl = ql(4)*oneOrhol
        qnl = ul*nx + vl*ny + wl*nz
        
        pl = (gm1)*(ql(5) - 0.5_dp*rhol*(ul*ul + vl*vl + wl*wl))
        al = sqrt(gamma*pl/rhol)
        Hl = al*al/(gm1) + half*(ul*ul+vl*vl+wl*wl)
        
        
        !  Right state
        rhor        = qr(1)
        oneOrhor    = 1.0_dp/rhor
        
        ur = qr(2)*oneOrhor
        vr = qr(3)*oneOrhor
        wr = qr(4)*oneOrhor
        qnr = ur*nx + vr*ny + wr*nz
        
        pr = (gm1)*(qr(5) - 0.5_dp*rhor*(ur*ur + vr*vr + wr*wr))
        ar = sqrt(gamma*pr/rhor)
        Hr = ar*ar/(gm1) + half*(ur*ur+vr*vr+wr*wr)
        

        !Compute the Roe-averaged quantities
        RT = sqrt(rhor/rhol)
        rho = RT*rhol                                   !Roe-averaged density
        u = (ul + RT*ur)/(one + RT)                     !Roe-averaged x-velocity
        v = (vl + RT*vr)/(one + RT)                     !Roe-averaged y-velocity
        w = (wl + RT*wr)/(one + RT)                     !Roe-averaged z-velocity
        H = (Hl + RT*Hr)/(one + RT)                     !Roe-averaged total enthalpy
        a = sqrt( (gm1)*(H-half*(u*u + v*v + w*w)) )    !Roe-averaged speed of sound
        qn = u*nx + v*ny + w*nz                         !Roe-averaged face-normal velocity


        !Wave Strengths

        drho = rhor - rhol      !Density difference
        dpress   = pr - pl      !Pressure difference
        dqn  = qnr - qnl        !Normal velocity difference

        LdU(1) = (dpress - rho*a*dqn )/(two*a*a)    !Left-moving acoustic wave strength
        LdU(2) =  drho - dpress/(a*a)               !Entropy wave strength
        LdU(3) = (dpress + rho*a*dqn )/(two*a*a)    !Right-moving acoustic wave strength
        LdU(4) = rho                                !Shear wave strength (not really, just a factor)

        
        !Absolute values of the wave Speeds
        ws(1) = abs(qn-a) !Left-moving acoustic wave
        ws(2) = abs(qn)   !Entropy wave
        ws(3) = abs(qn+a) !Right-moving acoustic wave
        ws(4) = abs(qn)   !Shear waves


        !Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
        !NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

        dws(1) = fifth
        if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
        dws(3) = fifth
        if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

        
        !Right Eigenvectors

        ! Left-moving acoustic wave
        R(1,1) = one    
        R(2,1) = u - a*nx
        R(3,1) = v - a*ny
        R(4,1) = w - a*nz
        R(5,1) = H - a*qn

        ! Entropy wave
        R(1,2) = one
        R(2,2) = u
        R(3,2) = v 
        R(4,2) = w
        R(5,2) = half*(u*u + v*v + w*w)

        ! Right-moving acoustic wave
        R(1,3) = one
        R(2,3) = u + a*nx
        R(3,3) = v + a*ny
        R(4,3) = w + a*nz
        R(5,3) = H + a*qn

        ! Two shear wave components combined into one (wave strength incorporated).
        du = ur - ul
        dv = vr - vl
        dw = wr - wl
        R(1,4) = zero
        R(2,4) = du - dqn*nx
        R(3,4) = dv - dqn*ny
        R(4,4) = dw - dqn*nz
        R(5,4) = u*du + v*dv + w*dw - qn*dqn

        !Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]
        diss(:) = ws(1)*LdU(1)*R(:,1) + &
                  ws(2)*LdU(2)*R(:,2) + &
                  ws(3)*LdU(3)*R(:,3) + &
                  ws(4)*LdU(4)*R(:,4)


        !Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)
        fluxl(1) = rhol*qnl
        fluxl(2) = rhol*qnl * ul + pl*nx
        fluxl(3) = rhol*qnl * vl + pl*ny
        fluxl(4) = rhol*qnl * wl + pl*nz
        fluxl(5) = rhol*qnl * Hl

        fluxr(1) = rhor*qnr
        fluxr(2) = rhor*qnr * ur + pr*nx
        fluxr(3) = rhor*qnr * vr + pr*ny
        fluxr(4) = rhor*qnr * wr + pr*nz
        fluxr(5) = rhor*qnr * Hr

        
        ! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]
        flux = half*(fluxl + fluxr - diss)
        flux = flux*mag

        
    end subroutine
    
    
    
end module