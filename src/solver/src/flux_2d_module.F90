!     
! \file   flux_module.F90
! \author akirby
!
! Created on August 27, 2018, 3:20 PM
!

module flux_2d_module
    implicit none
    contains

    
    subroutine flux_bc_2d(q,normal,flux)
        use my_kinddefs
        implicit none
        
        real(dp),intent(in) :: q(4)
        real(dp),intent(in) :: normal(2)
        real(dp),intent(out):: flux(4)
        
        
        real(dp) :: fluxLoc(4,2)

        
        call flux_native_2d(q,fluxLoc)
        
        flux = fluxLoc(:,1)*normal(1) + fluxLoc(:,2)*normal(2)

    end subroutine


    subroutine flux_native_2d(q,flux)
        use my_kinddefs
        implicit none

        real(dp),intent(in) :: q(4) 
        real(dp),intent(out):: flux(4,2)
        

        real(dp) :: rho,oneOrho,u,v,p
        real(dp) :: gamma = 1.4_dp

        
        rho = q(1)
        oneOrho = 1.0_dp/rho
        u   = q(2)*oneOrho
        v   = q(3)*oneOrho
        p   = (gamma-1.0_dp)*(q(4) - 0.5_dp*rho*(u*u + v*v))

        flux(1,1) = q(2)
        flux(2,1) = q(2)*u + p
        flux(3,1) = q(2)*v
        flux(4,1) = u*(q(4) + p)

        flux(1,2) = q(3)
        flux(2,2) = q(3)*u
        flux(3,2) = q(3)*v + p
        flux(4,2) = v*(q(4) + p)

    end subroutine


    subroutine flux_roe_2d(ql,qr,norm,flux)
        use my_kinddefs
        implicit none

        real(dp),intent(in)  :: ql(4)
        real(dp),intent(in)  :: qr(4)
        real(dp),intent(in)  :: norm(2)
        real(dp),intent(out) :: flux(4)
        

        !---> Local Variables
        real(dp) :: rhol,oneOrhol,ul,vl,El,Hl,pl,unl
        real(dp) :: rhor,oneOrhor,ur,vr,Er,Hr,pr,unr
        real(dp) :: rhobar,ubar,vbar,Hbar,abar,unbar,Mnbar
        real(dp) :: nx,ny,mag
        real(dp) :: roe_fact,oneOroe_fact
        real(dp) :: b1, b2
        real(dp) :: fluxl,fluxr,disp,bcoeff,Mcoeff,delta1,eigp,eigm
        real(dp) :: eigp_max,eigm_min
        real(dp) :: gamma = 1.4_dp
        real(dp) :: gm1 = 0.4_dp


        
        !---> Geometry
        mag = sqrt(norm(1)*norm(1) + norm(2)*norm(2))
        nx = norm(1)/mag
        ny = norm(2)/mag

        !---> Left Primatives
        rhol = ql(1)
        oneOrhol = 1.0_dp/rhol
        ul = ql(2)*oneOrhol
        vl = ql(3)*oneOrhol
        El = ql(4)*oneOrhol
        pl = (gamma - 1._dp)*(ql(4) - half*rhol*(ul*ul + vl*vl) ) 
        Hl = El + pl*oneOrhol
        unl = ul*nx + vl*ny

        !---> Right Primatives
        rhor = qr(1)
        oneOrhor = 1.0_dp/rhor
        ur = qr(2)*oneOrhor
        vr = qr(3)*oneOrhor
        Er = qr(4)*oneOrhor
        pr = (gamma - 1._dp)*(qr(4) - half*rhor*(ur*ur + vr*vr) ) 
        Hr = Er + pr*oneOrhor
        unr = ur*nx + vr*ny

        !---> Roe Averages
        roe_fact = sqrt(rhol) + sqrt(rhor)
        oneOroe_fact = 1.0_dp / roe_fact
        rhobar = sqrt(rhol*rhor)
        ubar = (ul*sqrt(rhol) + ur*sqrt(rhor))*oneOroe_fact
        vbar = (vl*sqrt(rhol) + vr*sqrt(rhor))*oneOroe_fact
        Hbar = (Hl*sqrt(rhol) + Hr*sqrt(rhor))*oneOroe_fact
        abar = sqrt(gm1*(Hbar - half*(ubar*ubar + vbar*vbar)))
        unbar = ubar*nx + vbar*ny
        Mnbar = unbar/abar

        !---> Coefficients b1 and b2
        eigp = unbar + abar
        eigm = unbar - abar

        eigp_max = unr + abar
        eigm_min = unl - abar
        !   b1 
        if( eigp > 0._dp .or. eigp_max > 0._dp) then
           b1 = max(eigp, eigp_max)
        else
           b1 = 0._dp
        end if

        !   b2 
        if( eigm < 0._dp .or. eigm_min < 0._dp) then 
           b2 = min(eigm, eigm_min)
        else
           b2 = 0._dp
        end if

        !---> Some extra coefficients
        bcoeff = b1*b2/(b1 - b2)
        delta1 = (rhor - rhol) - (pr - pl)/(abar*abar)
        Mcoeff = 1._dp/(1._dp + abs(Mnbar))

        !------------------------ Assemble the flux -------------------------------
        !   Row 1
        fluxl = ql(2)*nx + ql(3)*ny
        fluxr = qr(2)*nx + qr(3)*ny
        disp = bcoeff*(rhor - rhol) - bcoeff*Mcoeff*delta1

        flux(1) = (b1*fluxl - b2*fluxr)/(b1 - b2) + disp

        !   Row 2
        fluxl = (ql(2)*ul + pl)*nx + (ql(3)*ul)*ny
        fluxr = (qr(2)*ur + pr)*nx + (qr(3)*ur)*ny
        disp = bcoeff*(qr(2) - ql(2)) - bcoeff*Mcoeff*( delta1*ubar + &
               rhobar*(ur - ul  - nx*(unr - unl)) ) 

        flux(2) = (b1*fluxl - b2*fluxr)/(b1 - b2) + disp

        !   Row 3
        fluxl = (ql(2)*vl)*nx + (ql(3)*vl + pl)*ny
        fluxr = (qr(2)*vr)*nx + (qr(3)*vr + pr)*ny
        disp = bcoeff*(qr(3) - ql(3)) - bcoeff*Mcoeff*( delta1*vbar + &
               rhobar*(vr - vl - ny*(unr - unl)) )

        flux(3) = (b1*fluxl - b2*fluxr)/(b1 - b2) + disp

        !   Row 4
        fluxl = (ql(4) + pl)*unl
        fluxr = (qr(4) + pr)*unr  
        disp = bcoeff*(rhor*Hr - rhol*Hl) - bcoeff*Mcoeff*( delta1*Hbar + &
               rhobar*(Hr - Hl) ) 

        flux(4) = (b1*fluxl - b2*fluxr)/(b1 - b2) + disp

        flux = flux*mag

    end subroutine
    
    
end module