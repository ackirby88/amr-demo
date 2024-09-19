
!     
! \file   bc_3d_module.F90
! \author akirby
!
! Created on August 31, 2018, 12:32 PM
!


module bc_3d_module
    implicit none
    contains
    
    
    subroutine bc_3d(bc_type,normal,Q,Q_bc)
        use my_kinddefs
        implicit none
        
        integer(i4),intent(in) :: bc_type
        real(dp),   intent(in) :: normal(3)
        real(dp),   intent(in) :: Q(5)
        real(dp),   intent(out):: Q_bc(5)
        
        
        real(dp) :: mag
        real(dp) :: nx,ny,nz
        real(dp) :: rho,oneOrho,u,v,w,p
        real(dp) :: ub,vb,wb
        real(dp) :: U_dot_n
        real(dp) :: gamma = 1.4_dp
        real(dp) :: gm1
        
        
        gm1 = gamma - 1.0_dp

        
        rho = q(1)
        oneOrho = 1.0_dp/rho
        u   = q(2)*oneOrho
        v   = q(3)*oneOrho
        w   = q(4)*oneOrho
        p   = gm1*(q(5) - 0.5_dp*rho*(u*u + v*v + w*w))
        
        
        
        select case(bc_type)
            
            case(0) !characteristic
                
                Q_bc = Q
                
            case(1) !slip wall
                
                mag = sqrt(normal(1)*normal(1) + normal(2)*normal(2) + normal(3)*normal(3))
                nx = normal(1)/mag
                ny = normal(2)/mag
                nz = normal(3)/mag
                
                U_dot_n = u*nx + v*ny + w*nz
                
                ub = u - U_dot_n*nx
                vb = v - U_dot_n*ny
                wb = w - U_dot_n*nz
                
                Q_bc(1) = Q(1)
                Q_bc(2) = Q(1)*ub
                Q_bc(3) = Q(1)*vb
                Q_bc(4) = Q(1)*wb
                Q_bc(5) = Q(5)
                
            case(2) !no slip wall
                
                Q_bc    = 0.0_dp
                Q_bc(1) = Q(1)
                Q_bc(5) = Q(5) - 0.5_dp*rho*(u*u + v*v + w*w)
                
        end select
        
        
    end subroutine
    
    
end module
