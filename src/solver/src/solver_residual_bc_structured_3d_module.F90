
!>     
!> \file   solver_residual_bc_structured_3d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 3D boundary condition residual
!>
!> Created on September 4, 2018, 10:10 AM
!>

module solver_residual_bc_structured_3d_module
    implicit none
    contains
    

    !> 3D boundary condition residual
    !> 
    !> @param nfields           Number of fields (unknowns) per grid point
    !> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
    !> @param side              Side id (face number) of AMR quad
    !> @param cell_type         Cell type
    !> @param level             AMR level
    !> @param bc                Boundary condition type
    !> @param h_base            Coarsest level element width (structured)
    !> @param xlo               x-coordinate of first node in AMR quad
    !> @param ylo               y-coordinate of first node in AMR quad
    !> @param zlo               z-coordinate of first node in AMR quad
    !> @param soln              Solution pointer for this quad
    !> @param res               Residual pointer for this quad
    !>
    subroutine solver_residual_bc_structured_3d(nfields,nelem_subgrid,side,cell_type,level,bc,h_base,xlo,ylo,zlo,soln,res)
        use my_kinddefs
        use bc_3d_module
        use flux_3d_module
        use normals_3d_module
        use integrate_3d_module
        use subgrid_index_3d_module
        implicit none
        
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: side
        integer(i4),intent(in)      :: cell_type
        integer(i4),intent(in)      :: level
        integer(i4),intent(in)      :: bc
        real(dp),   intent(in)      :: h_base
        real(dp),   intent(in)      :: xlo
        real(dp),   intent(in)      :: ylo
        real(dp),   intent(in)      :: zlo
        real(dp),   intent(in)      :: soln(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
        
        
        real(dp)    :: normal_vec(3)
        real(dp)    :: Q_bc(5)
        real(dp)    :: flux(5)
        real(dp)    :: cell_h
        real(dp)    :: face_area
        integer(i4) :: i,j,k
        integer(i4) :: n,m
        
        
        cell_h = h_base * (1.0_dp / (2.0_dp**(level))) / real(nelem_subgrid,dp)
        face_area = cell_h*cell_h
        
        
        normal_vec = 0.0_dp
        select case(side)
            case(0)
                normal_vec(1) = -face_area
            case(1)
                normal_vec(1) =  face_area
            case(2)
                normal_vec(2) = -face_area
            case(3)
                normal_vec(2) =  face_area
            case(4)
                normal_vec(3) = -face_area
            case(5)
                normal_vec(3) =  face_area
        end select
        
        
        do m = 1,nelem_subgrid
        do n = 1,nelem_subgrid
            
            call subgrid_index_3d(n,m,nelem_subgrid,side,i,j,k)
            call bc_3d(bc,normal_vec,soln(:,i,j,k),Q_bc)
            call flux_bc_3d(Q_bc,normal_vec,flux)
            call integrate_3d_pos(i,j,k,nfields,nelem_subgrid,flux,res)
        
        end do
        end do
        
        
    end subroutine

    
end module
