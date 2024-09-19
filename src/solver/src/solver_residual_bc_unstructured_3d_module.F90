
!>     
!> \file   solver_residual_bc_unstructured_3d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 3D boundary condition residual unstructured grid
!>
!> Created on August 31, 2018, 12:09 PM
!>

module solver_residual_bc_unstructured_3d_module
    implicit none
    contains
    
    
    !> 3D boundary condition residual unstructured grid
    !> 
    !> @param nfields           Number of fields (unknowns) per grid point
    !> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
    !> @param side              Side id (face number) of AMR quad
    !> @param cell_type         Cell type 
    !> @param bc                Boundary condition type
    !> @param geom              Cell geometry coordinates
    !> @param Q                 Solution pointer for this quad
    !> @param res               Residual pointer for this quad
    !>
    subroutine solver_residual_bc_unstructured_3d(nfields,nelem_subgrid,side,cell_type,bc,geom,Q,res)
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
        integer(i4),intent(in)      :: bc
        real(dp),   intent(in)      :: geom(3,8)
        real(dp),   intent(in)      :: Q(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)

        
        real(dp)    :: normal_vec(3,nelem_subgrid,nelem_subgrid)
        real(dp)    :: flux(5)
        real(dp)    :: Q_bc(5)
        integer(i4) :: i,j,k
        integer(i4) :: n,m
        
        
        call normals_3d(nelem_subgrid,side,geom,normal_vec)
        
        
        do m = 1,nelem_subgrid
        do n = 1,nelem_subgrid
            
            call subgrid_index_3d(n,m,nelem_subgrid,side,i,j,k)
            call bc_3d(bc,normal_vec(:,n,m),Q(:,i,j,k),Q_bc)
            call flux_bc_3d(Q_bc,normal_vec(:,n,m),flux)
            call integrate_3d_pos(i,j,k,nfields,nelem_subgrid,flux,res)
        
        end do
        end do
        
        
    end subroutine
    
    
end module
