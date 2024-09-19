
!>     
!> \file   solver_residual_bc_unstructured_2d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D boundary condition residual unstructured grid
!>
!> Created on August 31, 2018, 12:07 PM
!>

module solver_residual_bc_unstructured_2d_module
    implicit none
    contains
    
    
    !> 2D boundary condition residual unstructured grid
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
    subroutine solver_residual_bc_unstructured_2d(nfields,nelem_subgrid,side,cell_type,bc,geom,Q,res)
        use my_kinddefs
        use bc_2d_module
        use flux_2d_module
        use normals_2d_module
        use integrate_2d_module
        use subgrid_index_2d_module
        implicit none

        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: side
        integer(i4),intent(in)      :: cell_type
        integer(i4),intent(in)      :: bc
        real(dp),   intent(in)      :: geom(3,4)
        real(dp),   intent(in)      :: Q(nfields,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(inout)   :: res(nfields,nelem_subgrid,nelem_subgrid)


        real(dp)    :: normal_vec(2,nelem_subgrid)
        real(dp)    :: flux(4)
        real(dp)    :: Q_bc(4)
        integer(i4) :: i,j
        integer(i4) :: n
        
        
        call normals_2d(nelem_subgrid,side,geom,normal_vec)
        
        
        do n = 1,nelem_subgrid
            
            call subgrid_index_2d(n,nelem_subgrid,side,i,j)
            call bc_2d(bc,normal_vec(:,n),Q(:,i,j),Q_bc)
            call flux_bc_2d(Q_bc,normal_vec(:,n),flux)
            call integrate_2d_pos(i,j,nfields,nelem_subgrid,flux,res)
        
        end do
        
        
    end subroutine
    
    
end module
