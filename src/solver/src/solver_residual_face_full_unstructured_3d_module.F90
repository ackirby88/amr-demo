!>     
!> \file   solver_residual_face_full_unstructured_3d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 3D full face residual unstructured grid
!>
!> Created on August 31, 2018, 10:26 AM
!>


module solver_residual_face_full_unstructured_3d_module
    implicit none
    contains


    !> 3D full face residual unstructured grid
    !>
    !> @param nfields          Number of fields (unknowns) per grid point
    !> @param nelem_subgrid    Number of 1D subgrid elements in AMR quad
    !> @param side_l           Side id (face) of AMR left quad 
    !> @param side_r           Side id (face) of AMR right quad 
    !> @param cell_type_l      Cell type of left element
    !> @param cell_type_r      Cell type of right element
    !> @param ind_l_in         Solution index of left element
    !> @param ind_r_in         Solution index of right element
    !> @param dof              Number of real degrees of freedom
    !> @param geom_l           Left AMR cell geometry cell coordinates
    !> @param geom_r           Right AMR cell geometry cell coordinates
    !> @param ql               Solution pointer of left element
    !> @param qr               Solution pointer of right element
    !> @param res              Residual pointer (all elements)
    !>
    subroutine solver_residual_face_full_unstructured_3d(nfields,nelem_subgrid,side_l,side_r,&
                                                         cell_type_l,cell_type_r,ind_l_in,ind_r_in,dof,&
                                                         geom_l,geom_r,ql,qr,res)
        use my_kinddefs
        use flux_3d_module
        use normals_3d_module
        use integrate_3d_module
        use subgrid_index_3d_module
        implicit none

        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: side_l
        integer(i4),intent(in)      :: side_r
        integer(i4),intent(in)      :: cell_type_l
        integer(i4),intent(in)      :: cell_type_r
        integer(i4),intent(in)      :: ind_l_in
        integer(i4),intent(in)      :: ind_r_in
        integer(i4),intent(in)      :: dof
        real(dp),   intent(in)      :: geom_l(3,8)
        real(dp),   intent(in)      :: geom_r(3,8)
        real(dp),   intent(in)      :: ql(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(in)      :: qr(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(inout)   :: res(dof)


        real(dp)    :: normal_vec(3,nelem_subgrid,nelem_subgrid)
        real(dp)    :: flux(5)
        integer(i4) :: i_l,j_l,k_l
        integer(i4) :: i_r,j_r,k_r
        integer(i4) :: n,m
        
        integer(i4) :: ind_l
        integer(i4) :: ind_r

        !incoming index is c-based
        ind_l = ind_l_in + 1    
        ind_r = ind_r_in + 1
        
        
        call normals_3d(nelem_subgrid,side_l,geom_l,normal_vec)
        
        
        do m = 1,nelem_subgrid
        do n = 1,nelem_subgrid
            
            call subgrid_index_3d(n,m,nelem_subgrid,side_l,i_l,j_l,k_l)
            call subgrid_index_3d(n,m,nelem_subgrid,side_r,i_r,j_r,k_r)
            
            call flux_roe_3d(ql(:,i_l,j_l,k_l),qr(:,i_r,j_r,k_r),normal_vec(:,n,m),flux)
            if(ind_l < dof) call integrate_3d_pos(i_l,j_l,k_l,nfields,nelem_subgrid,flux,res(ind_l))
            if(ind_r < dof) call integrate_3d_neg(i_r,j_r,k_r,nfields,nelem_subgrid,flux,res(ind_r))
        
        end do
        end do
        

    end subroutine


end module
