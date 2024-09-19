!>     
!> \file   solver_residual_face_hang_full_structured_2d.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D hang-full face residual structured grid
!>
!> Created on September 5, 2018, 2:27 PM
!>


!> 2D hang-full face residual structured grid
!>
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR element
!> @param side              Side id (face) of AMR left element 
!> @param level_l           Level of AMR left element
!> @param level_r           Level of AMR right element
!> @param type_l1           Cell type of AMR left element #1
!> @param type_l2           Cell type of AMR left element #2
!> @param type_r            Cell type of AMR right element
!> @param ind_l1_in         Solution index of AMR left element #1 (c-based)
!> @param ind_l2_in         Solution index of AMR left element #2 (c-based)
!> @param ind_r_in          Solution index of AMR right element (c-based)
!> @param dof               Number of real degrees of freedom
!> @param h_base            Coarse level grid element length
!> @param xlo               Lower corner x-coordinate of AMR element face
!> @param ylo               Lower corner y-coordinate of AMR element face
!> @param zlo               Lower corner z-coordinate of AMR element face
!> @param ql_1              Solution pointer of left element #1
!> @param ql_2              Solution pointer of left element #2
!> @param qr                Solution pointer of right element
!> @param res               Residual pointer (all elements)
!>
subroutine solver_residual_face_hang_full_structured_2d(nfields,nelem_subgrid,side,level_l,level_r,&
                type_l1,type_l2,type_r,ind_l1_in,ind_l2_in,ind_r_in,dof,h_base,xlo,ylo,zlo,ql_1,ql_2,qr,res) bind(C)
    use iso_c_binding
    use my_kinddefs
    use flux_2d_module
    use normals_2d_module
    use integrate_2d_module
    implicit none
    
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side
    integer(c_int),intent(in)      :: level_l
    integer(c_int),intent(in)      :: level_r
    integer(c_int),intent(in)      :: type_l1
    integer(c_int),intent(in)      :: type_l2
    integer(c_int),intent(in)      :: type_r
    integer(c_int),intent(in)      :: ind_l1_in
    integer(c_int),intent(in)      :: ind_l2_in
    integer(c_int),intent(in)      :: ind_r_in
    integer(c_int),intent(in)      :: dof
    real(c_double),intent(in)      :: h_base
    real(c_double),intent(in)      :: xlo(2)
    real(c_double),intent(in)      :: ylo(2)
    real(c_double),intent(in)      :: zlo(2)
    real(c_double),intent(in)      :: ql_1(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: ql_2(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(inout)   :: res(dof)
    
    
    real(dp)    :: normal_vec(2)
    real(dp)    :: flux(4)
    real(dp)    :: cell_h
    integer(i4) :: i,j

    integer(i4) :: ind_l1
    integer(i4) :: ind_l2
    integer(i4) :: ind_r

    !incoming index is c-based
    ind_l1 = ind_l1_in + 1
    ind_l2 = ind_l2_in + 1
    ind_r = ind_r_in + 1

    
    ! use hanging side level to get the area correct
    cell_h = h_base * (1.0_dp / (2.0_dp**(level_l))) / real(nelem_subgrid,dp)
    
    
    select case(side)

        case(1) !xhi

            normal_vec(1) = cell_h
            normal_vec(2) = 0.0_dp

            do j = 1,nelem_subgrid

                !element 1
                call flux_roe_2d(ql_1(:,nelem_subgrid,j),qr(:,1,j),normal_vec,flux)
                if(ind_l1 < dof) call integrate_2d_pos(nelem_subgrid,j,nfields,nelem_subgrid,flux,res(ind_l1))
                if(ind_r  < dof) call integrate_2d_neg(1,j,            nfields,nelem_subgrid,flux,res(ind_r))
                
                !element 2
                call flux_roe_2d(ql_2(:,nelem_subgrid,j),qr(:,1,j),normal_vec,flux)
                if(ind_l2 < dof) call integrate_2d_pos(nelem_subgrid,j,nfields,nelem_subgrid,flux,res(ind_l2))
                if(ind_r  < dof) call integrate_2d_neg(1,j,            nfields,nelem_subgrid,flux,res(ind_r))
                
            end do

        case(3) !yhi

            normal_vec(1) = 0.0_dp
            normal_vec(2) = cell_h

            do i = 1,nelem_subgrid

                !element 1
                call flux_roe_2d(ql_1(:,i,nelem_subgrid),qr(:,i,1),normal_vec,flux)
                if(ind_l1 < dof) call integrate_2d_pos(i,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l1))
                if(ind_r  < dof) call integrate_2d_neg(i,1,            nfields,nelem_subgrid,flux,res(ind_r))
                
                !element 2
                call flux_roe_2d(ql_2(:,i,nelem_subgrid),qr(:,i,1),normal_vec,flux)
                if(ind_l2 < dof) call integrate_2d_pos(i,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l2))
                if(ind_r  < dof) call integrate_2d_neg(i,1,            nfields,nelem_subgrid,flux,res(ind_r))

            end do

            
    end select
        
    
end subroutine
