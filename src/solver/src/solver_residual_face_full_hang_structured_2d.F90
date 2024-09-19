!>     
!> \file   solver_residual_face_full_hang_structured_2d.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D full-hang face residual structured grid
!>
!> Created on September 5, 2018, 2:09 PM
!>


!> 2D full-hang face residual structured grid
!>
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
!> @param side              Side id (face) of AMR left quad 
!> @param level_l           Level of AMR left element
!> @param level_r           Level of AMR right element
!> @param type_l            Cell type of AMR left element
!> @param type_r1           Cell type of AMR right element #1
!> @param type_r2           Cell type of AMR right element #2
!> @param ind_l_in          Solution index of AMR left element
!> @param ind_r1_in         Solution index of AMR right element #1
!> @param ind_r2_in         Solution index of AMR right element #2
!> @param dof               Number of real degrees of freedom
!> @param h_base            Coarse level grid element length
!> @param xlo               lower corner x-coordinate of AMR quad face
!> @param ylo               lower corner y-coordinate of AMR quad face
!> @param zlo               lower corner z-coordinate of AMR quad face
!> @param ql                Solution pointer of left element
!> @param qr_1              Solution pointer of right element #1
!> @param qr_2              Solution pointer of right element #2
!> @param res               Residual pointer (all elements)
!>
subroutine solver_residual_face_full_hang_structured_2d(nfields,nelem_subgrid,side,level_l,level_r,&
                type_l,type_r1,type_r2,ind_l_in,ind_r1_in,ind_r2_in,dof,h_base,xlo,ylo,zlo,ql,qr_1,qr_2,res) bind(C)
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
    integer(c_int),intent(in)      :: type_l
    integer(c_int),intent(in)      :: type_r1
    integer(c_int),intent(in)      :: type_r2
    integer(c_int),intent(in)      :: ind_l_in
    integer(c_int),intent(in)      :: ind_r1_in
    integer(c_int),intent(in)      :: ind_r2_in
    integer(c_int),intent(in)      :: dof
    real(c_double),intent(in)      :: h_base
    real(c_double),intent(in)      :: xlo(2)
    real(c_double),intent(in)      :: ylo(2)
    real(c_double),intent(in)      :: zlo(2)
    real(c_double),intent(in)      :: ql(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_1(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_2(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(inout)   :: res(dof)
    
    
    real(dp)    :: normal_vec(2)
    real(dp)    :: flux(4)
    real(dp)    :: cell_h
    integer(i4) :: i,j
    
    integer(i4) :: ind_l
    integer(i4) :: ind_r1
    integer(i4) :: ind_r2
    
    !incoming index is c-based
    ind_l = ind_l_in + 1    
    ind_r1 = ind_r1_in + 1
    ind_r2 = ind_r2_in + 1
    

    ! use hanging side level to get the area correct
    cell_h = h_base * (1.0_dp / (2.0_dp**(level_r))) / real(nelem_subgrid,dp)
    
    
    select case(side)

        case(1) !xhi

            normal_vec(1) = cell_h
            normal_vec(2) = 0.0_dp

            do j = 1,nelem_subgrid

                !element 1
                call flux_roe_2d(ql(:,nelem_subgrid,j),qr_1(:,1,j),normal_vec,flux)
                if(ind_l  < dof) call integrate_2d_pos(nelem_subgrid,j,nfields,nelem_subgrid,flux,res(ind_l))
                if(ind_r1 < dof) call integrate_2d_neg(1,j,            nfields,nelem_subgrid,flux,res(ind_r1))
                
                !element 2
                call flux_roe_2d(ql(:,nelem_subgrid,j),qr_2(:,1,j),normal_vec,flux)
                if(ind_l  < dof) call integrate_2d_pos(nelem_subgrid,j,nfields,nelem_subgrid,flux,res(ind_l))
                if(ind_r2 < dof) call integrate_2d_neg(1,j,            nfields,nelem_subgrid,flux,res(ind_r2))
                
            end do

        case(3) !yhi

            normal_vec(1) = 0.0_dp
            normal_vec(2) = cell_h

            do i = 1,nelem_subgrid

                !element 1
                call flux_roe_2d(ql(:,i,nelem_subgrid),qr_1(:,i,1),normal_vec,flux)
                if(ind_l  < dof) call integrate_2d_pos(i,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l))
                if(ind_r1 < dof) call integrate_2d_neg(i,1,            nfields,nelem_subgrid,flux,res(ind_r1))
                
                !element 2
                call flux_roe_2d(ql(:,i,nelem_subgrid),qr_2(:,i,1),normal_vec,flux)
                if(ind_l  < dof) call integrate_2d_pos(i,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l))
                if(ind_r2 < dof) call integrate_2d_neg(i,1,            nfields,nelem_subgrid,flux,res(ind_r2))

            end do
            

    end select
    
    
end subroutine
