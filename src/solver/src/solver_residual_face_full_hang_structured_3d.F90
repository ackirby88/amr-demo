!>     
!> \file   solver_residual_face_full_hang_structured_3d.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 3D full-hang face residual structured grid
!>
!> Created on September 5, 2018, 3:14 PM
!>


!> 3D full-hang face residual structured grid
!>
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
!> @param side              Side id (face) of AMR left quad 
!> @param level_l           Level of AMR left element
!> @param level_r           Level of AMR right element
!> @param type_l            Cell type of AMR left element
!> @param type_r1           Cell type of AMR right element #1
!> @param type_r2           Cell type of AMR right element #2
!> @param type_r3           Cell type of AMR right element #3
!> @param type_r4           Cell type of AMR right element #4
!> @param ind_l_in          Solution index of AMR left element
!> @param ind_r1_in         Solution index of AMR right element #1
!> @param ind_r2_in         Solution index of AMR right element #2
!> @param ind_r3_in         Solution index of AMR right element #3
!> @param ind_r4_in         Solution index of AMR right element #4
!> @param dof               Number of real degrees of freedom
!> @param h_base            Coarse level grid element length
!> @param xlo               Lower corner x-coordinate of AMR element face
!> @param ylo               Lower corner y-coordinate of AMR element face
!> @param zlo               Lower corner z-coordinate of AMR element face
!> @param ql                Solution pointer of AMR left element
!> @param qr_1              Solution pointer of AMR right element #1
!> @param qr_2              Solution pointer of AMR right element #2
!> @param qr_3              Solution pointer of AMR right element #3
!> @param qr_4              Solution pointer of AMR right element #4
!> @param res               Residual pointer (all elements)
!>
subroutine solver_residual_face_full_hang_structured_3d(nfields,nelem_subgrid,side,level_l,level_r,&
                type_l,type_r1,type_r2,type_r3,type_r4,ind_l_in,ind_r1_in,ind_r2_in,ind_r3_in,ind_r4_in,&
                dof,h_base,xlo,ylo,zlo,ql,qr_1,qr_2,qr_3,qr_4,res) bind(C)
    use iso_c_binding
    use my_kinddefs
    use flux_3d_module
    use normals_3d_module
    use integrate_3d_module
    implicit none
    
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side
    integer(c_int),intent(in)      :: level_l
    integer(c_int),intent(in)      :: level_r
    integer(c_int),intent(in)      :: type_l
    integer(c_int),intent(in)      :: type_r1
    integer(c_int),intent(in)      :: type_r2
    integer(c_int),intent(in)      :: type_r3
    integer(c_int),intent(in)      :: type_r4
    integer(c_int),intent(in)      :: ind_l_in
    integer(c_int),intent(in)      :: ind_r1_in
    integer(c_int),intent(in)      :: ind_r2_in
    integer(c_int),intent(in)      :: ind_r3_in
    integer(c_int),intent(in)      :: ind_r4_in
    integer(c_int),intent(in)      :: dof
    real(c_double),intent(in)      :: h_base
    real(c_double),intent(in)      :: xlo(2)
    real(c_double),intent(in)      :: ylo(2)
    real(c_double),intent(in)      :: zlo(2)
    real(c_double),intent(in)      :: ql(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_1(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_2(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_3(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_4(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(inout)   :: res(dof)
    
    
    real(dp)    :: normal_vec(3)
    real(dp)    :: flux(5)
    real(dp)    :: cell_h
    real(dp)    :: face_area
    integer(i4) :: i,j,k
    
    integer(i4) :: ind_l
    integer(i4) :: ind_r1
    integer(i4) :: ind_r2
    integer(i4) :: ind_r3
    integer(i4) :: ind_r4
    
    !incoming index is c-based
    ind_l = ind_l_in + 1    
    ind_r1 = ind_r1_in + 1
    ind_r2 = ind_r2_in + 1
    ind_r3 = ind_r3_in + 1
    ind_r4 = ind_r4_in + 1


    ! use hanging side level to get the area correct
    cell_h = h_base * (1.0_dp / (2.0_dp**(level_r))) / real(nelem_subgrid,dp)
    face_area = cell_h*cell_h
    
    
    select case(side)

        case(1) !xhi

                normal_vec(1) = face_area
                normal_vec(2) = 0.0_dp
                normal_vec(3) = 0.0_dp

                do k = 1,nelem_subgrid
                do j = 1,nelem_subgrid

                    !element 1
                    call flux_roe_3d(ql(:,nelem_subgrid,j,k),qr_1(:,1,j,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(nelem_subgrid,j,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r1 < dof) call integrate_3d_neg(1,j,k,            nfields,nelem_subgrid,flux,res(ind_r1))
                    
                    !element 2
                    call flux_roe_3d(ql(:,nelem_subgrid,j,k),qr_2(:,1,j,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(nelem_subgrid,j,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r2 < dof) call integrate_3d_neg(1,j,k,            nfields,nelem_subgrid,flux,res(ind_r2))
                    
                    !element 3
                    call flux_roe_3d(ql(:,nelem_subgrid,j,k),qr_3(:,1,j,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(nelem_subgrid,j,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r3 < dof) call integrate_3d_neg(1,j,k,            nfields,nelem_subgrid,flux,res(ind_r3))
                    
                    !element 4
                    call flux_roe_3d(ql(:,nelem_subgrid,j,k),qr_4(:,1,j,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(nelem_subgrid,j,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r4 < dof) call integrate_3d_neg(1,j,k,            nfields,nelem_subgrid,flux,res(ind_r4))

                end do
                end do

            case(3) !yhi

                normal_vec(1) = 0.0_dp
                normal_vec(2) = face_area
                normal_vec(3) = 0.0_dp

                do k = 1,nelem_subgrid
                do i = 1,nelem_subgrid

                    !element 1
                    call flux_roe_3d(ql(:,i,nelem_subgrid,k),qr_1(:,i,1,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,nelem_subgrid,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r1 < dof) call integrate_3d_neg(i,1,k,            nfields,nelem_subgrid,flux,res(ind_r1))

                    !element 2
                    call flux_roe_3d(ql(:,i,nelem_subgrid,k),qr_2(:,i,1,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,nelem_subgrid,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r2 < dof) call integrate_3d_neg(i,1,k,            nfields,nelem_subgrid,flux,res(ind_r2))
                    
                    !element 3
                    call flux_roe_3d(ql(:,i,nelem_subgrid,k),qr_3(:,i,1,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,nelem_subgrid,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r3 < dof) call integrate_3d_neg(i,1,k,            nfields,nelem_subgrid,flux,res(ind_r3))
                    
                    !element 4
                    call flux_roe_3d(ql(:,i,nelem_subgrid,k),qr_4(:,i,1,k),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,nelem_subgrid,k,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r4 < dof) call integrate_3d_neg(i,1,k,            nfields,nelem_subgrid,flux,res(ind_r4))
                    
                end do
                end do

            case(5) !zhi

                normal_vec(1) = 0.0_dp
                normal_vec(2) = 0.0_dp
                normal_vec(3) = face_area

                do j = 1,nelem_subgrid
                do i = 1,nelem_subgrid

                    !element 1
                    call flux_roe_3d(ql(:,i,j,nelem_subgrid),qr_1(:,i,j,1),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,j,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r1 < dof) call integrate_3d_neg(i,j,1,            nfields,nelem_subgrid,flux,res(ind_r1))

                    !element 2
                    call flux_roe_3d(ql(:,i,j,nelem_subgrid),qr_2(:,i,j,1),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,j,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r2 < dof) call integrate_3d_neg(i,j,1,            nfields,nelem_subgrid,flux,res(ind_r2))
                    
                    !element 3
                    call flux_roe_3d(ql(:,i,j,nelem_subgrid),qr_3(:,i,j,1),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,j,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r3 < dof) call integrate_3d_neg(i,j,1,            nfields,nelem_subgrid,flux,res(ind_r3))
                    
                    !element 4
                    call flux_roe_3d(ql(:,i,j,nelem_subgrid),qr_4(:,i,j,1),normal_vec,flux)
                    if(ind_l  < dof) call integrate_3d_pos(i,j,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r4 < dof) call integrate_3d_neg(i,j,1,            nfields,nelem_subgrid,flux,res(ind_r4))
                    
                end do
                end do
                
    end select
    
    
end subroutine
