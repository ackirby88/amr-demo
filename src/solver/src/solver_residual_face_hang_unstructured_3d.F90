!>     
!> \file   solver_residual_face_hang_unstructured_3d.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 3D face hang residual unstructured grid
!>
!> Created on September 10, 2018, 9:01 AM
!>

!> 3D face hang residual unstructured grid
!>
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR element
!> @param side_l            Side of AMR left element
!> @param side_r1           Side of AMR right element #1
!> @param side_r2           Side of AMR right element #2
!> @param side_r3           Side of AMR right element #3
!> @param side_r4           Side of AMR right element #4
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
!> @param dof               Total number of real degrees of freedom
!> @param geom_l            Geometry coordinates of AMR left element
!> @param geom_r1           Geometry coordinates of AMR right element #1
!> @param geom_r2           Geometry coordinates of AMR right element #2
!> @param geom_r3           Geometry coordinates of AMR right element #3
!> @param geom_r4           Geometry coordinates of AMR right element #4
!> @param ql                Solution pointer of AMR left element
!> @param qr_1              Solution pointer of AMR right element #1
!> @param qr_2              Solution pointer of AMR right element #2
!> @param qr_3              Solution pointer of AMR right element #3
!> @param qr_4              Solution pointer of AMR right element #4
!> @param res               Residual pointer
!>
subroutine solver_residual_face_hang_unstructured_3d(nfields,nelem_subgrid,side_l,side_r1,side_r2,side_r3,side_r4,&
                                                     type_l,type_r1,type_r2,type_r3,type_r4,&
                                                     ind_l_in,ind_r1_in,ind_r2_in,ind_r3_in,ind_r4_in,dof, &
                                                     geom_l,geom_r1,geom_r2,geom_r3,geom_r4,&
                                                     ql,qr_1,qr_2,qr_3,qr_4,res) bind(C)
    use iso_c_binding
    use my_kinddefs
    use flux_3d_module
    use normals_3d_module
    use integrate_3d_module
    use subgrid_index_3d_module
    implicit none
    
    
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side_l
    integer(c_int),intent(in)      :: side_r1
    integer(c_int),intent(in)      :: side_r2
    integer(c_int),intent(in)      :: side_r3
    integer(c_int),intent(in)      :: side_r4
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
    real(c_double),intent(in)      :: geom_l(3,8)
    real(c_double),intent(in)      :: geom_r1(3,8)
    real(c_double),intent(in)      :: geom_r2(3,8)
    real(c_double),intent(in)      :: geom_r3(3,8)
    real(c_double),intent(in)      :: geom_r4(3,8)
    real(c_double),intent(in)      :: ql(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_1(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_2(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_3(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_4(nfields,nelem_subgrid,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(inout)   :: res(dof)    
    
    
    real(dp)    :: normal_vec_r1(3,nelem_subgrid,nelem_subgrid)
    real(dp)    :: normal_vec_r2(3,nelem_subgrid,nelem_subgrid)
    real(dp)    :: normal_vec_r3(3,nelem_subgrid,nelem_subgrid)
    real(dp)    :: normal_vec_r4(3,nelem_subgrid,nelem_subgrid)
    real(dp)    :: flux(5)
    integer(i4) :: i_l,j_l,k_l
    integer(i4) :: i_r,j_r,k_r
    integer(i4) :: m,n
    
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
    
    call normals_3d(nelem_subgrid,side_r1,geom_r1,normal_vec_r1)
    call normals_3d(nelem_subgrid,side_r2,geom_r2,normal_vec_r2)
    call normals_3d(nelem_subgrid,side_r3,geom_r3,normal_vec_r3)
    call normals_3d(nelem_subgrid,side_r4,geom_r4,normal_vec_r4)
    
    
    !NOTE: here we are using the hanging elements (qr_*) 
    !      as the left element since the normal is facing 
    !      outward from the hanging element to account for 
    !      the integration surface area
    
    
    do m = 1,nelem_subgrid
    do n = 1,nelem_subgrid
    
        call subgrid_index_3d(n,m,nelem_subgrid,side_l,i_l,j_l,k_l)
        
        !element 1
        call subgrid_index_3d(n,m,nelem_subgrid,side_r1,i_r,j_r,k_r)
        call flux_roe_3d(qr_1(:,i_r,j_r,k_r),ql(:,i_l,j_l,k_l),normal_vec_r1(:,n,m),flux)
        if(ind_r1 < dof) call integrate_3d_pos(i_r,j_r,k_r,nfields,nelem_subgrid,flux,res(ind_r1))
        if(ind_l  < dof) call integrate_3d_neg(i_l,j_l,k_l,nfields,nelem_subgrid,flux,res(ind_l))
        
        !element 2
        call subgrid_index_3d(n,m,nelem_subgrid,side_r2,i_r,j_r,k_r)
        call flux_roe_3d(qr_2(:,i_r,j_r,k_r),ql(:,i_l,j_l,k_l),normal_vec_r2(:,n,m),flux)
        if(ind_r2 < dof) call integrate_3d_pos(i_r,j_r,k_r,nfields,nelem_subgrid,flux,res(ind_r2))
        if(ind_l  < dof) call integrate_3d_neg(i_l,j_l,k_l,nfields,nelem_subgrid,flux,res(ind_l))
        
        !element 3
        call subgrid_index_3d(n,m,nelem_subgrid,side_r3,i_r,j_r,k_r)
        call flux_roe_3d(qr_3(:,i_r,j_r,k_r),ql(:,i_l,j_l,k_l),normal_vec_r3(:,n,m),flux)
        if(ind_r3 < dof) call integrate_3d_pos(i_r,j_r,k_r,nfields,nelem_subgrid,flux,res(ind_r3))
        if(ind_l  < dof) call integrate_3d_neg(i_l,j_l,k_l,nfields,nelem_subgrid,flux,res(ind_l))
        
        !element 4
        call subgrid_index_3d(n,m,nelem_subgrid,side_r4,i_r,j_r,k_r)
        call flux_roe_3d(qr_4(:,i_r,j_r,k_r),ql(:,i_l,j_l,k_l),normal_vec_r4(:,n,m),flux)
        if(ind_r4 < dof) call integrate_3d_pos(i_r,j_r,k_r,nfields,nelem_subgrid,flux,res(ind_r4))
        if(ind_l  < dof) call integrate_3d_neg(i_l,j_l,k_l,nfields,nelem_subgrid,flux,res(ind_l))
        
    end do
    end do
        
    
end subroutine
