!>     
!> \file   solver_residual_face_hang_unstructured_2d.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D face hang residual unstructured grid
!>
!> Created on September 5, 2018, 4:08 PM
!>

!> 2D face hang residual unstructured grid
!>
!> @param nfields           Number of fields (unknowns) per grid point
!> @param nelem_subgrid     Number of 1D subgrid elements in AMR element
!> @param side_l            Side of AMR left element
!> @param side_r1           Side of AMR right element #1
!> @param side_r2           Side of AMR right element #2
!> @param type_l            Cell type of AMR left element
!> @param type_r1           Cell type of AMR right element #1
!> @param type_r2           Cell type of AMR right element #2
!> @param ind_l_in          Solution index of AMR left element
!> @param ind_r1_in         Solution index of AMR right element #1
!> @param ind_r2_in         Solution index of AMR right element #2
!> @param dof               Total number of real degrees of freedom
!> @param geom_l            Geometry coordinates of AMR left element
!> @param geom_r1           Geometry coordinates of AMR right element #1
!> @param geom_r2           Geometry coordinates of AMR right element #2
!> @param ql                Solution pointer of AMR left element
!> @param qr_1              Solution pointer of AMR right element #1
!> @param qr_2              Solution pointer of AMR right element #2
!> @param res               Residual pointer
!>
subroutine solver_residual_face_hang_unstructured_2d(nfields,nelem_subgrid,side_l,side_r1,side_r2,type_l,type_r1,type_r2,&
            ind_l_in,ind_r1_in,ind_r2_in,dof,geom_l,geom_r1,geom_r2,ql,qr_1,qr_2,res) bind(C)
    use iso_c_binding
    use my_kinddefs
    use flux_2d_module
    use normals_2d_module
    use integrate_2d_module
    use subgrid_index_2d_module
    implicit none
    
    
    integer(c_int),intent(in)      :: nfields
    integer(c_int),intent(in)      :: nelem_subgrid
    integer(c_int),intent(in)      :: side_l
    integer(c_int),intent(in)      :: side_r1
    integer(c_int),intent(in)      :: side_r2
    integer(c_int),intent(in)      :: type_l
    integer(c_int),intent(in)      :: type_r1
    integer(c_int),intent(in)      :: type_r2
    integer(c_int),intent(in)      :: ind_l_in
    integer(c_int),intent(in)      :: ind_r1_in
    integer(c_int),intent(in)      :: ind_r2_in
    integer(c_int),intent(in)      :: dof
    real(c_double),intent(in)      :: geom_l(3,4)
    real(c_double),intent(in)      :: geom_r1(3,4)
    real(c_double),intent(in)      :: geom_r2(3,4)
    real(c_double),intent(in)      :: ql(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_1(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(in)      :: qr_2(nfields,nelem_subgrid,nelem_subgrid)
    real(c_double),intent(inout)   :: res(dof)
    
    
    
    real(dp)    :: normal_vec1(2,nelem_subgrid)
    real(dp)    :: normal_vec2(2,nelem_subgrid)
    real(dp)    :: flux(4)
    integer(i4) :: i_l,j_l
    integer(i4) :: i_r,j_r
    integer(i4) :: n
    
    integer(i4) :: ind_l
    integer(i4) :: ind_r1
    integer(i4) :: ind_r2

    !incoming index is c-based
    ind_l = ind_l_in + 1
    ind_r1 = ind_r1_in + 1
    ind_r2 = ind_r2_in + 1
    
    
    call normals_2d(nelem_subgrid,side_r1,geom_r1,normal_vec1)
    call normals_2d(nelem_subgrid,side_r2,geom_r2,normal_vec2)
    
    
    do n = 1,nelem_subgrid
    
        call subgrid_index_2d(n,nelem_subgrid,side_l,i_l,j_l)
        
        !element 1
        call subgrid_index_2d(n,nelem_subgrid,side_r1,i_r,j_r)
        call flux_roe_2d(qr_1(:,i_r,j_r),ql(:,i_l,j_l),normal_vec1(:,n),flux)
        if(ind_r1 < dof) call integrate_2d_pos(i_r,j_r,nfields,nelem_subgrid,flux,res(ind_r1))
        if(ind_l  < dof) call integrate_2d_neg(i_l,j_l,nfields,nelem_subgrid,flux,res(ind_l))
        
        !element 2
        call subgrid_index_2d(n,nelem_subgrid,side_r2,i_r,j_r)
        call flux_roe_2d(qr_2(:,i_r,j_r),ql(:,i_l,j_l),normal_vec2(:,n),flux)
        if(ind_r2 < dof) call integrate_2d_pos(i_r,j_r,nfields,nelem_subgrid,flux,res(ind_r2))
        if(ind_l  < dof) call integrate_2d_neg(i_l,j_l,nfields,nelem_subgrid,flux,res(ind_l))
        
    end do

        
    
end subroutine
