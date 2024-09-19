
!>     
!> \file   solver_residual_face_full_structured_2d_module.F90
!> \author akirby
!> \ingroup solver_group
!> \brief 2D full face residual structured grid
!>
!> Created on August 30, 2018, 8:32 AM
!>

module solver_residual_face_full_structured_2d_module
    implicit none
    contains

    
    !> 2D full face residual structured grid
    !>
    !> @param nfields           Number of fields (unknowns) per grid point
    !> @param nelem_subgrid     Number of 1D subgrid elements in AMR quad
    !> @param side              Side id (face) of AMR left quad 
    !> @param level             AMR level
    !> @param cell_type_l       Cell type of left element
    !> @param cell_type_r       Cell type of right element
    !> @param ind_l_in          Solution index of left element (c-based)
    !> @param ind_r_in          Solution index of right element (c-based)
    !> @param dof               Number of real degrees of freedom
    !> @param h_base_level      Coarse level grid element length
    !> @param xlo               lower corner x-coordinate of AMR quad
    !> @param ylo               lower corner y-coordinate of AMR quad
    !> @param zlo               lower corner z-coordinate of AMR quad
    !> @param ql                Solution pointer of left element
    !> @param qr                Solution pointer of right element
    !> @param res               Residual pointer (all elements)
    !>
    subroutine solver_residual_face_full_structured_2d(nfields,nelem_subgrid,side,level,&
                                                       cell_type_l,cell_type_r,ind_l_in,ind_r_in,dof,&
                                                       h_base_level,xlo,ylo,zlo,ql,qr,res)
        use my_kinddefs
        use flux_2d_module
        use normals_2d_module
        use integrate_2d_module
        implicit none

        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nelem_subgrid
        integer(i4),intent(in)      :: side
        integer(i4),intent(in)      :: level
        integer(i4),intent(in)      :: cell_type_l
        integer(i4),intent(in)      :: cell_type_r
        integer(i4),intent(in)      :: ind_l_in
        integer(i4),intent(in)      :: ind_r_in
        integer(i4),intent(in)      :: dof
        real(dp),   intent(in)      :: h_base_level
        real(dp),   intent(in)      :: xlo(2)
        real(dp),   intent(in)      :: ylo(2)
        real(dp),   intent(in)      :: zlo(2)
        real(dp),   intent(in)      :: ql(nfields,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(in)      :: qr(nfields,nelem_subgrid,nelem_subgrid)
        real(dp),   intent(inout)   :: res(dof)


        real(dp)    :: normal_vec(2)
        real(dp)    :: flux(4)
        real(dp)    :: cell_h
        integer(i4) :: i,j
        
        integer(i4) :: ind_l
        integer(i4) :: ind_r

        !incoming index is c-based
        ind_l = ind_l_in + 1    
        ind_r = ind_r_in + 1


        cell_h = h_base_level * (1.0_dp / (2.0_dp**(level))) / real(nelem_subgrid,dp)
        

        select case(side)

            case(1) !xhi

                normal_vec(1) = cell_h
                normal_vec(2) = 0.0_dp

                do j = 1,nelem_subgrid

                    call flux_roe_2d(ql(:,nelem_subgrid,j),qr(:,1,j),normal_vec,flux)
                    if(ind_l < dof) call integrate_2d_pos(nelem_subgrid,j,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r < dof) call integrate_2d_neg(1,j,            nfields,nelem_subgrid,flux,res(ind_r))

                end do

            case(3) !yhi

                normal_vec(1) = 0.0_dp
                normal_vec(2) = cell_h

                do i = 1,nelem_subgrid

                    call flux_roe_2d(ql(:,i,nelem_subgrid),qr(:,i,1),normal_vec,flux)
                    if(ind_l < dof) call integrate_2d_pos(i,nelem_subgrid,nfields,nelem_subgrid,flux,res(ind_l))
                    if(ind_r < dof) call integrate_2d_neg(i,1,            nfields,nelem_subgrid,flux,res(ind_r))

                end do

        end select


    end subroutine


end module
