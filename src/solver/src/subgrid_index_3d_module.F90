!     
! \file   subgrid_index_3d_module.F90
! \author akirby
!
! Created on August 31, 2018, 10:23 AM
!


module subgrid_index_3d_module
    implicit none
    contains


    subroutine subgrid_index_3d(ind1,ind2,end_ind,side,i,j,k)
        use my_kinddefs
        implicit none

        integer(i4),intent(in) :: ind1
        integer(i4),intent(in) :: ind2
        integer(i4),intent(in) :: end_ind
        integer(i4),intent(in) :: side
        integer(i4),intent(out):: i
        integer(i4),intent(out):: j
        integer(i4),intent(out):: k
        
        
        select case(side)
            
            case(0)
                i = 1
                j = ind1
                k = ind2
                
            case(1)
                i = end_ind
                j = ind1
                k = ind2
                
            case(2)
                i = ind1
                j = 1
                k = ind2
                
            case(3)
                i = ind1
                j = end_ind
                k = ind2
                
            case(4)
                i = ind1
                j = ind2
                k = 1
                
            case(5)
                i = ind1
                j = ind2
                k = end_ind
                
        end select


    end subroutine


end module
