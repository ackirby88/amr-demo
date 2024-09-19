!     
! \file   subgrid_index_2d_module.F90
! \author akirby
!
! Created on August 31, 2018, 10:18 AM
!


module subgrid_index_2d_module
    implicit none
    contains


    subroutine subgrid_index_2d(ind,end_ind,side,i,j)
        use my_kinddefs
        implicit none

        integer(i4),intent(in) :: ind
        integer(i4),intent(in) :: end_ind
        integer(i4),intent(in) :: side
        integer(i4),intent(out):: i
        integer(i4),intent(out):: j
        
        
        select case(side)
            
            case(0)
                i = 1
                j = ind
                
            case(1)
                i = end_ind
                j = ind
                
            case(2)
                i = ind
                j = 1
                
            case(3)
                i = ind
                j = end_ind
                
        end select
                

    end subroutine


end module
