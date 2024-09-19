
!     
! \file   initial_boundary_2d_module.F90
! \author akirby
!
! Created on August 31, 2018, 11:28 AM
!


module initial_boundary_2d_module
    implicit none
    contains


    subroutine initial_boundary_2d(unstructured,cell_geom,cell_type)
        use my_kinddefs
        implicit none

        integer(i4),intent(in)      :: unstructured
        real(dp),   intent(in)      :: cell_geom(3,4) !always 3 dimensional
        integer(i4),intent(inout)   :: cell_type
        
        
        cell_type = 0
        if(unstructured==1)then
            
            !assign slip boundary
            if(abs(cell_geom(1,1)) < 10.0_dp)then
                cell_type = 1
            end if
            
            
        end if
        
        
    end subroutine
    
    
end module
        
