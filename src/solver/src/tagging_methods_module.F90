!     
! \file   tagging_methods_module.F90
! \author akirby
!
! Created on September 10, 2018, 10:01 AM
!


module tagging_methods_module
    implicit none
    contains
    
    
    subroutine tag_density_below(tolerance,dim,unstructured,nfields,nsub_elem,cell_type,elem_h,elem_xyz,elem_soln,tag)
        use my_kinddefs
        implicit none
        
        real(dp),   intent(in)      :: tolerance
        integer(i4),intent(in)      :: dim
        integer(i4),intent(in)      :: unstructured
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nsub_elem
        integer(i4),intent(in)      :: cell_type
        real(dp),   intent(in)      :: elem_h
        real(dp),   intent(in)      :: elem_xyz(3,2**dim)
        real(dp),   intent(in)      :: elem_soln(nfields,nsub_elem**dim)
        integer(i4),intent(inout)   :: tag
        
        
        integer(i4) :: i
        
        
        do i = 1,nsub_elem**dim
            if(elem_soln(1,i) <= tolerance) then
                tag = 1
                return
            end if
        end do
        
        
    end subroutine
    
    
    subroutine tag_density_above(tolerance,dim,unstructured,nfields,nsub_elem,cell_type,elem_h,elem_xyz,elem_soln,tag)
        use my_kinddefs
        implicit none
        
        real(dp),   intent(in)      :: tolerance
        integer(i4),intent(in)      :: dim
        integer(i4),intent(in)      :: unstructured
        integer(i4),intent(in)      :: nfields
        integer(i4),intent(in)      :: nsub_elem
        integer(i4),intent(in)      :: cell_type
        real(dp),   intent(in)      :: elem_h
        real(dp),   intent(in)      :: elem_xyz(3,2**dim)
        real(dp),   intent(in)      :: elem_soln(nfields,nsub_elem**dim)
        integer(i4),intent(inout)   :: tag
        
        
        integer(i4) :: i
        
        
        do i = 1,nsub_elem**dim
            if(elem_soln(1,i) >= tolerance)then
                tag = 1
                return
            end if
        end do
        
        
    end subroutine
    
    
end module
