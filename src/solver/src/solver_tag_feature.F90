!     
! File:   solver_tag_feature.F90
! Author: akirby
!
! Created on September 10, 2018, 9:17 AM
!


subroutine solver_tag_feature(ntags,tag_methods,tag_tolerances, &
                             dim,unstructured,nfields,nsub_elem,&
                             cell_type,elem_h,elem_xyz,elem_soln,tag) bind(C)
    use iso_c_binding
    use tagging_methods_module
    implicit none
    
    integer(c_int),intent(in) :: ntags
    integer(c_int),intent(in) :: tag_methods(ntags)
    real(c_double),intent(in) :: tag_tolerances(ntags)
    integer(c_int),intent(in) :: dim
    integer(c_int),intent(in) :: unstructured
    integer(c_int),intent(in) :: nfields
    integer(c_int),intent(in) :: nsub_elem
    integer(c_int),intent(in) :: cell_type
    real(c_double),intent(in) :: elem_h
    real(c_double),intent(in) :: elem_xyz(3,2**dim)
    real(c_double),intent(in) :: elem_soln(nfields,nsub_elem**dim)
    integer(c_int),intent(out):: tag
    
    
    integer(c_int) :: n
    
    
    tag = 0
    do n = 1,ntags
        
        select case(tag_methods(n))
            case(1) !tag density below
                call tag_density_below(tag_tolerances(n),dim,unstructured,nfields,nsub_elem,cell_type,elem_h,elem_xyz,elem_soln,tag)
            case(2) !tag density above
                call tag_density_above(tag_tolerances(n),dim,unstructured,nfields,nsub_elem,cell_type,elem_h,elem_xyz,elem_soln,tag)
            case(3) !tag vortcity
                !call tag_vorticity(dim,unstructured,nfields,nsub_elem,cell_type,elem_h,elem_xyz,elem_soln,tag)
            case(4) !tag qcriterion
                !call tag_qcriterion(dim,unstructured,nfields,nsub_elem,cell_type,elem_h,elem_xyz,elem_soln,tag)                
            case default
                
        end select
        
    end do
    
    
end subroutine
