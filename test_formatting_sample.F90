!>
! Sample poorly-formatted Fortran file for testing
! Contains various formatting violations to test the formatter
!

module test_module
  implicit none
  
  ! BAD: Mixed indentation
  integer, parameter :: kind_phys = 8
   real(kind_phys), parameter :: cp = 1004.6_kind_phys
     real(kind_phys), parameter :: rd = 287.0_kind_phys
real(kind_phys), parameter :: rv = 461.6_kind_phys

  ! BAD: Multiple variables on one line with inconsistent spacing
  real(kind_phys) :: a, b,c,d , e
  integer :: i,j , k
  real(kind_phys) :: very_long_variable_name_that_exceeds_the_line_limit_when_combined_with_other_declarations
  
contains

  subroutine test_sub(x, y, z, &
      a, b, c, &
       d, e, f)
    ! BAD: Inconsistent indentation and parameter alignment
    real(kind_phys), intent(in) :: x,y,z
      real(kind_phys), intent(in) :: a, b, c
    real(kind_phys),intent(inout)::d,e,f
    
    ! BAD: Long lines
    real(kind_phys) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9, temp10
    
    ! BAD: Poor operator spacing
    x=y+z
    a =   b   -   c
    d= e* f
    
    ! BAD: Mixed continuation alignment
    temp1 = x * y + &
         z * a + &
       b * c + &
    d * e
    
    temp2 = function_call(x, y, z, &
      a, b, c, &
       d, e, f)
    
  end subroutine test_sub
  
  function function_call(p1,p2,p3,p4,p5,p6) result(res)
    ! BAD: Inconsistent spacing in declarations
    real(kind_phys),intent(in)::p1,p2,p3
    real(kind_phys) , intent(in) :: p4 , p5 , p6
    real(kind_phys)::res
    
    res = p1 + p2 + p3 + p4 + p5 + p6
    
  end function function_call
  
end module test_module
