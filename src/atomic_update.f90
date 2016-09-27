! The MIT License (MIT)
! 
! Copyright (c) 2016 odlomax
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

! atomic update subroutines
module m_atomic_update

   use m_kind_parameters
   implicit none
   
   contains
   
   ! add real b to real a
   subroutine atomic_real_add(a,b)
   
      ! argument declarations
      real(kind=rel_kind),intent(inout) :: a
      real(kind=rel_kind),intent(in) :: b
      
      !$omp atomic update
      a=a+b
      
      return
   
   end subroutine
   
   ! add real b to real a
   subroutine atomic_integer_add(a,b)
   
      ! argument declarations
      integer(kind=int_kind),intent(inout) :: a
      integer(kind=int_kind),intent(in) :: b
      
      !$omp atomic update
      a=a+b
      
      return
   
   end subroutine
   
end module