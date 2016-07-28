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

! triangular array
module m_triangular_array

   use m_kind_parameters
   
   implicit none
   
   private
   public :: triangular_array
   
   
   ! row of reals
   type :: row_element
   
      real(kind=rel_kind),allocatable :: column(:)              ! elements in a row
   
   end type
   
   ! triangular array of reals
   type :: triangular_array
   
      type(row_element),allocatable :: row(:)                   ! rows of array
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: destroy
   
   end type
   
   contains
   
   ! allocate array and set all values to 0
   ! note: n_column must be >= n_row
   pure subroutine initialise(self,n_column,n_row,top_right)
   
      ! argument declarations
      class(triangular_array),intent(inout) :: self             ! triangular array object
      integer(kind=int_kind),intent(in) :: n_column             ! number of columns
      integer(kind=int_kind),intent(in) :: n_row                ! number of rows
      logical(kind=log_kind),intent(in) :: top_right            ! top right or bottom left triangle?
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      ! allocate rows
      allocate(self%row(n_row))
      
      do i=1,size(self%row)
      
         ! allocate columns
         if (top_right) then
            ! set bounds i:n_row
            allocate(self%row(i)%column(i:n_column))
         else
            ! set bounds 1:n_row-(i-1)
            allocate(self%row(i)%column(n_column-i+1))
         end if
      
         self%row(i)%column=0._rel_kind
      
      end do
      
      return
   
   end subroutine
   
   ! deallocate triangular array
   pure subroutine destroy(self)
   
      ! argument declarations
      class(triangular_array),intent(inout) :: self             ! triangular array object
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      ! deallocate rows
      do i=1,size(self%row)
      
         deallocate(self%row(i)%column)
      
      end do
      
      deallocate(self%row)
      
      return
      
   end subroutine

end module