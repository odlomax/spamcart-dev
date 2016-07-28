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

! string manipulation module
module m_string

   use m_kind_parameters

   implicit none
   
   contains
   
   ! split space delineated string up into reals
   pure function string_to_real(string) result (value)
   
      ! argument declarations
      character(kind=chr_kind,len=string_length),intent(in) :: string     ! character string
      
      ! result declaration
      real(kind=rel_kind),allocatable :: value(:)                          ! array of reals
      
      ! variable declarations
      integer(kind=int_kind) :: j                                          ! counter
      integer(kind=int_kind) :: i_sub_string                               ! position of sub_string in string
      integer(kind=int_kind) :: n_value                                    ! number of entries in string
      character(kind=chr_kind,len=string_length) :: trimmed_string         ! trimmed version of string variable
      
      ! count number of sub_string delineated entries
      trimmed_string=adjustl(string)
      n_value=0
      do
      
         if (len(trim(adjustl(trimmed_string)))==0) exit
      
         i_sub_string=index(trimmed_string," ")
         n_value=n_value+1
         trimmed_string=adjustl(trimmed_string(i_sub_string+1:))
      
      end do
      
      ! read string and assign to value array
      allocate(value(n_value))
      trimmed_string=adjustl(string)
      do j=1,n_value
      
         i_sub_string=index(trimmed_string," ")
         read(trimmed_string(:i_sub_string-1),*) value(j)
         trimmed_string=adjustl(trimmed_string(i_sub_string+1:))
      
      end do
      
      return
      
   end function

end module