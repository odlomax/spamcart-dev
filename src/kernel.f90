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

module m_kernel

   use m_kind_parameters
   use m_triangular_array
   
   implicit none
   
   private
   public :: kernel
   
   ! abstract kernel class
   type,abstract :: kernel
      
      real(kind=rel_kind) :: r_support                          ! support radius, in units of h
      real(kind=rel_kind) :: norm                               ! normalisation constant
      real(kind=rel_kind) :: inv_dr                             ! 1/dr
      real(kind=rel_kind),allocatable :: r_array(:)             ! list of radius samples
      real(kind=rel_kind),allocatable :: w_array(:)             ! density lookup table
      real(kind=rel_kind),allocatable :: dw_dr_array(:)         ! density gradient lookup table
      real(kind=rel_kind),allocatable :: d2w_dr2_array(:)       ! density 2nd derivative lookup table
      type(triangular_array) :: sigma_array                     ! column density array
      
      contains
      
      procedure(initialise_virtual),deferred :: initialise
      procedure,non_overridable :: destroy
      procedure,non_overridable :: w
      procedure,non_overridable :: dw_dr
      procedure,non_overridable :: d2w_dr2
      procedure,non_overridable :: sigma
      
   end type
   
   abstract interface
      pure subroutine initialise_virtual(self)
      
         import :: kernel
         
         class(kernel),intent(inout) :: self                   ! kernel object
      
      end subroutine
   end interface

   contains
   
   ! destroy kernel object
   pure subroutine destroy(self)
   
      ! argument declarations
      class(kernel),intent(inout) :: self                       ! kernel object
      
      deallocate(self%r_array)
      deallocate(self%w_array)
      deallocate(self%dw_dr_array)
      call self%sigma_array%destroy()
      
      return
   
   end subroutine
   
   ! kernel density
   elemental function w(self,r) result (w_value)
   
      ! argument declarations
      class(kernel),intent(in) :: self                          ! kernel object
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: w_value                            ! kernel density
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! index
      real(kind=rel_kind) :: r_mod                              ! modified r value
      
      ! make sure r is within range
      r_mod=max(0._rel_kind,min(r,self%r_support))
      
      ! get i
      i=min(1+int(r_mod*self%inv_dr),size(self%r_array)-1)
      
      ! interpolate value from lookup table
      w_value=self%w_array(i)+&
         &(self%w_array(i+1)-self%w_array(i))*&
         &(r_mod-self%r_array(i))*self%inv_dr
      
      return
   
   end function
   
   ! kernel density gradient wrt r
   elemental function dw_dr(self,r) result (dw_dr_value)
   
      ! argument declarations
      class(kernel),intent(in) :: self                          ! kernel object
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: dw_dr_value                        ! kernel density gradient
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! index
      real(kind=rel_kind) :: r_mod                              ! modified r value
      
      ! make sure r is within range
      r_mod=max(0._rel_kind,min(r,self%r_support))
      
      ! get i
      i=min(1+int(r_mod*self%inv_dr),size(self%r_array)-1)
      
      ! interpolate value from lookup table
      dw_dr_value=self%dw_dr_array(i)+&
         &(self%dw_dr_array(i+1)-self%dw_dr_array(i))*&
         &(r_mod-self%r_array(i))*self%inv_dr
      
      return
   
   end function
   
   ! kernel density 2nd derivative wrt r
   elemental function d2w_dr2(self,r) result (d2w_dr2_value)
   
      ! argument declarations
      class(kernel),intent(in) :: self                          ! kernel object
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: d2w_dr2_value                        ! kernel density gradient
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! index
      real(kind=rel_kind) :: r_mod                              ! modified r value
      
      ! make sure r is within range
      r_mod=max(0._rel_kind,min(r,self%r_support))
      
      ! get i
      i=min(1+int(r_mod*self%inv_dr),size(self%r_array)-1)
      
      ! interpolate value from lookup table
      d2w_dr2_value=self%d2w_dr2_array(i)+&
         &(self%d2w_dr2_array(i+1)-self%d2w_dr2_array(i))*&
         &(r_mod-self%r_array(i))*self%inv_dr
      
      return
   
   end function
   
   ! kernel column density
   elemental function sigma(self,b,r) result (sigma_value)
   
      ! argument declarations
      class(kernel),intent(in) :: self                          ! kernel object
      real(kind=rel_kind),intent(in) :: b                       ! impact parameter
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: sigma_value                        ! kernel column density
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                             ! index
      integer(kind=int_kind) :: i_1                             ! index
      integer(kind=int_kind) :: j                               ! index
      real(kind=rel_kind) :: b_mod                              ! modified b value
      real(kind=rel_kind) :: r_mod                              ! modified r value
      real(kind=rel_kind) :: sigma_0                            ! interpolation point
      real(kind=rel_kind) :: sigma_1                            ! interpolation point
      
      ! first interpolation
      
      ! make sure b is within range
      b_mod=max(0._rel_kind,min(b,self%r_support))
      
      ! get j
      j=min(1+int(b_mod*self%inv_dr),size(self%r_array)-1)
      
      ! make sure r is within range
      r_mod=max(min(r,self%r_support),self%r_array(j))
      
      ! get i_0 and i_1
      i_0=max(1+int(r_mod*self%inv_dr),j)
      i_1=min(i_0+1,size(self%r_array))
      
      ! get sigma_0
      sigma_0=self%sigma_array%row(j)%column(i_0)+&
         &(self%sigma_array%row(j)%column(i_1)-self%sigma_array%row(j)%column(i_0))*&
         &(r_mod-self%r_array(i_0))*self%inv_dr
      
      ! second interpolation
      
      ! make sure r is within range
      r_mod=max(min(r,self%r_support),self%r_array(j+1))
      
      ! get i_0 and i_1
      i_0=max(1+int(r_mod*self%inv_dr),j+1)
      i_1=min(i_0+1,size(self%r_array))
      
      ! get sigma_1
      sigma_1=self%sigma_array%row(j+1)%column(i_0)+&
         &(self%sigma_array%row(j+1)%column(i_1)-self%sigma_array%row(j+1)%column(i_0))*&
         &(r_mod-self%r_array(i_0))*self%inv_dr

      ! third interpolation
         
      ! get sigma value
      sigma_value=sigma_0+&
         &(sigma_1-sigma_0)*&
         &(b_mod-self%r_array(j))*self%inv_dr
      
      return
   
   end function
   
end module