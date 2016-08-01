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


! particle module
module m_particle

   use m_kind_parameters
   use m_maths
   
   implicit none
   
   private
   public :: particle
   
   ! define particle class
   type :: particle
   
      real(kind=rel_kind) :: r(n_dim)           ! position of particle
      real(kind=rel_kind) :: v(n_dim)           ! velocity of particle
      real(kind=rel_kind) :: m                  ! mass of particle
      real(kind=rel_kind) :: inv_h              ! 1/h
      real(kind=rel_kind) :: inv_rho            ! 1/rho
      real(kind=rel_kind) :: h                  ! smoothing length of particle
      real(kind=rel_kind) :: rho                ! density of particle
      real(kind=rel_kind) :: a_dot              ! mass weighted energy absorbtion rate
      real(kind=rel_kind) :: a_dot_new          ! new a_dot
      real(kind=rel_kind),pointer,contiguous :: lambda_array(:)      ! wavelength array
      real(kind=rel_kind),pointer,contiguous :: a_dot_scatter_array(:)  ! scattered light intensity bins
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: destroy
      procedure,non_overridable :: normalise_a
      procedure,non_overridable :: reset_a_dot_scatter
      procedure,non_overridable :: a_dot_scatter
   
   end type
   
   contains
   
   ! initialise particle
   pure subroutine initialise(self,r,v,m,a_dot,lambda_min,lambda_max,n_bins)
   
      ! argument declarations
      class(particle),intent(inout) :: self     ! particle object
      real(kind=rel_kind),intent(in) :: r(n_dim)! position
      real(kind=rel_kind),intent(in) :: v(n_dim)! velocity
      real(kind=rel_kind),intent(in) :: m       ! mass
      real(kind=rel_kind),intent(in) :: a_dot   ! energy absorption rate per unit mass
      real(kind=rel_kind),intent(in),optional :: lambda_min     ! minimum wavelength
      real(kind=rel_kind),intent(in),optional :: lambda_max     ! maximum wavelength
      integer(kind=int_kind),intent(in),optional :: n_bins      ! number of intensity bins
      
      ! set quantities
      self%r=r
      self%v=v
      self%m=m
      self%h=0._rel_kind
      self%a_dot=a_dot
      self%a_dot_new=0._rel_kind
      
      if (present(lambda_min).and.present(lambda_max).and.present(n_bins)) then
         allocate(self%lambda_array(n_bins))
         allocate(self%a_dot_scatter_array(n_bins-1))
         self%lambda_array=log_lin_space(lambda_min,lambda_max,n_bins)
         self%a_dot_scatter_array=0._rel_kind
      else
         self%lambda_array=>null()
         self%a_dot_scatter_array=>null()
      end if
      
      return
     
   end subroutine
   
   elemental subroutine destroy(self)
   
      ! argument declarations
      class(particle),intent(inout) :: self     ! particle object
      
      if (associated(self%lambda_array)) deallocate(self%lambda_array)
      if (associated(self%a_dot_scatter_array)) deallocate(self%a_dot_scatter_array)
   
      return
   
   end subroutine
   
   ! normalise a_dot and a_dot_scatter and also calc temperature
   elemental subroutine normalise_a(self)
   
      ! argument declarations
      class(particle),intent(inout) :: self  ! particle object
      
      self%a_dot=self%a_dot_new/self%m
      self%a_dot_new=0._rel_kind
      
      if (associated(self%a_dot_scatter_array)) then
      
         self%a_dot_scatter_array=self%a_dot_scatter_array/&
            &((self%lambda_array(2:)-self%lambda_array(:size(self%lambda_array)-1))*self%m)
      
      end if
      
      return
   
   end subroutine
   
   elemental subroutine reset_a_dot_scatter(self)
   
      ! argument declarations
      class(particle),intent(inout) :: self        ! particle object
      
      if (associated(self%a_dot_scatter_array)) self%a_dot_scatter_array=0._rel_kind
   
   end subroutine
   
   elemental function a_dot_scatter(self,lambda) result (value)
   
      ! argument declarations
      class(particle),intent(in) :: self        ! particle object
      real(kind=rel_kind),intent(in) :: lambda  ! wavelength
      
      ! result declaration
      real(kind=rel_kind) :: value              ! a_dot_scatter value
      
   
      value=0._rel_kind
      if (associated(self%lambda_array)) then
      
         if (lambda>self%lambda_array(1).and.lambda<self%lambda_array(size(self%lambda_array))) then

            value=value+self%a_dot_scatter_array(binary_search(lambda,self%lambda_array))    
               
         end if   
         
      end if
      
      return
   
   end function
   
end module