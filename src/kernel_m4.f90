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

! m4 cubic spline kernel
module m_kernel_m4

   use m_kind_parameters
   use m_sph_parameters
   use m_constants_parameters
   use m_kernel
   
   implicit none
   
   private
   public :: kernel_m4
   
   type,extends(kernel) :: kernel_m4
   
      private
   
      contains
      
      procedure :: initialise
      procedure,non_overridable,nopass,private :: w_exact
      procedure,non_overridable,nopass,private :: dw_dr_exact
      procedure,non_overridable,nopass,private :: dw_dh_exact
      procedure,non_overridable,nopass,private :: sigma_exact
   
   end type

   contains
   
   ! initialise m4 kernel
   pure subroutine initialise(self)
   
      ! argument declarations
      class(kernel_m4),intent(inout) :: self                    ! m4 kernel object
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      ! set support radius
      self%r_support=2._rel_kind
      
      ! set r interval
      self%inv_dr=real(n_lookup-1,rel_kind)/self%r_support
      
      ! set normalisation constant
      select case (n_dim)
         case(1)
            self%norm=2._rel_kind/3._rel_kind
         case(2)
            self%norm=10._rel_kind/(7._rel_kind*pi)
         case(3)
            self%norm=1._rel_kind/pi
      end select
      
      ! allocate arrays
      allocate(self%r_array(n_lookup))
      allocate(self%w_array(n_lookup))
      allocate(self%dw_dr_array(n_lookup))
      allocate(self%dw_dh_array(n_lookup))
      call self%sigma_array%initialise(n_lookup,n_lookup,.true._log_kind)
      
      
      ! populate arrays
      self%r_array=(/(real(i,rel_kind),i=0,n_lookup-1)/)/self%inv_dr
      
      self%w_array=self%norm*self%w_exact(self%r_array)
      self%dw_dr_array=self%norm*self%dw_dr_exact(self%r_array)
      self%dw_dh_array=self%norm*self%dw_dh_exact(self%r_array)
      
      do i=1,n_lookup
            
         self%sigma_array%row(i)%column=self%norm*&
            &self%sigma_exact(self%r_array(i),self%r_array(i:))

      end do
      
      return
   
   end subroutine
   
   ! kernel density
   elemental function w_exact(r) result (w_value)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: w_value                            ! kernel value
      
      w_value=&
         &0.25_rel_kind*max(2._rel_kind-r,0._rel_kind)**3-&
         &max(1._rel_kind-r,0._rel_kind)**3
      
      return
   
   end function
   
   ! kernel density gradient wrt r
   elemental function dw_dr_exact(r) result (dw_dr_value)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: dw_dr_value                        ! kernel value
      
      dw_dr_value=&
         &-0.75_rel_kind*max(2._rel_kind-r,0._rel_kind)**2+&
         &3._rel_kind*max(1._rel_kind-r,0._rel_kind)**2
      
      return
   
   end function
   
   ! kernel density gradient wrt h
   elemental function dw_dh_exact(r) result (dw_dh_value)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: dw_dh_value                        ! kernel value
      
      dw_dh_value=&
         &1.5_rel_kind*(r-1._rel_kind)*max(2._rel_kind-r,0._rel_kind)**2-&
         &3._rel_kind*(2._rel_kind*r-1._rel_kind)*max(1._rel_kind-r,0._rel_kind)**2
      
      return
   
   end function
   
   ! kernel column density
   elemental function sigma_exact(b,r) result(sigma_value)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: b                       ! impact parameter
      real(kind=rel_kind),intent(in) :: r                       ! radius
      
      ! result declaration
      real(kind=rel_kind) :: sigma_value                        ! column density
      
      ! variable declarations
      real(kind=rel_kind) :: b_mod                              ! modified impact parameter
      real(kind=rel_kind) :: r_mod                              ! modified radius
   
      ! integrate through outer kernel
      b_mod=max(min(b,2._rel_kind),epsilon(0._rel_kind))
      r_mod=max(min(r,2._rel_kind),b_mod)
      
      sigma_value=&
         &sqrt(r_mod**2-b_mod**2)*&
         &(-0.0625_rel_kind*r_mod**3+&
         &0.5_rel_kind*r_mod**2-&
         &0.09375_rel_kind*(b_mod**2+16._rel_kind)*r_mod+&
         &b_mod**2+2._rel_kind)-&
         &0.09375_rel_kind*(b_mod**2+16._rel_kind)*b_mod**2*&
         &log((r_mod+sqrt(r_mod**2-b_mod**2))/b_mod)
         
      ! integrate through inner kernel
      b_mod=max(min(b,1._rel_kind),epsilon(0._rel_kind))
      r_mod=max(min(r,1._rel_kind),b_mod)
      
      sigma_value=sigma_value-&
         &sqrt(r_mod**2-b_mod**2)*&
         &(-0.25_rel_kind*r_mod**3+&
         &r_mod**2-&
         &0.375_rel_kind*(b_mod**2+4._rel_kind)*r_mod+&
         &2._rel_kind*b_mod**2+1._rel_kind)+&
         &0.375_rel_kind*(b_mod**2+4._rel_kind)*b_mod**2*&
         &log((r_mod+sqrt(r_mod**2-b_mod**2))/b_mod)
         
      return
   
   end function
    
end module