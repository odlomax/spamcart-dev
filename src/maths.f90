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

! module containing some useful maths function and simple computer science procedures
module m_maths

   use m_kind_parameters
   use m_constants_parameters

   implicit none

   ! all routines are public
   public

   ! generic interfaces
   
   interface quicksort
      procedure :: quicksort_1
      procedure :: quicksort_2
   end interface
   
   interface swap
      procedure :: swap_real
      procedure :: swap_integer
   end interface
   
   interface lerp
      procedure :: lerp_1
      procedure :: lerp_2
   end interface
   
   interface gerp
      procedure :: gerp_1
      procedure :: gerp_2
   end interface
   
   interface lookup_and_interpolate
      procedure :: lookup_and_interpolate_1
      procedure :: lookup_and_interpolate_2
      procedure :: lookup_and_interpolate_3
   end interface
   
   interface lookup_and_geo_interpolate
      procedure :: lookup_and_geo_interpolate_1
      procedure :: lookup_and_geo_interpolate_2
   end interface
   
   interface rotate_vector
      procedure :: rotate_vector_3d
   end interface

   contains
   
!-------------------------------------------------------------
!-------------------------------------------------------------

! MATHEMATICAL FUNCTIONS


!-------------------------------------------------------------
!-------------------------------------------------------------

   ! return the nth real root of x
   elemental function nth_rt(x,n) result (nth_rt_x)
   
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x       ! variate
      integer(kind=int_kind),intent(in) :: n    ! nth root
      
      ! result declaration
      real(kind=rel_kind) :: nth_rt_x           ! nth root of x
      
      ! different result if number is odd or even
      if (mod(n,2_int_kind)==0) then
      
         ! return NaN if x is negative
         nth_rt_x=x**(1._rel_kind/real(n,rel_kind))
      
      else
      
         ! return real number if x is negative
         nth_rt_x=sign(abs(x)**(1._rel_kind/real(n,rel_kind)),x)
      
      end if
      
      return
   
   end function
   
   ! returns the real cube root of real number
   elemental function real_cbrt(a)
   
      real(kind=rel_kind) :: real_cbrt                   ! cube root of a
      real(kind=rel_kind),intent(in) :: a
      
      real_cbrt=sign(abs(a)**one_third,a)
      
      return
   
   end function
   
   ! return roots of quadratic equation
   pure function quadratic(a,b,c)                
   
      ! argument declarations
      complex(kind=cpx_kind) :: quadratic(2)          ! roots
      real(kind=rel_kind),intent(in) :: a             ! ax**2+bx+c=0
      real(kind=rel_kind),intent(in) :: b
      real(kind=rel_kind),intent(in) :: c
      
      ! variable declarations
      real(kind=rel_kind) :: inv2a                    ! some useful quantities
      complex(kind=cpx_kind) :: sqrtdelta
      
      !calculate useful quantities
      inv2a=0.5_rel_kind/a
      sqrtdelta=sqrt(cmplx(b**2-4._rel_kind*a*c,0._rel_kind,cpx_kind))
      
      !calculate roots
      
      quadratic(1)=-(b+sqrtdelta)*inv2a
      quadratic(2)=-(b-sqrtdelta)*inv2a
      
      return
   
   end function
   
   ! return roots of cubic equation (method from Wolfram Alpha)
   pure function cubic(a,b,c,d)                       
   
      complex(kind=cpx_kind) :: cubic(3)              ! roots
      
      real(kind=rel_kind),intent(in) :: a             ! ax**3+bx**2+cx+d=0
      real(kind=rel_kind),intent(in) :: b
      real(kind=rel_kind),intent(in) :: c
      real(kind=rel_kind),intent(in) :: d
   
      real(kind=rel_kind) :: a2                       ! some useful quantities
      real(kind=rel_kind) :: a1
      real(kind=rel_kind) :: a0
      real(kind=rel_kind) :: inva
      real(kind=rel_kind) :: big_q
      real(kind=rel_kind) :: big_r
      real(kind=rel_kind) :: big_d
      real(kind=rel_kind) :: big_s
      real(kind=rel_kind) :: big_t
      real(kind=rel_kind) :: sqrtbig_d
      real(kind=rel_kind) :: sqrtbig_q
      real(kind=rel_kind) :: theta
      
      !calculate useful quantities
      
      inva=1._rel_kind/a
      a2=b*inva
      a1=c*inva
      a0=d*inva
      
      big_q=(3._rel_kind*a1-a2**2)*(1._rel_kind/9._rel_kind)
      big_r=(9._rel_kind*a2*a1-27._rel_kind*a0-2._rel_kind*a2**3)*(1._rel_kind/54._rel_kind)
      
      big_d=big_q**3+big_r**2
      
      !calculate roots
      
      if (big_d<0._rel_kind) then                      ! three distinct real roots
         sqrtbig_q=sqrt(abs(big_q))
         theta=acos(big_r/sqrtbig_q**3)
         cubic(1)=2._rel_kind*sqrtbig_q*cos(theta*one_third)-one_third*a2
         cubic(2)=2._rel_kind*sqrtbig_q*cos((theta+2._rel_kind*pi)*one_third)-one_third*a2
         cubic(3)=2._rel_kind*sqrtbig_q*cos((theta+4._rel_kind*pi)*one_third)-one_third*a2
      else                                            ! two identical or complex roots
         sqrtbig_d=sqrt(big_d)
         big_s=real_cbrt(big_r+sqrtbig_d)
         big_t=real_cbrt(big_r-sqrtbig_d)
         cubic(1)=-one_third*a2+(big_s+big_t)
         cubic(2)=&
            &-one_third*a2-0.5_rel_kind*(big_s+big_t)+0.5_rel_kind*(0._rel_kind,1._rel_kind)*root_three*(big_s-big_t)
         cubic(3)=&
            &-one_third*a2-0.5_rel_kind*(big_s+big_t)-0.5_rel_kind*(0._rel_kind,1._rel_kind)*root_three*(big_s-big_t)
      end if
   
      return
   
   end function
   
   ! return roots of quartic equation (method from Wolfram Alpha)
   pure function quartic(a,b,c,d,e)
   
      complex(kind=cpx_kind) :: quartic(4)                  ! roots
      
      real(kind=rel_kind),intent(in) :: a                   ! ax**4+bx**3+cx**2+dx+e=0
      real(kind=rel_kind),intent(in) :: b
      real(kind=rel_kind),intent(in) :: c
      real(kind=rel_kind),intent(in) :: d
      real(kind=rel_kind),intent(in) :: e
      
      real(kind=rel_kind) :: a3                             ! some useful quantites
      real(kind=rel_kind) :: a2
      real(kind=rel_kind) :: a1
      real(kind=rel_kind) :: a0
      real(kind=rel_kind) :: inva
      
      complex(kind=cpx_kind) :: y(3)
      complex(kind=cpx_kind) :: big_r
      complex(kind=cpx_kind) :: big_d
      complex(kind=cpx_kind) :: big_e
      complex(kind=cpx_kind) :: y1term
      
      !calculate useful quantities
      
      inva=1._rel_kind/a
      a3=b*inva
      a2=c*inva
      a1=d*inva
      a0=e*inva
      
      y=cubic(1._rel_kind,-a2,a1*a3-4._rel_kind*a0,4._rel_kind*a2*a0-a1**2-a3**2*a0)
      
      big_r=sqrt(0.25_rel_kind*a3**2-a2+y(1))
      
      if (abs(big_r)>0._rel_kind) then
         y1term=0.25_rel_kind*(4._rel_kind*a3*a2-8._rel_kind*a1-a3**3)/big_r
         big_d=sqrt(0.75_rel_kind*a3**2-big_r**2-2._rel_kind*a2+y1term)
         big_e=sqrt(0.75_rel_kind*a3**2-big_r**2-2._rel_kind*a2-y1term)
      else
         y1term=2._rel_kind*sqrt(y(1)**2-4._rel_kind*a0)
         big_d=sqrt(0.75_rel_kind*a3**2-2._rel_kind*a2+y1term)
         big_e=sqrt(0.75_rel_kind*a3**2-2._rel_kind*a2-y1term)
      end if
      
      !calculate roots
      
      quartic(1)=-0.25_rel_kind*a3+0.5_rel_kind*big_r+0.5_rel_kind*big_d
      quartic(2)=-0.25_rel_kind*a3+0.5_rel_kind*big_r-0.5_rel_kind*big_d
      quartic(3)=-0.25_rel_kind*a3-0.5_rel_kind*big_r+0.5_rel_kind*big_e
      quartic(4)=-0.25_rel_kind*a3-0.5_rel_kind*big_r-0.5_rel_kind*big_e
      
      return
   
   end function

   
!-------------------------------------------------------------
!-------------------------------------------------------------

! RANDOM NUMBERS

!-------------------------------------------------------------
!-------------------------------------------------------------
   
   
   subroutine gaussian_random_number(r_num)

      ! result declaration
      real(kind=rel_kind),intent(out) :: r_num(2)  ! gaussian random numbers
      
      ! variable declarations
      real(kind=rel_kind) :: w                     ! polar magnitude sqrd

      do
         call random_number(r_num)
         r_num=(r_num*2._rel_kind)-1._rel_kind
         w=sum(r_num**2)
         if (w<=1._rel_kind) exit
      end do

      w=sqrt(-2._rel_kind*log(w)/w)
      r_num=r_num*w

      return

   end subroutine
   
   subroutine gaussian_random_array(r_num)
   
      ! argument declaration
      real(kind=rel_kind),intent(out) :: r_num(:)       ! gaussian random number
      
      ! variable declarations
      integer(kind=int_kind) :: i                       ! counter
      real(kind=rel_kind) :: r_num_temp(2)
      
      do i=1,size(r_num)/2
      
         call gaussian_random_number(r_num_temp)
         r_num(2*i-1:2*i)=r_num_temp
         
      end do
      
      if (mod(size(r_num),2)==1) then
      
         call gaussian_random_number(r_num_temp)
         r_num(size(r_num))=r_num_temp(1)
         
      end if
      
      return
   
   end subroutine
   
   ! calculates the unit vector of a random direction 
   subroutine random_direction(direction)

      ! argument declarations
      real(kind=rel_kind),intent(out) :: direction(:)

      ! variable declarations
      call gaussian_random_array(direction)
      direction=direction/norm2(direction)
      
      return

   end subroutine
   
   
   ! calculates the unit vector of a direction scattered from a normal by a cosine
   subroutine scattered_direction_3d(direction,norm,cos_theta)

      ! argument declarations
      real(kind=rel_kind),intent(out) :: direction(3)   ! direction vector
      real(kind=rel_kind),intent(in) :: norm(3)         ! normal vector
      real(kind=rel_kind),intent(in) :: cos_theta       ! cosine between norm and direction
      
      ! variable declaration
      real(kind=rel_kind) :: r_num                      ! random number
      real(kind=rel_kind) :: sin_theta
      real(kind=rel_kind) :: phi
      real(kind=rel_kind) :: cos_theta_s                ! norm angles
      real(kind=rel_kind) :: sin_theta_s
      real(kind=rel_kind) :: phi_s
      real(kind=rel_kind) :: cos_phi_s
      real(kind=rel_kind) :: sin_phi_s
      real(kind=rel_kind) :: r_matrix(3,3)              ! rotation matrix
      
      sin_theta=sqrt(1._rel_kind-cos_theta**2)
      
      ! draw phi from interval [0:2pi]
      call random_number(r_num)
      phi=2._rel_kind*pi*r_num
      
      ! make direction vector
      direction(1)=sin_theta*cos(phi)
      direction(2)=sin_theta*sin(phi)
      direction(3)=cos_theta
      
      ! rotate direction into plane of norm
      cos_theta_s=norm(3)
      sin_theta_s=sin(acos(cos_theta_s))
      phi_s=atan2(norm(2),norm(1))
      cos_phi_s=cos(phi_s)
      sin_phi_s=sin(phi_s)
      
      r_matrix(1,:)=(/cos_phi_s*cos_theta_s,-sin_phi_s,cos_phi_s*sin_theta_s/)
      r_matrix(2,:)=(/sin_phi_s*cos_theta_s,cos_phi_s,sin_phi_s*sin_theta_s/)
      r_matrix(3,:)=(/-sin_theta_s,0._rel_kind,cos_theta_s/)
      
      direction=matmul(r_matrix,direction) 

      return

   end subroutine
   
!-------------------------------------------------------------
!-------------------------------------------------------------

! PHYSICAL QUANTITIES

!-------------------------------------------------------------
!-------------------------------------------------------------
   
   ! Planck function [erg s^-1 sr^-1 cm^-2 micron^-1]
   elemental function planck_lambda(lambda,t) result(b)
   
      ! argument delcarations
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength in microns
      real(kind=rel_kind),intent(in) :: t                     ! temperature in K
      
      ! result declaration
      real(kind=rel_kind) :: b                                ! planck function
                                                        ! [erg s^-1 sr^-1 cm^-2 micron^-1]
      ! variable declaration
      real(kind=rel_kind) :: lambda_si                        ! wavelength in m
      real(kind=rel_kind) :: exp_part                   ! exponential part
      real(kind=rel_kind) :: denominator
      
      lambda_si=1e-6_rel_kind*lambda
      exp_part=exp(h_planck*c_light/(lambda_si*k_boltz*t))
   
      if (exp_part<huge(0._rel_kind)) then
         ! calculate planck function
         
         denominator=lambda_si**5*(exp_part-1._rel_kind)
         
         if (denominator<huge(0._rel_kind)) then
         
            b=2._rel_kind*h_planck*c_light**2/denominator
            
         else
         
            b=0._rel_kind
            
         end if
               
      else
      
         b=0._rel_kind
         
      end if
         
      ! normalize to cgs
      b=1e-3_rel_kind*b
         
      return
   
   end function
   
!-------------------------------------------------------------
!-------------------------------------------------------------

! INTERPOLATION, INTEGRATION, EXTRAPOLATION, ETC

!-------------------------------------------------------------
!-------------------------------------------------------------

   ! 1D lookup and interpolate
   pure function lookup_and_interpolate_1(x,x_array,f_array) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in),contiguous :: x_array(:) ! array of x_values
      real(kind=rel_kind),intent(in),contiguous :: f_array(:) ! array of f(x)
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x)
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation indices
      real(kind=rel_kind) :: x_adj                            ! adjusted value of x
      
      ! enforce that x is within range of x_array
      x_adj=min(x,x_array(size(x_array)))
      x_adj=max(x_adj,x_array(1))
      
      
      ! lookup indices
      i_0=binary_search(x,x_array)

      ! interpolate for f
      f=lerp(x_adj,x_array(i_0),x_array(i_0+1),f_array(i_0),f_array(i_0+1))
      
      return
   
   end function
   
   ! 1D lookup and interpolate two quantities
   pure function lookup_and_interpolate_2(x,x_array,f_array_1,f_array_2) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in),contiguous :: x_array(:) ! array of x_values
      real(kind=rel_kind),intent(in),contiguous :: f_array_1(:) ! array of f(x)
      real(kind=rel_kind),intent(in),contiguous :: f_array_2(:) ! array of f(x)
      
      ! result declaration
      real(kind=rel_kind) :: f(2)                             ! f(x)
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation indices
      real(kind=rel_kind) :: x_adj                            ! adjusted value of x
      
      ! enforce that x is within range of x_array
      x_adj=min(x,x_array(size(x_array)))
      x_adj=max(x_adj,x_array(1))
      
      
      ! lookup indices
      i_0=binary_search(x,x_array)

      ! interpolate for f
      f=lerp(x_adj,x_array(i_0),x_array(i_0+1),&
         &(/f_array_1(i_0),f_array_2(i_0)/),(/f_array_1(i_0+1),f_array_2(i_0+1)/))
      
      return
   
   end function
   
   ! lookup and inverse quad interpolate
   pure function lookup_and_inv_interpolate(f,f_array,x_array,dfdx_array) result (x)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: f                     ! value of f(x)
      real(kind=rel_kind),intent(in),contiguous :: f_array(:) ! array of f(x)
      real(kind=rel_kind),intent(in),contiguous :: x_array(:) ! array of x
      real(kind=rel_kind),intent(in),contiguous :: dfdx_array(:)        ! array of df(x)/dx
      
      ! result declaration
      real(kind=rel_kind) :: x                                ! x
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation indices
      real(kind=rel_kind) :: f_adj                            ! adjusted value of f(x)
      
      ! enforce that x is within range of x_array
      f_adj=min(f,f_array(size(f_array)))
      f_adj=max(f_adj,f_array(1))
      
      ! lookup indices
      i_0=binary_search(f,f_array)
      
      ! interpolate for f
      x=inv_querp(f_adj,f_array(i_0),f_array(i_0+1),x_array(i_0),x_array(i_0+1),&
         &dfdx_array(i_0),dfdx_array(i_0+1))
      
      return
   
   end function
   
   ! 2D lookup and interpolate
   pure function lookup_and_interpolate_3(x,y,x_array,y_array,f_array) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in) :: y                     ! value of y
      real(kind=rel_kind),intent(in),contiguous :: x_array(:) ! array of x values
      real(kind=rel_kind),intent(in),contiguous :: y_array(:) ! array of y values
      real(kind=rel_kind),intent(in),contiguous :: f_array(:,:) ! array of f(x,y)
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x,y)
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation indices
      integer(kind=int_kind) :: j_0
      real(kind=rel_kind) :: x_adj                            ! adjusted value of x
      real(kind=rel_kind) :: y_adj                            ! adjusted value of y
      
      ! enforce that x and y are within range
      x_adj=min(x,x_array(size(x_array)))
      x_adj=max(x_adj,x_array(1))
      y_adj=min(y,y_array(size(y_array)))
      y_adj=max(y_adj,y_array(1))
      
      ! lookup indices
      i_0=binary_search(x_adj,x_array)
      j_0=binary_search(y_adj,y_array)
      
      ! interpolate for f
      f=lerp(x_adj,y_adj,x_array(i_0),x_array(i_0+1),y_array(j_0),y_array(j_0+1),&
         &f_array(i_0,j_0),f_array(i_0,j_0+1),f_array(i_0+1,j_0),f_array(i_0+1,j_0+1))
      
      return
   
   end function
   
   ! look and geometically interpolate quantity in one dimension
   pure function lookup_and_geo_interpolate_1(x,x_array,f_array) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in),contiguous :: x_array(:) ! array of x_values
      real(kind=rel_kind),intent(in),contiguous :: f_array(:) ! array of f(x)
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x)
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation indices
      real(kind=rel_kind) :: x_adj                            ! adjusted value of x
      
      ! enforce that x is within range of x_array
      x_adj=min(x,x_array(size(x_array)))
      x_adj=max(x_adj,x_array(1))
      
      
      ! lookup indices
      i_0=binary_search(x,x_array)

      ! interpolate for f
      f=gerp(x_adj,x_array(i_0),x_array(i_0+1),f_array(i_0),f_array(i_0+1))
      
      return
      
   end function
   
   ! look and geometically interpolate quantity in two dimensions
   pure function lookup_and_geo_interpolate_2(x,y,x_array,y_array,f_array) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in) :: y                     ! value of y
      real(kind=rel_kind),intent(in),contiguous :: x_array(:) ! array of x values
      real(kind=rel_kind),intent(in),contiguous :: y_array(:) ! array of y values
      real(kind=rel_kind),intent(in),contiguous :: f_array(:,:) ! array of f(x,y)
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x,y)
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation indices
      integer(kind=int_kind) :: j_0
      real(kind=rel_kind) :: x_adj                            ! adjusted value of x
      real(kind=rel_kind) :: y_adj                            ! adjusted value of y
      
      ! enforce that x and y are within range
      x_adj=min(x,x_array(size(x_array)))
      x_adj=max(x_adj,x_array(1))
      y_adj=min(y,y_array(size(y_array)))
      y_adj=max(y_adj,y_array(1))
      
      ! lookup indices
      i_0=binary_search(x_adj,x_array)
      j_0=binary_search(y_adj,y_array)
      
      ! interpolate for f
      f=gerp(x_adj,y_adj,x_array(i_0),x_array(i_0+1),y_array(j_0),y_array(j_0+1),&
         &f_array(i_0,j_0),f_array(i_0,j_0+1),f_array(i_0+1,j_0),f_array(i_0+1,j_0+1))
      
      return
   
   end function
  
   ! linear interpolation 
   elemental function lerp_1(x,x_0,x_1,f_0,f_1) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in) :: x_0                   ! x interpolation points
      real(kind=rel_kind),intent(in) :: x_1
      real(kind=rel_kind),intent(in) :: f_0                   ! f(x) interpolation points
      real(kind=rel_kind),intent(in) :: f_1
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x)
      
      ! solve for f
      f=f_0+(f_1-f_0)*(x-x_0)/(x_1-x_0)
      
      return
   
   end function
   
   ! bilinear interpolation
   ! three linear interpolations
   elemental function lerp_2(x,y,x_0,x_1,y_0,y_1,f_00,f_01,f_10,f_11) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in) :: y                     ! value of y
      real(kind=rel_kind),intent(in) :: x_0                   ! x interpolation points
      real(kind=rel_kind),intent(in) :: x_1
      real(kind=rel_kind),intent(in) :: y_0                   ! y interpolation points
      real(kind=rel_kind),intent(in) :: y_1
      real(kind=rel_kind),intent(in) :: f_00                  ! f(x,y) values
      real(kind=rel_kind),intent(in) :: f_01
      real(kind=rel_kind),intent(in) :: f_10
      real(kind=rel_kind),intent(in) :: f_11
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x,y)
      
      ! solve for f
      f=lerp(y,y_0,y_1,lerp(x,x_0,x_1,f_00,f_10),lerp(x,x_0,x_1,f_01,f_11))
      
      return
      
   end function
   
   ! geometic interpolation
   elemental function gerp_1(x,x_0,x_1,f_0,f_1) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in) :: x_0                   ! x interpolation points
      real(kind=rel_kind),intent(in) :: x_1
      real(kind=rel_kind),intent(in) :: f_0                   ! f(x) interpolation points
      real(kind=rel_kind),intent(in) :: f_1
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x)
      
      ! variable declarations
      real(kind=rel_kind) :: weight                           ! (x-x_0)/(x_1-x_0)
      
      ! solve for f
      
      weight=(x-x_0)/(x_1-x_0)
      f=f_0**(1._rel_kind-weight)*f_1**(weight)
      
      return
   
   end function
   
   ! bi-geometic interpolation
   ! three geometric interpolations
   elemental function gerp_2(x,y,x_0,x_1,y_0,y_1,f_00,f_01,f_10,f_11) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in) :: y                     ! value of y
      real(kind=rel_kind),intent(in) :: x_0                   ! x interpolation points
      real(kind=rel_kind),intent(in) :: x_1
      real(kind=rel_kind),intent(in) :: y_0                   ! y interpolation points
      real(kind=rel_kind),intent(in) :: y_1
      real(kind=rel_kind),intent(in) :: f_00                  ! f(x,y) values
      real(kind=rel_kind),intent(in) :: f_01
      real(kind=rel_kind),intent(in) :: f_10
      real(kind=rel_kind),intent(in) :: f_11
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x,y)
      
      ! solve for f
      f=gerp(y,y_0,y_1,gerp(x,x_0,x_1,f_00,f_10),gerp(x,x_0,x_1,f_01,f_11))
      
      return
      
   end function
   
   ! linear interpolation (with x forced into interval)
   elemental function lerp_forced(x,x_0,x_1,f_0,f_1) result (f)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                     ! value of x
      real(kind=rel_kind),intent(in) :: x_0                   ! x interpolation points
      real(kind=rel_kind),intent(in) :: x_1
      real(kind=rel_kind),intent(in) :: f_0                   ! f(x) interpolation points
      real(kind=rel_kind),intent(in) :: f_1
      
      ! result declaration
      real(kind=rel_kind) :: f                                ! f(x)
      
      ! solve for f
      f=f_0+(f_1-f_0)*(min(max(x,x_0),x_1)-x_0)/(x_1-x_0)
      
      return
   
   end function
   
   ! inverse quadratic interpolation 
   elemental function inv_querp(f,f_0,f_1,x_0,x_1,dfdx_0,dfdx_1) result(x)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: f                     ! value of f(x)
      real(kind=rel_kind),intent(in) :: f_0                   ! f(x) interpolation points
      real(kind=rel_kind),intent(in) :: f_1
      real(kind=rel_kind),intent(in) :: x_0                   ! x interpolation points
      real(kind=rel_kind),intent(in) :: x_1
      real(kind=rel_kind),intent(in) :: dfdx_0                ! df(x)/dx intepolation points
      real(kind=rel_kind),intent(in) :: dfdx_1

      ! result declaration
      real(kind=rel_kind) :: x                                ! x
      
      ! variable declarations
      real(kind=rel_kind) :: delta_f                          ! f-f_0
      real(kind=rel_kind) :: delta_x                          ! x-x_0
      real(kind=rel_kind) :: dfdx                             ! simple df(x)/dx
      real(kind=rel_kind) :: d2fdx2                           ! d^2f(x)/dx^2
      real(kind=rel_kind) :: discriminant                     ! dfdx_0**2+2._rel_kind*d2fdx2*delta_f
      
      delta_f=f_0-f
      dfdx=(f_1-f_0)/(x_1-x_0)
      d2fdx2=(dfdx_1-dfdx_0)/(x_1-x_0)
      discriminant=dfdx_0**2-2._rel_kind*delta_f*d2fdx2
      
      if (discriminant>0._rel_kind) then
      
         ! Halley's irrational formula
         delta_x=-2._rel_kind*delta_f/(dfdx_0+sign(sqrt(discriminant),dfdx_0))
         
      else if (abs(dfdx)>0._rel_kind) then
      
         ! Newton Raphson
         delta_x=-delta_f/dfdx
         
      else
      
         delta_x=0._rel_kind
      
      end if
      
      x=x_0+delta_x
        
      return
      
   end function

   ! makes a 1D histogram from a list of values and and array of bins
   pure function make_hist(values,bins,weight_in) result (hist)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: values(:)             ! array of values
      real(kind=rel_kind),intent(in) :: bins(:)               ! array of bin values
      real(kind=rel_kind),intent(in),optional :: weight_in(:) ! weight of values
      
      ! result declaration
      real(kind=rel_kind) :: hist(size(bins)-1)               ! histogram values
      
      ! variable declarations
      integer(kind=int_kind) :: i                             ! counter
      integer(kind=int_kind) :: j
      real(kind=rel_kind) :: weight(size(values))             ! weight of values
      
      if (present(weight_in)) then
         weight=weight_in
      else
         weight=1._rel_kind
      end if
      
      ! initialise histogram
      hist=0._rel_kind
      
      ! loop over values and update histogram
      do i=1,size(values)
      
         j=binary_search(values(i),bins)
         hist(j)=hist(j)+weight(i)
      
      end do
      
      ! normalise bins
      hist=hist/(real(size(values),rel_kind)*(bins(2:)-bins(:size(bins)-1)))
   
      return
   
   end function
   
   ! produce an array of linearly spaced values
   pure function lin_space(x_min,x_max,n_x) result(x)
   
      ! argument declarations
      integer(kind=int_kind),intent(in) :: n_x        ! number of values
      real(kind=rel_kind),intent(in) :: x_min         ! min value
      real(kind=rel_kind),intent(in) :: x_max         ! max value
      
      ! result declaration
      real(kind=rel_kind) :: x(n_x)                   ! arrray of n_x values
      
      ! variable declarations
      integer(kind=int_kind) :: i                     ! counter
      real(kind=rel_kind) :: dx                       ! x incriment
      
      dx=(x_max-x_min)/real(n_x-1,rel_kind)
      x=(/(x_min+real(i-1,rel_kind)*dx,i=1,n_x)/)
   
      return
   
   end function
   
   ! produce an array of logarithmically spaced values
   pure function log_lin_space(x_min,x_max,n_x) result(x)
      
      ! argument declarations
      integer(kind=int_kind),intent(in) :: n_x        ! number of values
      real(kind=rel_kind),intent(in) :: x_min         ! min value
      real(kind=rel_kind),intent(in) :: x_max         ! max value
      
      ! result declaration
      real(kind=rel_kind) :: x(n_x)                   ! arrray of n_x values
      
      ! variable declarations
      integer(kind=int_kind) :: i                     ! counter
      real(kind=rel_kind) :: log_x_min                ! min log x
      real(kind=rel_kind) :: log_x_max                ! max log x
      real(kind=rel_kind) :: dlog_x                   ! log x incriment
      
      log_x_min=log(x_min)
      log_x_max=log(x_max)
      
      dlog_x=(log_x_max-log_x_min)/real(n_x-1,rel_kind)
      x=(/(exp(log_x_min+real(i-1,rel_kind)*dlog_x),i=1,n_x)/)
   
      return 
   
   end function
   
   ! perform a trapezoidal integration
   pure function trapz_intgr(x_values,f_values) result (value)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x_values(:)   ! array of x
      real(kind=rel_kind),intent(in) :: f_values(:)   ! array of f(x)
      
      ! result declaration
      real(kind=rel_kind) :: value                    ! integral of f(x)dx
      
      ! variable declarations
      integer(kind=int_kind) :: n_x                   ! number of x_values
      
      n_x=size(x_values)
      
      value=0.5_rel_kind*sum((f_values(1:n_x-1)+f_values(2:n_x))*&
         &(x_values(2:n_x)-x_values(1:n_x-1)))
      
      return
   
   end function
   
   ! create cumulative dist of f(x)
   pure function cum_dist_func(x_values,f_values,normalise) result (values)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x_values(:)   ! array of x
      real(kind=rel_kind),intent(in) :: f_values(:)   ! array of f(x)
      logical(kind=log_kind),intent(in),optional :: normalise  ! normalise distribution to one
      
      ! result declaration
      real(kind=rel_kind) :: values(size(x_values))   ! integral of f(x)dx
      
      ! variable declaration
      integer(kind=int_kind) :: i                     ! counter
      integer(kind=int_kind) :: n_x                   ! number of x_values
      
      n_x=size(x_values)
      
      ! first element is zero by definition
      values(1)=0._rel_kind
      
      ! loop over array and calculate cumulative dist
      do i=2,n_x
         values(i)=values(i-1)+&
            &0.5_rel_kind*(f_values(i-1)+f_values(i))*&
            &(x_values(i)-x_values(i-1))
      end do
      ! normalise distribution to one
      if (present(normalise)) then
         if (normalise) values=values/values(n_x)
      end if
      
      return
   
   end function
   
   

!-------------------------------------------------------------
!-------------------------------------------------------------

! ARRAY SEARCHES, SORTS, ETC

!-------------------------------------------------------------
!-------------------------------------------------------------

   ! binary search
   ! find right-most element in monotonically increasing array with value less than x
   ! x must be in range array(i) <= x <= array(size(array))
   pure function binary_search(x,array) result(i_0)
   
      ! argument declations
      real(kind=rel_kind),intent(in) :: x                     ! target value
      real(kind=rel_kind),intent(in),contiguous :: array(:)   ! array of values
      
      ! result declaration
      integer(kind=int_kind) :: i_0                           ! right most lower index
      
      ! variable declarations
      integer(kind=int_kind) :: i_1                           ! left most upper index
      integer(kind=int_kind) :: i_mid                         ! mid point between i_0 and i_1
      
      ! set initial array bounds
      i_0=1
      i_1=size(array)
      
      ! find i_0
      do 
         ! set i_0 to array mid point
         i_mid=i_0+(i_1-i_0)/2
         
         ! exit or change array boundaries
         if (i_1-i_0<=1) then
            exit
         else if (array(i_mid)<=x) then
            i_0=i_mid
         else
            i_1=i_mid
         end if

      end do
      
   end function
   
   ! produce an array of integers equal to the indices of sorted array
   pure function array_ranking(array) result (indices)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: array(:)              ! array (obviously)
      
      ! result declaration
      integer(kind=int_kind) :: indices(size(array))          ! ranked array indices
      
      ! variable declarations
      integer(kind=int_kind) :: i                             ! counter
      real(kind=rel_kind) :: temp_array(size(array))
      
      ! make local copy of array
      temp_array=array
   
      ! initialise ranking
      indices=(/(i,i=1,size(indices))/)
      
      call quicksort(temp_array,indices)
      
      return
   
   end function
   
   ! quick sort algorithm
   pure recursive subroutine quicksort_1(array)

      !argument declarations
      real(kind=rel_kind),intent(inout) :: array(:)    !array to be sorted
      
      !variable declarations
      integer(kind=int_kind) :: i                      !counter
      integer(kind=int_kind) :: n                      !size of array
      integer(kind=int_kind) :: pivot                  !pivot index
      
      !Quicksort algorithm
      
      !return if array has one or zero elements
      n=size(array)
      if(n<=1) return

      !set pivot to central array element
      pivot=1+(n-1)/2
      
      !move pivot to end of array
      call swap(array(pivot),array(n))
      
      !loop over array and place elements
      !lower than the pivot value left of pivot
      pivot=1
      do i=1,n-1
      
         if(array(i)<=array(n)) then
            call swap(array(pivot),array(i))
            pivot=pivot+1
         end if
      
      end do
      
      !move pivot back from end of array
      call swap(array(pivot),array(n))
      
      !sort the array slices either side of pivot
      call quicksort(array(:pivot-1))
      call quicksort(array(pivot+1:))
         
      return

   end subroutine
   
   ! quick sort algorithm
   ! also maintains list of original indices
   pure recursive subroutine quicksort_2(array,indices)

      !argument declarations
      integer(kind=int_kind),intent(inout) :: indices(:)  !original indices of array elements
      real(kind=rel_kind),intent(inout) :: array(:)    !array to be sorted
      
      !variable declarations
      integer(kind=int_kind) :: i                      !counter
      integer(kind=int_kind) :: n                      !size of array
      integer(kind=int_kind) :: pivot                  !pivot index
      
      
      !Quicksort algorithm
      
      !return if array has one or zero elements
      n=size(array)
      if(n<=1) return

      !set pivot to central array element
      pivot=1+(n-1)/2
      
      !move pivot to end of array
      call swap(array(pivot),array(n))
      call swap(indices(pivot),indices(n))
      
      !loop over array and place elements
      !lower than the pivot value left of pivot
      pivot=1
      do i=1,n-1
      
         if(array(i)<=array(n)) then
            call swap(array(pivot),array(i))
            call swap(indices(pivot),indices(i))
            pivot=pivot+1
         end if
      
      end do
      
      !move pivot back from end of array
      call swap(array(pivot),array(n))
      call swap(indices(pivot),indices(n))
      
      !sort the array slices either side of pivot
      call quicksort(array(:pivot-1),indices(:pivot-1))
      call quicksort(array(pivot+1:),indices(pivot+1:))
         
      return

   end subroutine
   
   pure subroutine swap_real(a,b)
   
      ! argument declarations
      real(kind=rel_kind),intent(inout) :: a
      real(kind=rel_kind),intent(inout) :: b
      
      ! variable declarations
      real(kind=rel_kind) :: temp
      
      temp=a
      a=b
      b=temp
      
      return
   
   end subroutine
   
   pure subroutine swap_integer(a,b)
   
      ! argument declarations
      integer(kind=int_kind),intent(inout) :: a
      integer(kind=int_kind),intent(inout) :: b
      
      ! variable declarations
      integer(kind=int_kind) :: temp
      
      temp=a
      a=b
      b=temp
      
      return
   
   end subroutine
   
   pure function rotate_vector_3d(x,u,theta) result(x_new)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(3)       ! vector to rotate
      real(kind=rel_kind),intent(in) :: u(3)       ! axis of rotation unit vector
      real(kind=rel_kind),intent(in) :: theta      ! angle of rotation
      
      ! result declaration
      real(kind=rel_kind) :: x_new(3)              ! rotated x vector
      
      x_new=matmul(rotation_matrix_3d(u,theta),x)
      
      return   
   
   end function
   
   ! construct a 3d rotation matrix for an angle about a unit vector
   pure function rotation_matrix_3d(u,theta) result(r)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: u(3)    ! unit vector
      real(kind=rel_kind),intent(in) :: theta   ! angle of rotation
      
      ! result declaration
      real(kind=rel_kind) :: r(3,3)             ! rotation matrix
      
      ! variable declarations
      real(kind=rel_kind) :: cos_theta             ! cos(theta)
      real(kind=rel_kind) :: sin_theta             ! sin(theta)
      
      cos_theta=cos(theta)
      sin_theta=sin(theta)
      
      r(:,1)=(/cos_theta+u(1)**2*(1._rel_kind-cos_theta),&
         &u(2)*u(1)*(1._rel_kind-cos_theta)+u(3)*sin_theta,&
         &u(3)*u(1)*(1._rel_kind-cos_theta)-u(2)*sin_theta/)
      r(:,2)=(/u(1)*u(2)*(1._rel_kind-cos_theta)-u(3)*sin_theta,&
         &cos_theta+u(2)**2*(1._rel_kind-cos_theta),&
         &u(3)*u(2)*(1._rel_kind-cos_theta)+u(1)*sin_theta/)
      r(:,3)=(/u(1)*u(3)*(1._rel_kind-cos_theta)+u(2)*sin_theta,&
         &u(2)*u(3)*(1._rel_kind-cos_theta)-u(1)*sin_theta,&
         &cos_theta+u(3)**2*(1._rel_kind-cos_theta)/)
      
      return
   
   end function
   
 
!-------------------------------------------------------------
!-------------------------------------------------------------

! POINT DISTRIBUTIONS

!-------------------------------------------------------------
!-------------------------------------------------------------  

   ! cube of edge length 1 rotated through euler angles 
   pure function rotated_cube_vertices(angles) result(vertices)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: angles(3)  ! three Euler angles

      ! result declaration
      real(kind=rel_kind) :: vertices(3,8)        ! eight three-dimensional positions
      
      ! variable declarations
      integer(kind=int_kind) :: i                  ! counter
      real(kind=rel_kind) :: unit_vectors(3,3)     ! three unit vectors
      real(kind=rel_kind) :: rotation_matrix(3,3)       ! rotation matrix
      
      ! set initial vector/vertices
      vertices(:,1)=(/-0.5_rel_kind,-0.5_rel_kind,-0.5_rel_kind/)
      vertices(:,2)=(/-0.5_rel_kind,-0.5_rel_kind,0.5_rel_kind/)
      vertices(:,3)=(/-0.5_rel_kind,0.5_rel_kind,-0.5_rel_kind/)
      vertices(:,4)=(/-0.5_rel_kind,0.5_rel_kind,0.5_rel_kind/)
      vertices(:,5)=(/0.5_rel_kind,-0.5_rel_kind,-0.5_rel_kind/)
      vertices(:,6)=(/0.5_rel_kind,-0.5_rel_kind,0.5_rel_kind/)
      vertices(:,7)=(/0.5_rel_kind,0.5_rel_kind,-0.5_rel_kind/)
      vertices(:,8)=(/0.5_rel_kind,0.5_rel_kind,0.5_rel_kind/)
      unit_vectors(:,1)=(/1._rel_kind,0._rel_kind,0._rel_kind/)
      unit_vectors(:,2)=(/0._rel_kind,1._rel_kind,0._rel_kind/)
      unit_vectors(:,3)=(/0._rel_kind,0._rel_kind,1._rel_kind/)
      
      ! perform rotation about x axis
      rotation_matrix=rotation_matrix_3d(unit_vectors(:,1),angles(1))
      ! rotate vertices
      do i=1,8
         vertices(:,i)=matmul(rotation_matrix,vertices(:,i))
      end do
      ! rotate unit vectors
      unit_vectors(:,2)=matmul(rotation_matrix,unit_vectors(:,2))
      unit_vectors(:,3)=matmul(rotation_matrix,unit_vectors(:,3))
      
      ! perform rotation about y axis
      rotation_matrix=rotation_matrix_3d(unit_vectors(:,2),angles(2))
      ! rotate vertices
      do i=1,8
         vertices(:,i)=matmul(rotation_matrix,vertices(:,i))
      end do
      ! rotate unit vectors
      unit_vectors(:,1)=matmul(rotation_matrix,unit_vectors(:,1))
      unit_vectors(:,3)=matmul(rotation_matrix,unit_vectors(:,3))
      
      ! perform rotation about x axis
      rotation_matrix=rotation_matrix_3d(unit_vectors(:,3),angles(3))
      ! rotate vertices
      do i=1,8
         vertices(:,i)=matmul(rotation_matrix,vertices(:,i))
      end do

      return

   end function
   
      ! cube of edge length 1 rotated through euler angles 
   pure function rotated_hcp_ball_vertices(angles) result(vertices)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: angles(3)  ! three Euler angles

      ! result declaration
      real(kind=rel_kind) :: vertices(3,13)        ! eight three-dimensional positions
      
      ! variable declarations
      integer(kind=int_kind) :: i                  ! counter
      real(kind=rel_kind) :: unit_vectors(3,3)     ! three unit vectors
      real(kind=rel_kind) :: rotation_matrix(3,3)       ! rotation matrix
      
      ! set initial vector/vertices
      vertices(:,1)=(/0._rel_kind,sqrt(3._rel_kind)/3._rel_kind,-sqrt(6._rel_kind)/3._rel_kind/)
      vertices(:,2)=(/-0.5_rel_kind,-sqrt(3._rel_kind)/6._rel_kind,-sqrt(6._rel_kind)/3._rel_kind/)
      vertices(:,3)=(/0.5_rel_kind,-sqrt(3._rel_kind)/6._rel_kind,-sqrt(6._rel_kind)/3._rel_kind/)
      vertices(:,4)=(/-0._rel_kind,-sqrt(3._rel_kind)/3._rel_kind,sqrt(6._rel_kind)/3._rel_kind/)
      vertices(:,5)=(/0.5_rel_kind,sqrt(3._rel_kind)/6._rel_kind,sqrt(6._rel_kind)/3._rel_kind/)
      vertices(:,6)=(/-0.5_rel_kind,sqrt(3._rel_kind)/6._rel_kind,sqrt(6._rel_kind)/3._rel_kind/)
      vertices(:,7)=(/-0.5_rel_kind,0.5_rel_kind*sqrt(3._rel_kind),0._rel_kind/)
      vertices(:,8)=(/0.5_rel_kind,0.5_rel_kind*sqrt(3._rel_kind),0._rel_kind/)
      vertices(:,9)=(/0.5_rel_kind,-0.5_rel_kind*sqrt(3._rel_kind),0._rel_kind/)
      vertices(:,10)=(/-0.5_rel_kind,-0.5_rel_kind*sqrt(3._rel_kind),0._rel_kind/)
      vertices(:,11)=(/-1._rel_kind,0._rel_kind,0._rel_kind/)
      vertices(:,12)=(/1._rel_kind,0._rel_kind,0._rel_kind/)
      vertices(:,13)=0._rel_kind
      
      
      unit_vectors(:,1)=(/1._rel_kind,0._rel_kind,0._rel_kind/)
      unit_vectors(:,2)=(/0._rel_kind,1._rel_kind,0._rel_kind/)
      unit_vectors(:,3)=(/0._rel_kind,0._rel_kind,1._rel_kind/)
      
      ! perform rotation about x axis
      rotation_matrix=rotation_matrix_3d(unit_vectors(:,1),angles(1))
      ! rotate vertices
      do i=1,12
         vertices(:,i)=matmul(rotation_matrix,vertices(:,i))
      end do
      ! rotate unit vectors
      unit_vectors(:,2)=matmul(rotation_matrix,unit_vectors(:,2))
      unit_vectors(:,3)=matmul(rotation_matrix,unit_vectors(:,3))
      
      ! perform rotation about y axis
      rotation_matrix=rotation_matrix_3d(unit_vectors(:,2),angles(2))
      ! rotate vertices
      do i=1,12
         vertices(:,i)=matmul(rotation_matrix,vertices(:,i))
      end do
      ! rotate unit vectors
      unit_vectors(:,1)=matmul(rotation_matrix,unit_vectors(:,1))
      unit_vectors(:,3)=matmul(rotation_matrix,unit_vectors(:,3))
      
      ! perform rotation about x axis
      rotation_matrix=rotation_matrix_3d(unit_vectors(:,3),angles(3))
      ! rotate vertices
      do i=1,12
         vertices(:,i)=matmul(rotation_matrix,vertices(:,i))
      end do

      return

   end function
   
   
!-------------------------------------------------------------
!-------------------------------------------------------------

! DISTRIBUTION MOMENTS

!-------------------------------------------------------------
!-------------------------------------------------------------

   ! distribution moment
   pure function dist_moment(x,c,n,w_in) result(mu_n)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(:)                ! x values
      real(kind=rel_kind),intent(in) :: c                   ! centre
      integer(kind=int_kind),intent(in) :: n                ! order
      real(kind=rel_kind),intent(in),optional :: w_in(:)    ! weights
      
      ! result declaration
      real(kind=rel_kind) :: mu_n                           ! moment
      
      ! variable declarations
      real(kind=rel_kind) :: w(size(x))                     ! weights
      
      
      if (present(w_in)) then
         w=w_in
      else
         w=1._rel_kind
      end if
      
      mu_n=sum(w*(x-c)**n)/sum(w)
      
      return
   
   end function

   ! arithmetic mean
   pure function mean(x,w) result(mu)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(:)                ! x values
      real(kind=rel_kind),intent(in),optional :: w(:)       ! weights
      
      ! variable declarations
      real(kind=rel_kind) :: mu                             ! mean
      
      mu=dist_moment(x,0._rel_kind,1_int_kind,w)
      
      return
   
   end function
   
   ! central moment of distribution
   pure function central_moment(x,n,w) result(mu_n)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(:)                ! x values
      integer(kind=int_kind),intent(in) :: n                ! order
      real(kind=rel_kind),intent(in),optional :: w(:)       ! weights
      
      ! result declaration
      real(kind=rel_kind) :: mu_n                           ! moment
      
      ! variable declarations
      real(kind=rel_kind) :: mu                             ! mean
      
      mu=mean(x,w)
      mu_n=dist_moment(x,mu,n,w)
      
      return
   
   end function
   
   ! arithmetic standard deviation
   pure function standard_deviation(x,w) result(sigma)

      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(:)                ! x values
      real(kind=rel_kind),intent(in),optional :: w(:)       ! weights
      
      ! result declaration
      real(kind=rel_kind) :: sigma                          ! standard deviation
      
      sigma=sqrt(central_moment(x,2_int_kind,w))
   
      return
   
   end function
   
   ! arithmetic skewness
   pure function skewness(x,w) result(gamma)

      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(:)                ! x values
      real(kind=rel_kind),intent(in),optional :: w(:)       ! weights
      
      ! result declaration
      real(kind=rel_kind) :: gamma                          ! skewness
      
      ! variable declarations
      real(kind=rel_kind) :: mu                             ! mean
      
      mu=mean(x,w)
      gamma=dist_moment(x,mu,3_int_kind,w)/sqrt(dist_moment(x,mu,2_int_kind,w))**3
      
      return
      
   end function
   
   ! arithmetic kurtosis
   pure function kurtosis(x,w) result(beta)

      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(:)                ! x values
      real(kind=rel_kind),intent(in),optional :: w(:)       ! weights
      
      ! result declaration
      real(kind=rel_kind) :: beta                           ! kurtosis
      
      ! variable declarations
      real(kind=rel_kind) :: mu                             ! mean

      mu=mean(x,w)
      beta=dist_moment(x,mu,4_int_kind,w)/dist_moment(x,mu,2_int_kind,w)**2
      
      return
      
   end function
   
end module
