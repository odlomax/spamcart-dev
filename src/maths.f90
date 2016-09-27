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
   
   interface make_hist
      procedure :: make_hist_1d
   end interface
   
   interface rotate_vector
      procedure :: rotate_vector_3d
   end interface

   contains

   
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
      
      lambda_si=1e-6_rel_kind*lambda
      exp_part=exp(h_planck*c_light/(lambda_si*k_boltz*t))
   
      if (exp_part<huge(0._rel_kind)) then
         ! calculate planck function
         b=2._rel_kind*h_planck*c_light**2/&
            (lambda_si**5*(exp_part-1._rel_kind))
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
   pure function make_hist_1d(values,bins) result (hist)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: values(:)             ! array of values
      real(kind=rel_kind),intent(in) :: bins(:)               ! array of bin values
      
      ! result declaration
      real(kind=rel_kind) :: hist(size(bins)-1)               ! histogram values
      
      ! variable declarations
      integer(kind=int_kind) :: i                             ! counter
      integer(kind=int_kind) :: j
      
      ! initialise histogram
      hist=0._rel_kind
      
      ! loop over values and update histogram
      do i=1,size(values)
      
         j=binary_search(values(i),bins)
         hist(j)=hist(j)+1._rel_kind
      
      end do
      
      ! normalise bins
      forall (i=1:size(hist))
         hist(i)=hist(i)/(real(size(values),rel_kind)*(bins(i+1)-bins(i)))
      end forall
   
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
      
      forall (i=1:n_x)
         x(i)=x_min+real(i-1,rel_kind)*dx
      end forall
   
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
      
      forall (i=1:n_x)
         x(i)=exp(log_x_min+real(i-1,rel_kind)*dlog_x)
      end forall
   
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
            &(f_values(i-1)+f_values(i))*&
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
      do i=1,size(indices)
         indices(i)=i
      end do
      
      call quicksort(temp_array,indices)
      
      return
   
   end function
   
   ! quick sort alorithm
   pure recursive subroutine quicksort_1(array)

      !argument declarations
      real(kind=rel_kind),intent(inout) :: array(:)    !array to be sorted
      
      !variable declarations
      integer(kind=int_kind) :: i                      !counter
      integer(kind=int_kind) :: n                      !size of array
      integer(kind=int_kind) :: pivot                  !pivot index
      real(kind=rel_kind) :: temp                      !temporary element
      
      !Quicksort algorithm
      
      !return if array has one or zero elements
      n=size(array)
      if(n<=1) return

      !set pivot to central array element
      pivot=1+(n-1)/2
      
      !move pivot to end of array
      temp=array(pivot)
      array(pivot)=array(n)
      array(n)=temp
      
      !loop over array and place elements
      !lower than the pivot value left of pivot
      pivot=1
      do i=1,n-1
      
         if(array(i)<=array(n)) then
            temp=array(pivot)
            array(pivot)=array(i)
            array(i)=temp
            pivot=pivot+1
         end if
      
      end do
      
      !move pivot back from end of array
      temp=array(pivot)
      array(pivot)=array(n)
      array(n)=temp
      
      !sort the array slices either side of pivot
      call quicksort(array(1:pivot-1))
      call quicksort(array(pivot+1:n))
         
      return

   end subroutine
   
   ! quick sort alorithm
   ! also maintains list of original indices
   pure recursive subroutine quicksort_2(array,indices)

      !argument declarations
      integer(kind=int_kind),intent(inout) :: indices(:)  !original indices of array elements
      real(kind=rel_kind),intent(inout) :: array(:)    !array to be sorted
      
      !variable declarations
      integer(kind=int_kind) :: i                      !counter
      integer(kind=int_kind) :: n                      !size of array
      integer(kind=int_kind) :: pivot                  !pivot index
      integer(kind=int_kind) :: temp_i                 !temporary index
      real(kind=rel_kind) :: temp_r                    !temporary element
      
      
      !Quicksort algorithm
      
      !return if array has one or zero elements
      n=size(array)
      if(n<=1) return

      !set pivot to central array element
      pivot=1+(n-1)/2
      
      !move pivot to end of array
      temp_r=array(pivot)
      array(pivot)=array(n)
      array(n)=temp_r
      temp_i=indices(pivot)
      indices(pivot)=indices(n)
      indices(n)=temp_i
      
      !loop over array and place elements
      !lower than the pivot value left of pivot
      pivot=1
      do i=1,n-1
      
         if(array(i)<=array(n)) then
            temp_r=array(pivot)
            array(pivot)=array(i)
            array(i)=temp_r
            temp_i=indices(pivot)
            indices(pivot)=indices(i)
            indices(i)=temp_i
            pivot=pivot+1
         end if
      
      end do
      
      !move pivot back from end of array
      temp_r=array(pivot)
      array(pivot)=array(n)
      array(n)=temp_r
      temp_i=indices(pivot)
      indices(pivot)=indices(n)
      indices(n)=temp_i
      
      !sort the array slices either side of pivot
      call quicksort(array(1:pivot-1),indices(1:pivot-1))
      call quicksort(array(pivot+1:n),indices(pivot+1:n))
         
      return

   end subroutine
   
   pure function rotate_vector_3d(x,u,theta) result(x_new)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(3)       ! vector to rotate
      real(kind=rel_kind),intent(in) :: u(3)       ! axis of rotation unit vector
      real(kind=rel_kind),intent(in) :: theta      ! angle of rotation
      
      ! result declaration
      real(kind=rel_kind) :: x_new(3)              ! rotated x vector
      
      ! variable declarations
      real(kind=rel_kind) :: cos_theta             ! cos(theta)
      real(kind=rel_kind) :: sin_theta             ! sin(theta)
      
      cos_theta=cos(theta)
      sin_theta=sin(theta)
      
      x_new=matmul(reshape((/cos_theta+u(1)**2*(1._rel_kind-cos_theta),&
         &u(2)*u(1)*(1._rel_kind-cos_theta)+u(3)*sin_theta,&
         &u(3)*u(1)*(1._rel_kind-cos_theta)-u(2)*sin_theta,&
         &u(1)*u(2)*(1._rel_kind-cos_theta)-u(3)*sin_theta,&
         &cos_theta+u(2)**2*(1._rel_kind-cos_theta),&
         &u(3)*u(2)*(1._rel_kind-cos_theta)+u(1)*sin_theta,&
         &u(1)*u(3)*(1._rel_kind-cos_theta)+u(2)*sin_theta,&
         &u(2)*u(3)*(1._rel_kind-cos_theta)-u(1)*sin_theta,&
         &cos_theta+u(3)**2*(1._rel_kind-cos_theta)/),(/3,3/)),x)
      
      return   
   
   end function
   
end module
