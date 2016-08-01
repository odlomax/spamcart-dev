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

! module containing compile-time parameters
module m_kind_parameters

   use iso_c_binding
   
   implicit none 
   
   ! variable types
   integer,parameter :: rel_kind=c_double               ! c double/float kind
   integer,parameter :: int_kind=c_int                  ! c int kind
   integer,parameter :: cpx_kind=c_double_complex       ! c double complex kind
   integer,parameter :: chr_kind=c_char                 ! c character kind
   integer,parameter :: log_kind=c_bool                 ! c boolean kind
   integer,parameter :: string_length=2048              ! default character length
   
   integer(kind=int_kind),parameter :: n_dim=3          ! number of dimensions
   
end module

module m_sph_parameters

   use m_kind_parameters
   
   implicit none

   ! sph parameters
   integer(kind=int_kind),parameter :: n_leaf=2**n_dim          ! max particles per leaf cell
   integer(kind=int_kind),parameter :: h_iterations=10          ! maximum number of h iterations
   integer(kind=int_kind),parameter :: n_lookup=64              ! number of entries in kernel lookup table
   real(kind=rel_kind),parameter :: delta_h=0.01_rel_kind       ! smoothing length convergence parameter
   
   ! ray parameters
   integer(kind=int_kind),parameter :: n_ray=128                ! initial size of ray_particle array (can be resized)
   real(kind=rel_kind),parameter :: ray_eta=1.2_rel_kind        ! ray length over estimation factor
   real(kind=rel_kind),parameter :: ray_conv=0.01_rel_kind      ! ray convergence parameter
   
end module

module m_constants_parameters

   use m_kind_parameters
   
   implicit none
   
   ! numbers
   real(kind=rel_kind),parameter :: pi=2._rel_kind*acos(0._rel_kind)    ! 3.1412...
   real(kind=rel_kind),parameter :: one_third=1._rel_kind/3._rel_kind   ! 1/3
   real(kind=rel_kind),parameter :: arcsec_rad=pi/6.48e+5_rel_kind     ! 1 arcsec in radians
   real(kind=rel_kind),parameter :: fwhm_sigma=2._rel_kind*sqrt(2._rel_kind*log(2._rel_kind))   ! FWHM in units of sigma
   
   ! physical constants
   real(kind=rel_kind),parameter :: k_boltz=1.3806488e-23_rel_kind      ! Boltzmann constant [J/K]
   real(kind=rel_kind),parameter :: c_light=2.99792458e+8_rel_kind      ! speed of light [m/s]
   real(kind=rel_kind),parameter :: h_planck=6.6260696e-34_rel_kind     ! Planck constant [J/Hz]

end module