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

! blackbody external radiation field module
module m_source_external_bb

   use m_kind_parameters
   use m_constants_parameters
   use m_maths
   use m_source
   
   implicit none
   
   private
   public :: source_external_bb
   
   ! blackbody extended source type
   type,extends(source_external) :: source_external_bb
   
      contains
      
      procedure :: initialise
      
   end type
   
   contains
   
   ! set up diluted blackbody radiation field 
   subroutine initialise(self,position,velocity,&
      &temperature,luminosity,beta,wavelength_min,wavelength_max,n_wavelength,&
      &radius,dilution,gal_r,gal_z,g_0,add_cmb,mass,age,metallicity,luminosity_file)
   
      ! argument declarations
      class(source_external_bb),intent(inout) :: self             ! source object
      real(kind=rel_kind),intent(in) :: position(n_dim)           ! position of source
      real(kind=rel_kind),intent(in) :: velocity(n_dim)           ! velocity of source
      real(kind=rel_kind),intent(in),optional :: temperature      ! temperature of source (required)
      real(kind=rel_kind),intent(in),optional :: luminosity       ! luminosity of source
      real(kind=rel_kind),intent(in),optional :: beta             ! greybody exponent
      real(kind=rel_kind),intent(in),optional :: wavelength_min   ! minimum wavelength (required)
      real(kind=rel_kind),intent(in),optional :: wavelength_max   ! maximum wavelength (required)
      integer(kind=int_kind),intent(in),optional :: n_wavelength  ! number of wavelength points (required)
      real(kind=rel_kind),intent(in),optional :: radius           ! radius of source (required)
      real(kind=rel_kind),intent(in),optional :: dilution         ! blackbody dilution factor (required)
      real(kind=rel_kind),intent(in),optional :: gal_r            ! galactic radius
      real(kind=rel_kind),intent(in),optional :: gal_z            ! galactic height
      real(kind=rel_kind),intent(in),optional :: g_0              ! optional g_0 normalisation
      logical(kind=log_kind),intent(in),optional :: add_cmb       ! add cosmic microwave background (required)
      real(kind=rel_kind),intent(in),optional :: mass             ! mass of star
      real(kind=rel_kind),intent(in),optional :: age              ! age of star
      real(kind=rel_kind),intent(in),optional :: metallicity      ! metallicity of star
      character(kind=chr_kind,len=string_length),intent(in),optional :: luminosity_file   ! mono luminosity file name
      
      ! set variables
      self%position=position
      self%velocity=velocity
      self%radius=radius
      
      ! set wavelength array
      allocate(self%wavelength_array(n_wavelength))
      allocate(self%intensity_array(n_wavelength))
      allocate(self%norm_intensity_array(n_wavelength))
      allocate(self%cum_intensity_array(n_wavelength))
      self%wavelength_array=log_lin_space(wavelength_min,wavelength_max,n_wavelength)
      
      ! set intensity array
      self%intensity_array=dilution*planck_lambda(self%wavelength_array,temperature)
      
      ! add cosmic microwave background
      if (add_cmb) self%intensity_array=self%intensity_array+planck_lambda(self%wavelength_array,2.72548_rel_kind)
      
      ! calculate g_0
      self%g_0=trapz_intgr(pack(self%wavelength_array,self%wavelength_array<0.20664032_rel_kind),&
         &pack(self%intensity_array,self%wavelength_array<0.20664032_rel_kind))/1.2e-4_rel_kind
      
      ! calculate ave. intensity etc
      self%bolometric_intensity=trapz_intgr(self%wavelength_array,self%intensity_array)
      self%luminosity=self%bolometric_intensity*4._rel_kind*pi**2*radius**2
      
      ! normalise intensity and calculate cumulative distribution
      self%norm_intensity_array=self%intensity_array/self%bolometric_intensity
      self%cum_intensity_array=cum_dist_func(self%wavelength_array,self%intensity_array,.true._log_kind)
      
      
      return
   
   end subroutine
      
end module