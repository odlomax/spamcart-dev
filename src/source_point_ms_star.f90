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

! Main Sequence star point source module
module m_source_point_ms_star

   use m_kind_parameters
   use m_constants_parameters
   use m_maths
   use m_source
   use m_ms_star_sed
   
   implicit none
   
   private
   public :: source_point_ms_star
   
   ! blackbody point source
   type,extends(source_point) :: source_point_ms_star
   
      contains
      
      procedure :: initialise
   
   end type
   
   contains
      
   subroutine initialise(self,position,velocity,&
      &temperature,luminosity,beta,wavelength_min,wavelength_max,n_wavelength,&
      &radius,dilution,gal_r,gal_z,g_0,add_cmb,mass,age,metallicity,luminosity_file)
      
      ! argument declarations
      class(source_point_ms_star),intent(inout) :: self           ! source object
      real(kind=rel_kind),intent(in) :: position(n_dim)           ! position of source
      real(kind=rel_kind),intent(in) :: velocity(n_dim)           ! velocity of source
      real(kind=rel_kind),intent(in),optional :: temperature      ! temperature of source
      real(kind=rel_kind),intent(in),optional :: luminosity       ! luminosity of source
      real(kind=rel_kind),intent(in),optional :: beta             ! greybody exponent
      real(kind=rel_kind),intent(in),optional :: wavelength_min   ! minimum wavelength
      real(kind=rel_kind),intent(in),optional :: wavelength_max   ! maximum wavelength
      integer(kind=int_kind),intent(in),optional :: n_wavelength  ! number of wavelength points
      real(kind=rel_kind),intent(in),optional :: radius           ! radius of source
      real(kind=rel_kind),intent(in),optional :: dilution         ! blackbody dilution factor
      real(kind=rel_kind),intent(in),optional :: gal_r            ! galactic radius (kpc)
      real(kind=rel_kind),intent(in),optional :: gal_z            ! galactic height (kpc)
      real(kind=rel_kind),intent(in),optional :: g_0              ! optional g_0 normalisation
      logical(kind=log_kind),intent(in),optional :: add_cmb       ! add cosmic microwave background
      real(kind=rel_kind),intent(in),optional :: mass             ! mass of star
      real(kind=rel_kind),intent(in),optional :: age              ! age of star
      real(kind=rel_kind),intent(in),optional :: metallicity      ! metallicity of star
      character(kind=chr_kind,len=string_length),intent(in),optional :: luminosity_file   ! mono luminosity file name
      
      ! variable declarations
      type(ms_star_sed) :: star_sed                               ! star sed object
      
      ! initialise star sed object
      call star_sed%initialise(mass,age,metallicity)
      
      ! allocate arrays
      allocate(self%wavelength_array(size(star_sed%wavelength_array)))
      allocate(self%intensity_array(size(star_sed%wavelength_array)))
      allocate(self%norm_intensity_array(size(star_sed%wavelength_array)))
      allocate(self%cum_intensity_array(size(star_sed%wavelength_array)))
      
      self%position=position
      self%velocity=velocity
      
      ! radius and g_0 undefined for point source
      self%radius=0._rel_kind
      self%g_0=0._rel_kind
      
      ! populate arrays
      self%wavelength_array=star_sed%wavelength_array
      self%intensity_array=star_sed%intensity_array
      
      ! calculate luminosity
      self%luminosity=star_sed%bol_luminosity
      self%bolometric_intensity=trapz_intgr(self%wavelength_array,self%intensity_array)
      
      ! normalise intensity and calculate cumulative distribution
      self%norm_intensity_array=self%intensity_array/&
         &trapz_intgr(self%wavelength_array,self%intensity_array)
      self%cum_intensity_array=cum_dist_func(self%wavelength_array,self%intensity_array,.true._log_kind)
      
      ! destroy star sed object
      call star_sed%destroy()
      
      return
      
   end subroutine
   
end module