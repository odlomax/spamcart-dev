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

! custom point source module
module m_source_point_custom

   use m_kind_parameters
   use m_constants_parameters
   use m_maths
   use m_source
   
   implicit none
   
   private
   public :: source_point_custom
   
   ! custom point source
   type,extends(source_point) :: source_point_custom
   
      contains
      
      procedure :: initialise
   
   end type
   
   contains
      
   subroutine initialise(self,position,velocity,&
      &temperature,luminosity,beta,wavelength_min,wavelength_max,n_wavelength,&
      &radius,dilution,gal_r,gal_z,g_0,add_cmb,mass,age,metallicity,luminosity_file)
      
      ! argument declarations
      class(source_point_custom),intent(inout) :: self            ! source object
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
      integer(kind=int_kind) :: i                                 ! counter
      integer(kind=int_kind) :: n_line                            ! number of lines in file
      integer(kind=int_kind) :: read_status                       ! read status of file
      
      ! set position
      self%position=position
      self%velocity=velocity
      
      ! radius and g_0 undefined for point source
      self%radius=0._rel_kind
      self%g_0=0._rel_kind
      
      
      ! read size of luminosity file
      open(1,file=trim(luminosity_file))
      n_line=0
      do
      
         read(1,*,iostat=read_status)
         
         if (read_status<0) exit
         
         n_line=n_line+1
         
      end do
      ! allocate array
      allocate(self%wavelength_array(n_line))
      allocate(self%intensity_array(n_line))
      allocate(self%norm_intensity_array(n_line))
      allocate(self%cum_intensity_array(n_line))
      
      ! read luminosity file
      rewind(1)
      do i=1,n_line
         read(1,*) self%wavelength_array(i), self%intensity_array(i)
      end do
      close(1)
      
      ! calculate luminosity
      self%luminosity=trapz_intgr(self%wavelength_array,self%intensity_array)
      self%bolometric_intensity=self%luminosity
      
      ! normalise intensity and calculate cumulative distribution
      self%norm_intensity_array=self%intensity_array/&
         &trapz_intgr(self%wavelength_array,self%intensity_array)
      self%cum_intensity_array=cum_dist_func(self%wavelength_array,self%intensity_array,.true._log_kind)
      
      return
      
   end subroutine
   
end module