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

! Porter and Strong 2005 external radiation field module
module m_source_external_ps05

   use m_kind_parameters
   use m_constants_parameters
   use m_maths
   use m_source
   
   implicit none
   
   private
   public :: source_external_ps05
   
   ! porter and strong 2005 radiation field
   type,extends(source_external) :: source_external_ps05
   
      contains
      
      procedure :: initialise
      
   end type
   
   contains

   ! set up Porter and Strong 2005 interstellar radiation field
   subroutine initialise(self,position,velocity,&
      &temperature,luminosity,wavelength_min,wavelength_max,n_wavelength,&
      &radius,dilution,gal_r,gal_z,g_0,add_cmb,mass,age,metallicity)
   
      ! argument declarations
      class(source_external_ps05),intent(inout) :: self           ! source object
      real(kind=rel_kind),intent(in) :: position(n_dim)           ! position of source
      real(kind=rel_kind),intent(in) :: velocity(n_dim)           ! velocity of source
      real(kind=rel_kind),intent(in),optional :: temperature      ! temperature of source
      real(kind=rel_kind),intent(in),optional :: luminosity       ! luminosity of source
      real(kind=rel_kind),intent(in),optional :: wavelength_min   ! minimum wavelength
      real(kind=rel_kind),intent(in),optional :: wavelength_max   ! maximum wavelength
      integer(kind=int_kind),intent(in),optional :: n_wavelength  ! number of wavelength points
      real(kind=rel_kind),intent(in),optional :: radius           ! radius of source (required)
      real(kind=rel_kind),intent(in),optional :: dilution         ! blackbody dilution factor
      real(kind=rel_kind),intent(in),optional :: gal_r            ! galactic radius (required)
      real(kind=rel_kind),intent(in),optional :: gal_z            ! galactic height (required)
      real(kind=rel_kind),intent(in),optional :: g_0              ! optional g_0 normalisation (applicable)
      logical(kind=log_kind),intent(in),optional :: add_cmb       ! add cosmic microwave background (required)
      real(kind=rel_kind),intent(in),optional :: mass             ! mass of star
      real(kind=rel_kind),intent(in),optional :: age              ! age of star
      real(kind=rel_kind),intent(in),optional :: metallicity      ! metallicity of star
      
      ! variable declarations
      integer(kind=int_kind) :: i                                       ! counter
      integer(kind=int_kind) :: j                                       ! counter
      integer(kind=int_kind) :: k                                       ! counter
      integer(kind=int_kind) :: n_wavelength_file                       ! number of wavelengths
      real(kind=rel_kind) :: dummy                                      ! dummy variable
      real(kind=rel_kind),allocatable :: gal_r_array(:)                 ! galactic radius array
      real(kind=rel_kind),allocatable :: gal_z_array(:)                 ! galactic height array
      real(kind=rel_kind),allocatable :: energy_density_array(:,:,:)    ! energy density from file
      character(kind=chr_kind,len=3) :: r_string(9)                     ! file name substrings
      character(kind=chr_kind,len=3) :: z_string(10)                    ! file name substrings

      
      ! set parameters
      n_wavelength_file=505
      r_string=(/"0  ","2  ","4  ","6  ","8.5","12 ","16 ","20 ","30 "/)
      z_string=(/"0  ","0.1","0.2","0.5","1  ","2  ","5  ","10 ","20 ","30 "/)
      
            
      ! set variables
      self%position=position
      self%velocity=velocity
      self%radius=radius
         
      ! allocate arrays
      allocate(self%wavelength_array(n_wavelength_file))
      allocate(self%intensity_array(n_wavelength_file))
      allocate(self%norm_intensity_array(n_wavelength_file))
      allocate(self%cum_intensity_array(n_wavelength_file))
      allocate(gal_r_array(size(r_string)))
      allocate(gal_z_array(size(z_string)))
      
      ! read in wavelength array
      open(1,file="./isrfdata/Standard_0_0_0_Flux.dat")
      do i=1,n_wavelength_file
         read(1,*) self%wavelength_array(i)
      end do
      close(1)
      
      ! read in gal_r_array
      do i=1,size(r_string)
         read(r_string(i),*) gal_r_array(i)
      end do
      
      ! read in gal_z_array
      do i=1,size(z_string)
         read(z_string(i),*) gal_z_array(i)
      end do
      
      ! allocate flux array
      allocate(energy_density_array(size(r_string),size(z_string),n_wavelength_file))
      
      ! read in fluxes (funny order, so that interpolation uses contiguous array slice)
      do j=1,size(z_string)
         do i=1,size(r_string)
         
            open(1,file="./isrfdata/Standard_"//trim(r_string(i))//"_0_"//trim(z_string(j))//"_Flux.dat")
            do k=1,n_wavelength_file
               read(1,*) dummy,energy_density_array(i,j,k)
            end do
            close(1)
         
         end do
      end do
      
      ! build intensity array
      do i=1,n_wavelength_file
      
         ! interpolate energy density (lambda eV cm^-3) and convert to intensity (erg s^-1 sr^-1 cm^-2 micron^-1)
         self%intensity_array(i)=&
            &lookup_and_geo_interpolate(gal_r,gal_z,gal_r_array,gal_z_array,energy_density_array(:,:,i))*&
            &0.0038222687_rel_kind/self%wavelength_array(i)
      
      end do
      
      ! calculate g_0
      self%g_0=trapz_intgr(pack(self%wavelength_array,self%wavelength_array<0.20664032_rel_kind),&
         &pack(self%intensity_array,self%wavelength_array<0.20664032_rel_kind))/1.2e-4_rel_kind
         
      ! override g_0, if argument present
      if (present(g_0)) then
      
         self%intensity_array=self%intensity_array*g_0/self%g_0
         self%g_0=g_0
      
      end if
      
      ! add cosmic microwave background
      if (add_cmb) self%intensity_array=self%intensity_array+planck_lambda(self%wavelength_array,2.72548_rel_kind)
      
      ! calculate ave. intensity etc
      self%bolometric_intensity=trapz_intgr(self%wavelength_array,self%intensity_array)
      self%luminosity=self%bolometric_intensity*4._rel_kind*pi**2*radius**2
      
      ! normalise intensity and calculate cumulative distribution
      self%norm_intensity_array=self%intensity_array/self%bolometric_intensity
      self%cum_intensity_array=cum_dist_func(self%wavelength_array,self%intensity_array,.true._log_kind)
   
      return
   
   end subroutine
   
end module