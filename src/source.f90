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

! module defining point and external sources
module m_source

   use m_kind_parameters
   use m_constants_parameters
   use m_maths
   
   implicit none
   
   private
   public :: source
   public :: source_point
   public :: source_external
   
   ! abstract types
   
   ! radiation source type
   type,abstract :: source
   
      real(kind=rel_kind) :: luminosity                                 ! luminosity of source
      real(kind=rel_kind) :: position(n_dim)                            ! position of source
      real(kind=rel_kind) :: velocity(n_dim)                            ! velocity of source
      real(kind=rel_kind) :: radius                                     ! radius of closed surface
      real(kind=rel_kind) :: g_0                                        ! ionising intensity in Harbing units
      real(kind=rel_kind) :: bolometric_intensity                       ! intensity over all wavelength
      real(kind=rel_kind),allocatable :: wavelength_array(:)            ! array of wavelengths
      real(kind=rel_kind),allocatable :: intensity_array(:)             ! array of intensity
      real(kind=rel_kind),allocatable :: norm_intensity_array(:)        ! normalised intensity array
      real(kind=rel_kind),allocatable :: cum_intensity_array(:)         ! cumulative intensity (normalised)
      
      contains
      
      procedure(initialise_virtual),deferred :: initialise
      procedure,non_overridable :: destroy
      procedure,non_overridable :: intensity
      procedure,non_overridable :: random_wavelength
      procedure,non_overridable :: random_position
      procedure(random_direction_virtual),deferred :: random_direction
   
   end type
   
   ! point source type
   type,extends(source),abstract :: source_point
   
      contains
      
      procedure :: random_direction=>random_direction_source_point
   
   end type
   
   ! 3d external source (radiation field)
   type,extends(source),abstract :: source_external
   
      contains
      
      procedure :: random_direction=>random_direction_source_external
   
   end type
   
   ! abstract interfaces
   
   abstract interface
      subroutine initialise_virtual(self,position,velocity,&
         &temperature,luminosity,wavelength_min,wavelength_max,n_wavelength,&
         &radius,dilution,gal_r,gal_z,g_0,add_cmb,mass,age,metallicity)
         
         import :: source,int_kind,rel_kind,log_kind,n_dim
         
         ! argument declarations
         class(source),intent(inout) :: self                         ! source object
         real(kind=rel_kind),intent(in) :: position(n_dim)           ! position of source
         real(kind=rel_kind),intent(in) :: velocity(n_dim)           ! velocity of source
         real(kind=rel_kind),intent(in),optional :: temperature      ! temperature of source
         real(kind=rel_kind),intent(in),optional :: luminosity       ! luminosity of source
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
         
      end subroutine
   end interface
   
   abstract interface
      function random_direction_virtual(self,position) result (value)
   
         import :: source,rel_kind,n_dim
   
         ! argument declarations
         class(source),intent(in) :: self                            ! radiation source object
         real(kind=rel_kind),intent(in),optional :: position(n_dim)  ! position on sphere
   
         ! result declaration
         real(kind=rel_kind) :: value(n_dim)
   
      end function
   end interface
   
   contains
   
   elemental subroutine destroy(self)
   
      ! argument declarations
      class(source),intent(inout) :: self                               ! radiation source object
      
      deallocate(self%wavelength_array)
      deallocate(self%intensity_array)
      deallocate(self%norm_intensity_array)
      deallocate(self%cum_intensity_array)
      
      return
   
   end subroutine
   
   elemental function intensity(self,wavelength_ob,v_rec) result (value)
   
      ! argument declarations
      class(source),intent(in) :: self                                  ! radiation source object
      real(kind=rel_kind),intent(in) :: wavelength_ob                   ! wavelength in observer's frame
      real(kind=rel_kind),intent(in) :: v_rec                           ! recession velocity
      
      ! result declaration
      real(kind=rel_kind) :: value                                      ! intensity value (cgs-microns)
      
      ! variable declarations
      real(kind=rel_kind) :: wavelength_em                              ! wavelength in emitting frame
      
      wavelength_em=wavelength_ob*(v_rec/c_light+1._rel_kind)
      
      value=lookup_and_interpolate(wavelength_em,self%wavelength_array,self%intensity_array)
      
      return
   
   end function
   
   function random_wavelength(self) result (value)
   
      ! argument declarations
      class(source),intent(in) :: self                                  ! radiation source object
      
      ! result declaration
      real(kind=rel_kind) :: value                                      ! wavelength value (microns)
      
      ! variable declarations
      real(kind=rel_kind) :: r_num                                      ! random_number
      
      call random_number(r_num)
      
      value=lookup_and_inv_interpolate(r_num,self%cum_intensity_array,&
         &self%wavelength_array,self%norm_intensity_array)

      return
         
   end function
      
   function random_position(self) result (value)
   
      ! argument declarations
      class(source),intent(in) :: self              ! source object
   
      ! result declaration
      real(kind=rel_kind) :: value(n_dim)
      
      call random_direction(value)
      value=self%position+self%radius*value
      
      return
   
   end function
   
   function random_direction_source_point(self,position) result (value)
   
      ! argument declarations
      class(source_point),intent(in) :: self                    ! radiation source object
      real(kind=rel_kind),intent(in),optional :: position(n_dim)! position on sphere
   
      ! result declaration
      real(kind=rel_kind) :: value(n_dim)
      
      call random_direction(value)
      
      return
   
   end function
   
   function random_direction_source_external(self,position) result (value)
   
      ! argument declarations
      class(source_external),intent(in) :: self                 ! radiation source object
      real(kind=rel_kind),intent(in),optional :: position(n_dim)! position on sphere (required)
   
      ! result declaration
      real(kind=rel_kind) :: value(n_dim)
      
      ! variable declaration
      real(kind=rel_kind) :: cos_theta                          ! cosine between norm and value
      real(kind=rel_kind) :: surface_norm(n_dim)                ! surface normal vector
      
      call random_number(cos_theta)
      cos_theta=sqrt(cos_theta)
      surface_norm=(position-self%position)
      surface_norm=-surface_norm/norm2(surface_norm)
      
      call scattered_direction_3d(value,surface_norm,cos_theta)
      
      return
   
   end function
   
end module