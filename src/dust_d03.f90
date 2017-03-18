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

! Draine 2003 Dust module.
module m_dust_d03 

   use m_kind_parameters
   use m_constants_parameters
   use m_maths
   use m_dust
   
   implicit none
   
   ! entities private by default
   private
   
   ! public entities
   public :: dust_d03
   
   ! define Draine (2003) dust class
   type,extends(dust) :: dust_d03
   
      contains
      
      procedure :: initialise
   
   end type
   
   contains
   
   ! initialise Draine 2003 dust class
   subroutine initialise(self,r_v,temperature_min,temperature_max,n_temperature,gas_to_dust_mass,iso_scatter)
   
      ! argument declarations
      class(dust_d03),intent(inout) :: self                             ! dust object
      integer(kind=int_kind),intent(in) :: n_temperature                ! number of temperatures
      real(kind=rel_kind),intent(in) :: temperature_min                 ! minimum temperature
      real(kind=rel_kind),intent(in) :: temperature_max                 ! maximum temperature
      real(kind=rel_kind),intent(in) :: r_v                             ! ratio of absolute to selective extinction
      real(kind=rel_kind),intent(in),optional :: gas_to_dust_mass       ! M_dust/M_gas
      logical(kind=log_kind),intent(in),optional :: iso_scatter         ! is dust scattering isotropic?
      
      ! variable declarations
      integer(kind=int_kind) :: i                                       ! counter
      integer(kind=int_kind) :: j                                       ! counter
      integer(kind=int_kind) :: i_gas_to_dust_mass                      ! line with dust to gas mass ratio in file
      integer(kind=int_kind) :: i_start                                 ! starting line in file
      integer(kind=int_kind) :: n_wavelength                            ! number of wavelengths
      real(kind=rel_kind) :: dummy                                      ! dummy variable
      real(kind=rel_kind) :: model_gas_to_dust_mass                     ! dust to gas mass ratio from model
      real(kind=rel_kind) :: temp_blackbody                             ! integrated blackbody spectrum
      real(kind=rel_kind),allocatable :: temp_dust_mass_abs(:,:)        ! temp dust mass absorption coeff
      real(kind=rel_kind),allocatable :: temp_albedo(:,:)               ! temp albedo
      real(kind=rel_kind),allocatable :: temp_mean_cos_theta(:,:)       ! temp mean cosine
      real(kind=rel_kind),allocatable :: temp_gas_to_dust_mass(:)       ! temp dust to gas mass ratio
      real(kind=rel_kind),allocatable :: temp_blackbody_array(:)        ! temporary blackbody spectrum array
      real(kind=rel_kind),allocatable :: iso_dust_mass_ext_array(:)     ! isotropic scattering dust mass ext
      real(kind=rel_kind),allocatable :: iso_albedo_array(:)            ! isotropic albedo array
      real(kind=rel_kind),allocatable :: r_v_value(:)                   ! r_v value array
      character(kind=chr_kind,len=3) :: r_v_string(3)                   ! r_v string array
      
      ! set parameters
      n_wavelength=1064
      i_gas_to_dust_mass=76
      i_start=81
      r_v_string=(/"3.1","4.0","5.5"/)
      
      ! allocate arrays
      allocate(self%temperature_array(n_temperature))
      allocate(self%wavelength_array(n_wavelength))
      allocate(self%dust_mass_ext_array(n_wavelength))
      allocate(self%albedo_array(n_wavelength))
      allocate(self%mean_cos_theta_array(n_wavelength))
      allocate(self%mono_mass_emissivity_array(n_wavelength,n_temperature))
      allocate(self%cum_mono_mass_emissivity(n_wavelength,n_temperature))
      allocate(self%norm_mono_mass_emissivity(n_wavelength,n_temperature))
      allocate(self%bol_mass_emissivity_array(n_temperature))
      allocate(self%planck_albedo_array(n_temperature))
      allocate(r_v_value(size(r_v_string)))
      allocate(temp_dust_mass_abs(size(r_v_string),n_wavelength))
      allocate(temp_albedo(size(r_v_string),n_wavelength))
      allocate(temp_mean_cos_theta(size(r_v_string),n_wavelength))
      allocate(temp_gas_to_dust_mass(size(r_v_string)))
      allocate(temp_blackbody_array(n_wavelength))
      
      ! get temperatures
      self%temperature_array=log_lin_space(temperature_min,temperature_max,n_temperature)
      
      ! get wavelengths
      open(1,file="./dustdata/draine_rv3.1.dat")
      
      ! get to end of file
      do i=1,i_start+n_wavelength-1
         read(1,*)
      end do
      backspace(1)
      
      ! read in wavelengths backwards
      do i=1,n_wavelength
         read(1,*) self%wavelength_array(i)
         backspace(1)
         backspace(1)
      end do
      
      close(1)
      
      ! get r_v values
      do i=1,size(r_v_string)
         read(r_v_string(i),*) r_v_value(i)
      end do
      
      ! populate temp arrays
      do i=1,size(r_v_string)
      
         open(1,file="./dustdata/draine_rv"//r_v_string(i)//".dat")
         
         ! get temp dust to mass ratio
         do j=1,i_gas_to_dust_mass-1
            read(1,*)
         end do
         read(1,*) temp_gas_to_dust_mass(i)
         rewind(1)
         
         ! get to end of file
         do j=1,i_start+n_wavelength-1
            read(1,*)
         end do
         backspace(1)
         
         ! read in values backwards
         do j=1,n_wavelength
            read(1,*) dummy,temp_albedo(i,j),temp_mean_cos_theta(i,j),dummy,temp_dust_mass_abs(i,j)
            backspace(1)
            backspace(1)
         end do
         
         close(1)
      
      end do
      
      ! interpolate tabulated values
      do i=1,n_wavelength
      
         self%albedo_array(i)=lookup_and_geo_interpolate(r_v,r_v_value,temp_albedo(:,i))
         self%mean_cos_theta_array(i)=lookup_and_interpolate(r_v,r_v_value,temp_mean_cos_theta(:,i))
         self%dust_mass_ext_array(i)=&
            &lookup_and_geo_interpolate(r_v,r_v_value,temp_dust_mass_abs(:,i))/&
            &(1._rel_kind-self%albedo_array(i))
      
      end do
      model_gas_to_dust_mass=lookup_and_geo_interpolate(r_v,r_v_value,temp_gas_to_dust_mass)

      ! normalise mass extinction array to dust-to-gas-mass ratio
      ! override model value if gas_to_dust_mass argument is present
      if (present(gas_to_dust_mass)) model_gas_to_dust_mass=gas_to_dust_mass
      self%dust_mass_ext_array=&
         &self%dust_mass_ext_array/(model_gas_to_dust_mass+1._rel_kind)
         
      ! deallocate temp arrays
      deallocate(temp_albedo)
      deallocate(temp_mean_cos_theta)
      deallocate(temp_dust_mass_abs)
      deallocate(temp_gas_to_dust_mass)
      deallocate(r_v_value)
      
      ! modify values if using isotropic scattering
      ! phase function is sum of isotropic and delta function components
      if (present(iso_scatter)) then
         if (iso_scatter) then
         
            ! allocate arrays
            allocate(iso_dust_mass_ext_array(n_wavelength))
            allocate(iso_albedo_array(n_wavelength))
            
            ! remove delta function component from scattering opacity
            iso_dust_mass_ext_array=&
               &(1._rel_kind-self%albedo_array)*self%dust_mass_ext_array+&
               &self%albedo_array*self%dust_mass_ext_array*(1._rel_kind-self%mean_cos_theta_array)
            
            iso_albedo_array=(self%albedo_array-self%mean_cos_theta_array*self%albedo_array)/&
               &(1._rel_kind-self%mean_cos_theta_array*self%albedo_array)
               
            ! set isotropic phase function
            self%mean_cos_theta_array=0._rel_kind
               
            call move_alloc(iso_dust_mass_ext_array,self%dust_mass_ext_array)
            call move_alloc(iso_albedo_array,self%albedo_array)
         
         end if
      end if
      
      ! calculate temperature dependent values
      do i=1,n_temperature
      
         temp_blackbody_array=planck_lambda(self%wavelength_array,self%temperature_array(i))
         temp_blackbody=trapz_intgr(self%wavelength_array,temp_blackbody_array)
      
         self%mono_mass_emissivity_array(:,i)=temp_blackbody_array*self%dust_mass_ext_array*(1._rel_kind-self%albedo_array)
         self%cum_mono_mass_emissivity(:,i)=&
            &cum_dist_func(self%wavelength_array,self%mono_mass_emissivity_array(:,i),.true._log_kind)
         self%bol_mass_emissivity_array(i)=trapz_intgr(self%wavelength_array,self%mono_mass_emissivity_array(:,i))
         self%norm_mono_mass_emissivity(:,i)=self%mono_mass_emissivity_array(:,i)/&
            &self%bol_mass_emissivity_array(i)
            
         do j=1,n_wavelength
            self%planck_albedo_array(i)=trapz_intgr(self%wavelength_array,temp_blackbody_array*self%albedo_array)/temp_blackbody
         end do   
      
      end do
      
      deallocate(temp_blackbody_array)
      
      return
   
   end subroutine
   
end module