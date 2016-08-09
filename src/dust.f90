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

! Dust module.
! Contains dust lookup tables
module m_dust 

   use m_kind_parameters
   use m_constants_parameters
   use m_maths
   
   implicit none
   
   ! entities private by default
   private
   
   ! public entities
   public :: dust
   
   ! define dust class
   type,abstract :: dust
   
      ! standard dust variables
   
      real(kind=rel_kind),allocatable :: temperature_array(:)              ! array of temperatures
      real(kind=rel_kind),allocatable :: wavelength_array(:)               ! array of wavelength
      real(kind=rel_kind),allocatable :: dust_mass_ext_array(:)            ! dust mass extinction as function of lambda
      real(kind=rel_kind),allocatable :: albedo_array(:)                   ! albedo as function of lambda
      real(kind=rel_kind),allocatable :: mean_cos_theta_array(:)           ! mean scattering cosine
      real(kind=rel_kind),allocatable :: mono_mass_emissivity_array(:,:)   ! monochromatic luminosity per unit dust mass
      real(kind=rel_kind),allocatable :: cum_mono_mass_emissivity(:,:)     ! cumulative mono_mass_emissivity along lambda
      real(kind=rel_kind),allocatable :: norm_mono_mass_emissivity(:,:)    ! normalised mono_mass_emissivity along lambda
      real(kind=rel_kind),allocatable :: bol_mass_emissivity_array(:)      ! luminosity per unit dust mass
      
      ! Modified Random Walk variables
      
      real(kind=rel_kind),allocatable :: mrw_zeta_array(:)                 ! mrw zeta parameter
      real(kind=rel_kind),allocatable :: mrw_y_array(:)                    ! mrw y parameter
      real(kind=rel_kind),allocatable :: mrw_planck_abs_array(:)           ! planck absorption coefficient
      real(kind=rel_kind),allocatable :: mrw_inv_planck_ext_array(:)       ! inverse planck extinction coefficient
      real(kind=rel_kind),allocatable :: mrw_planck_sca_array(:,:)         ! planck scattering coefficient (definite integral)
      real(kind=rel_kind),allocatable :: mrw_planck_function_array(:,:)    ! planck function
      real(kind=rel_kind),allocatable :: mrw_bol_planck_function_array(:)  ! integrated planck function
      real(kind=rel_kind),allocatable :: mrw_cum_planck_function(:,:)      ! cumulative planck function
      real(kind=rel_kind),allocatable :: mrw_norm_planck_function(:,:)     ! normalised planck function
            
      contains
      
      ! standard dust procedures
      
      procedure(initialise_virtual),deferred :: initialise
      procedure,non_overridable :: destroy
      procedure,non_overridable :: dust_mass_ext
      procedure,non_overridable :: albedo
      procedure,non_overridable :: dust_mass_ext_and_albedo
      procedure,non_overridable :: mean_scatter_cosine
      procedure,non_overridable :: mono_mass_emissivity
      procedure,non_overridable :: bol_mass_emissivity
      procedure,non_overridable :: phase_func
      procedure,non_overridable :: dust_temperature
      procedure,non_overridable :: random_wavelength
      procedure,non_overridable :: random_scatter_cosine
      procedure,non_overridable,nopass :: random_emission_direction
      procedure,non_overridable :: random_scatter_direction
      
      ! Modified Random Walk procedures
      
      procedure,non_overridable :: mrw_initialise
      procedure,non_overridable :: mrw_random_wavelength
      procedure,non_overridable :: mrw_planck_abs
      procedure,non_overridable :: mrw_inv_planck_ext
      procedure,non_overridable :: mrw_planck_sca
      procedure,non_overridable :: mrw_distance
      procedure,non_overridable,nopass :: mrw_random_direction
   
   end type

   abstract interface
      subroutine initialise_virtual(self,r_v,temperature_min,temperature_max,n_temperature,gas_to_dust_mass,iso_scatter)
      
      import :: dust,int_kind,rel_kind,log_kind
      
      ! argument declarations
      class(dust),intent(inout) :: self                                 ! dust object
      integer(kind=int_kind),intent(in) :: n_temperature                ! number of temperatures
      real(kind=rel_kind),intent(in) :: temperature_min                 ! minimum temperature
      real(kind=rel_kind),intent(in) :: temperature_max                 ! maximum temperature
      real(kind=rel_kind),intent(in) :: r_v                             ! ratio of absolute to selective extinction
      real(kind=rel_kind),intent(in),optional :: gas_to_dust_mass       ! M_dust/M_gas
      logical(kind=log_kind),intent(in),optional :: iso_scatter         ! is dust scattering isotropic?
      
      end subroutine
   end interface
   
   contains
   
   ! deallocate all arrays
   pure subroutine destroy(self)
   
      ! argument declarations
      class(dust),intent(inout) :: self                       ! dust object
      
      deallocate(self%temperature_array)
      deallocate(self%wavelength_array)
      deallocate(self%dust_mass_ext_array)
      deallocate(self%albedo_array)
      deallocate(self%mean_cos_theta_array)
      deallocate(self%mono_mass_emissivity_array)
      deallocate(self%cum_mono_mass_emissivity)
      deallocate(self%norm_mono_mass_emissivity)
      deallocate(self%bol_mass_emissivity_array)
      if (allocated(self%mrw_zeta_array)) deallocate(self%mrw_zeta_array)
      if (allocated(self%mrw_y_array)) deallocate(self%mrw_y_array)
      if (allocated(self%mrw_planck_abs_array)) deallocate(self%mrw_planck_abs_array)
      if (allocated(self%mrw_inv_planck_ext_array)) deallocate(self%mrw_inv_planck_ext_array)
      if (allocated(self%mrw_planck_sca_array)) deallocate(self%mrw_planck_sca_array)
      if (allocated(self%mrw_planck_function_array)) deallocate(self%mrw_planck_function_array)
      if (allocated(self%mrw_bol_planck_function_array)) deallocate(self%mrw_bol_planck_function_array)
      if (allocated(self%mrw_cum_planck_function)) deallocate(self%mrw_cum_planck_function)
      if (allocated(self%mrw_norm_planck_function)) deallocate(self%mrw_norm_planck_function)

      return
   
   end subroutine
   
   ! gives the specific dust mass extinction (cm^2/g) given wavelength (micron)
   elemental function dust_mass_ext(self,lambda) result (value)
   
      ! argument declaration
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength (micron)
      
      ! result declaration
      real(kind=rel_kind) :: value                            ! dust mass extinction (cm^2/g)
      
      value=lookup_and_interpolate(lambda,self%wavelength_array,self%dust_mass_ext_array)
      
      return
   
   end function
   
   ! gives the dust albedo given wavelength (micron)
   elemental function albedo(self,lambda) result (value)
   
      ! argument declaration
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength (micron)
      
      ! result declaration
      real(kind=rel_kind) :: value                            ! dust albedo
      
      value=lookup_and_interpolate(lambda,self%wavelength_array,self%albedo_array)
      
      return
   
   end function
   
   elemental subroutine dust_mass_ext_and_albedo(self,lambda,kappa_ext,a)
   
      ! argument declaration
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength (micron)
      real(kind=rel_kind),intent(out) :: kappa_ext            ! dust opacity
      real(kind=rel_kind),intent(out) :: a                    ! albedo
      
      ! variable declaration
      real(kind=rel_kind) :: value(2)                         ! temp values
      
      value=lookup_and_interpolate(lambda,self%wavelength_array,&
      &self%dust_mass_ext_array(:),self%albedo_array)
      
      kappa_ext=value(1)
      a=value(2)
      
      return
   
   end subroutine
   
   ! gives the mean scattering cosine given wavelength (micron)
   elemental function mean_scatter_cosine(self,lambda) result (value)
   
      ! argument declaration
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength (micron)
      
      ! result declaration
      real(kind=rel_kind) :: value                            ! mean cosine

      value=lookup_and_interpolate(lambda,self%wavelength_array,self%mean_cos_theta_array)
      
      return
   
   end function
   
   ! gives the monochromatic mass emissivity (B(T,lambda)*kappa(lambda)) (erg s^-1 sr^-1 g^-1 micron^-1) given wavelength (micron) and abs_rate (erg s^-1 g^-1)
   elemental function mono_mass_emissivity(self,lambda,abs_rate) result (value)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength (micron)
      real(kind=rel_kind),intent(in) :: abs_rate              ! energy absorption rate
      
      ! result declaration
      real(kind=rel_kind) :: value                            ! luminosity

      value=lookup_and_interpolate(lambda,abs_rate*0.25_rel_kind/pi,&
         &self%wavelength_array,self%bol_mass_emissivity_array,self%mono_mass_emissivity_array)
      
      return
   
   end function
   
   ! gives the bolometric mass emissivity (B(T)*kappa_p) (erg s^-1 sr^-1 g^-1) given temperature (K)
   elemental function bol_mass_emissivity(self,t) result (value)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: t                     ! temperature (K)
      
      ! result declaration
      real(kind=rel_kind) :: value                            ! luminosty

      value=lookup_and_interpolate(t,self%temperature_array,self%bol_mass_emissivity_array)
      
      return
   
   end function
   
   ! returns phase function value, given scattering cosine and wavelength (micron)
   elemental function phase_func(self,cos_theta,lambda) result (value)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: cos_theta             ! cos theta
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength (micron)
      
      ! result declaration
      real(kind=rel_kind) :: value                            ! phase function
      
      ! variable declaration
      real(kind=rel_kind) :: mean_cos_theta                   ! mean cos theta value
      
      mean_cos_theta=self%mean_scatter_cosine(lambda)

      value=0.5_rel_kind*(1._rel_kind-mean_cos_theta**2)/&
         &(1._rel_kind+mean_cos_theta**2-2._rel_kind*mean_cos_theta*&
         &cos_theta)**(1.5_rel_kind)
      
      return
   
   end function
   
   ! gives the dust emission temperature (K) given local absorption rate (ergs s^-1 g^-1)
   elemental function dust_temperature(self,abs_rate) result(value)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: abs_rate              ! dust absorption rate
      
      ! result declaration
      real(kind=rel_kind) :: value                            ! emission temperature
      
      value=lookup_and_interpolate(abs_rate*0.25_rel_kind/pi,&
         &self%bol_mass_emissivity_array,self%temperature_array)
      
      return
   
   end function
   
   ! gives a random wavelength (micron) from dust spectral radiance given absorption rate (erg s^-1 g^-1)
   function random_wavelength(self,abs_rate) result(lambda)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: abs_rate              ! energy absorption rate
      
      ! result declaration
      real(kind=rel_kind) :: lambda                           ! wavelength (micron)
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation index
      real(kind=rel_kind) :: r_num                            ! random number
      real(kind=rel_kind) :: temp_abs_rate                    ! temporary energy absorption rate
      real(kind=rel_kind) :: temp_norm_emissivity_array(size(self%wavelength_array))
      real(kind=rel_kind) :: temp_cum_emissivity_array(size(self%wavelength_array))
      
      ! enforce that abs_rate/4pi is in range of bol_mass_emissivity_array
      temp_abs_rate=max(min(abs_rate*0.25_rel_kind/pi,&
         &self%bol_mass_emissivity_array(size(self%bol_mass_emissivity_array))),&
         &self%bol_mass_emissivity_array(1))
      
      ! find temperature indices
      i_0=binary_search(temp_abs_rate,self%bol_mass_emissivity_array)
      
      
      temp_norm_emissivity_array=lerp(temp_abs_rate,&
         &self%bol_mass_emissivity_array(i_0),self%bol_mass_emissivity_array(i_0+1),&
         &self%norm_mono_mass_emissivity(:,i_0),self%norm_mono_mass_emissivity(:,i_0+1))   
      temp_cum_emissivity_array=lerp(temp_abs_rate,&
         &self%bol_mass_emissivity_array(i_0),self%bol_mass_emissivity_array(i_0+1),&
         &self%cum_mono_mass_emissivity(:,i_0),self%cum_mono_mass_emissivity(:,i_0+1))
      
      ! get random wavelength for t_0 and t_1
      call random_number(r_num)
      lambda=lookup_and_inv_interpolate(r_num,temp_cum_emissivity_array,&
         &self%wavelength_array,temp_norm_emissivity_array)
      
      return
   
   end function
   
   ! gives a random scattering cosine from phase function, given wavelength (micron
   function random_scatter_cosine(self,lambda) result(cos_theta)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength (microm)
      
      ! result declaration
      real(kind=rel_kind) :: cos_theta                        ! wavelength scattering cosine
      
      ! variable declarations
      real(kind=rel_kind) :: mean_cos_theta                   ! mean cos theta value
      real(kind=rel_kind) :: r_num                            ! random number
      
      ! get mean cos theta
      mean_cos_theta=self%mean_scatter_cosine(lambda)
      call random_number(r_num)
      
      if (abs(mean_cos_theta)>0._rel_kind) then
      
         ! Henyey Greenstein
         cos_theta=(1._rel_kind+mean_cos_theta**2-&
            &((1._rel_kind-mean_cos_theta**2)/(1._rel_kind+mean_cos_theta*&
            &(2._rel_kind*r_num-1._rel_kind)))**2)/(2._rel_kind*mean_cos_theta)
      else
         ! isotropic
         cos_theta=2._rel_kind*r_num-1._rel_kind
      end if
      
      return
   
   end function
   
   ! generate an isotropic random direction
   function random_emission_direction() result (new_direction)
   
      ! result declaration
      real(kind=rel_kind) :: new_direction(n_dim)     ! new direction after reemission
      
      call random_direction(new_direction)
      
      return
   
   end function
   
   
   ! gives the direction of a ray after a dust scattering event
   function random_scatter_direction(self,old_direction,lambda) result (new_direction)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: old_direction(n_dim)  ! direction before scatter event
      real(kind=rel_kind),intent(in) :: lambda                ! wavelength in m
      
      ! result declaration
      real(kind=rel_kind) :: new_direction(n_dim)             ! direction after scattering event
      
      ! variable declarations
      real(kind=rel_kind) :: cos_theta                        ! scattering angles
      
      
      
      ! get a random cos theta from phase function
      cos_theta=self%random_scatter_cosine(lambda)

      call scattered_direction_3d(new_direction,old_direction,cos_theta)
      
      return
   
   end function
   
   ! initialise Modified Random Walk arrays
   subroutine mrw_initialise(self,n_mrw)
   
      ! argument declarations
      class(dust),intent(inout) :: self                        ! dust object
      integer(kind=int_kind),intent(in) :: n_mrw               ! number of sample points
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                            ! counter
      
      ! allocate arrays
      allocate(self%mrw_zeta_array(n_mrw))
      allocate(self%mrw_y_array(n_mrw))
      allocate(self%mrw_planck_abs_array(size(self%temperature_array)))
      allocate(self%mrw_inv_planck_ext_array(size(self%temperature_array)))
      allocate(self%mrw_planck_sca_array(size(self%wavelength_array),size(self%temperature_array)))
      allocate(self%mrw_planck_function_array(size(self%wavelength_array),size(self%temperature_array)))
      allocate(self%mrw_bol_planck_function_array(size(self%temperature_array)))
      allocate(self%mrw_cum_planck_function(size(self%wavelength_array),size(self%temperature_array)))
      allocate(self%mrw_norm_planck_function(size(self%wavelength_array),size(self%temperature_array)))
            
      ! create y and zeta array
      self%mrw_y_array=lin_space(0._rel_kind,1._rel_kind,n_mrw)
      self%mrw_zeta_array(1)=0._rel_kind
      
      do i=2,n_mrw
         
         ! perform sum
         self%mrw_zeta_array(i)=0._rel_kind
         do j=1,digits(0._rel_kind)-1
            self%mrw_zeta_array(i)=self%mrw_zeta_array(i)+(-1._rel_kind)**(j+1)*self%mrw_y_array(i)**(j**2)
         end do
         
         ! average of odd and even sum
         self%mrw_zeta_array(i)=min(2._rel_kind*self%mrw_zeta_array(i)+&
            &(-1._rel_kind)**(digits(0._rel_kind)+1)*self%mrw_y_array(i)**(digits(0._rel_kind)**2),1._rel_kind)
            
         ! ensure values are monotonic
         self%mrw_zeta_array(i)=max(self%mrw_zeta_array(i),self%mrw_zeta_array(i-1))
         
      end do
      
      ! create remaining arrays
      do concurrent (i=1:size(self%temperature_array))
      
         ! planck function
         self%mrw_planck_function_array(:,i)=planck_lambda(self%wavelength_array,self%temperature_array(i))
         self%mrw_bol_planck_function_array(i)=trapz_intgr(self%wavelength_array,self%mrw_planck_function_array(:,i))
         self%mrw_cum_planck_function(:,i)=cum_dist_func(self%wavelength_array,self%mrw_planck_function_array(:,i))
         self%mrw_norm_planck_function(:,i)=self%mrw_planck_function_array(:,i)/self%mrw_bol_planck_function_array(i)
         
         ! planck abs/ext/sca
         self%mrw_planck_abs_array=trapz_intgr(self%wavelength_array,self%dust_mass_ext_array*(1._rel_kind-self%albedo_array)*&
            &self%mrw_planck_function_array(:,i))/self%mrw_bol_planck_function_array(i)
         
         self%mrw_inv_planck_ext_array=self%mrw_bol_planck_function_array(i)/&
            &trapz_intgr(self%wavelength_array,(1._rel_kind/self%dust_mass_ext_array)*self%mrw_planck_function_array(:,i))
         
         self%mrw_planck_sca_array(:,i)=cum_dist_func(self%wavelength_array,self%dust_mass_ext_array*self%albedo_array*&
            &self%mrw_planck_function_array(:,i),.false._log_kind)/self%mrw_bol_planck_function_array(i)
         
      end do
      
      return
   
   end subroutine
   
   ! sample random wavelength from planck function
   function mrw_random_wavelength(self,abs_rate) result (lambda)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: abs_rate              ! energy absorption rate
      
      ! result declaration
      real(kind=rel_kind) :: lambda                           ! wavelength (micron)
      
      ! variable declarations
      integer(kind=int_kind) :: i_0                           ! interpolation index
      real(kind=rel_kind) :: r_num                            ! random number
      real(kind=rel_kind) :: temp_abs_rate                    ! temporary energy absorption rate
      real(kind=rel_kind) :: temp_norm_emissivity_array(size(self%wavelength_array))
      real(kind=rel_kind) :: temp_cum_emissivity_array(size(self%wavelength_array))
      
      ! enforce that abs_rate/4pi is in range of bol_mass_emissivity_array
      temp_abs_rate=max(min(abs_rate*0.25_rel_kind/pi,&
         &self%bol_mass_emissivity_array(size(self%bol_mass_emissivity_array))),&
         &self%bol_mass_emissivity_array(1))
      
      ! find temperature indices
      i_0=binary_search(temp_abs_rate,self%bol_mass_emissivity_array)
      
      
      temp_norm_emissivity_array=lerp(temp_abs_rate,&
         &self%bol_mass_emissivity_array(i_0),self%bol_mass_emissivity_array(i_0+1),&
         &self%mrw_norm_planck_function(:,i_0),self%mrw_norm_planck_function(:,i_0+1))   
      temp_cum_emissivity_array=lerp(temp_abs_rate,&
         &self%bol_mass_emissivity_array(i_0),self%bol_mass_emissivity_array(i_0+1),&
         &self%mrw_cum_planck_function(:,i_0),self%mrw_cum_planck_function(:,i_0+1))
      
      ! get random wavelength for t_0 and t_1
      call random_number(r_num)
      lambda=lookup_and_inv_interpolate(r_num,temp_cum_emissivity_array,&
         &self%wavelength_array,temp_norm_emissivity_array)
         
      return
   
   end function
   
   ! get planck absorption coefficient
   elemental function mrw_planck_abs(self,abs_rate) result (planck_abs)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: abs_rate              ! energy absorption rate
      
      ! result declaration
      real(kind=rel_kind) :: planck_abs                       ! planck absorption coefficient
         
      planck_abs=lookup_and_interpolate(abs_rate,self%bol_mass_emissivity_array,self%mrw_planck_abs_array)
      
      return
   
   end function
   
   ! get inverse planck extinction coefficient
   elemental function mrw_inv_planck_ext(self,abs_rate) result (inv_planck_ext)
   
      ! argument declarations
      class(dust),intent(in) :: self                          ! dust object
      real(kind=rel_kind),intent(in) :: abs_rate              ! energy absorption rate
      
      ! result declaration
      real(kind=rel_kind) :: inv_planck_ext                   ! inverse planck extinction coefficient
         
      inv_planck_ext=lookup_and_interpolate(abs_rate,self%bol_mass_emissivity_array,self%mrw_inv_planck_ext_array)
      
      return
   
   end function
   
   ! get planck scattering coefficient between wavelengths
   elemental function mrw_planck_sca(self,abs_rate,wavelength_0,wavelength_1) result (planck_sca)
   
      ! argument declarations
      class(dust),intent(in) :: self                        ! dust object
      real(kind=rel_kind),intent(in) :: abs_rate            ! energy absorption
      real(kind=rel_kind),intent(in) :: wavelength_0        ! wavelength lower bound
      real(kind=rel_kind),intent(in) :: wavelength_1        ! wavelength upper bound
      
      ! result declaration
      real(kind=rel_kind) :: planck_sca                     ! planck scattering coefficient
      
      planck_sca=lookup_and_interpolate(wavelength_1,abs_rate*0.25_rel_kind/pi,&
         &self%wavelength_array,self%bol_mass_emissivity_array,self%mrw_planck_sca_array)-&
         &lookup_and_interpolate(wavelength_0,abs_rate*0.25_rel_kind/pi,&
         &self%wavelength_array,self%bol_mass_emissivity_array,self%mrw_planck_sca_array)
         
      return
      
   end function
   
   ! get mrw distance travelled through optically thick isothermal sphere
   function mrw_distance(self,abs_rate,radius,density) result (distance)
   
      ! argument declarations
      class(dust),intent(in) :: self                        ! dust object
      real(kind=rel_kind),intent(in) :: abs_rate            ! energy absorption
      real(kind=rel_kind),intent(in) :: radius              ! radius of isothermal sphere
      real(kind=rel_kind),intent(in) :: density             ! density of isothermal sphere
      
      ! result declaration
      real(kind=rel_kind) :: distance                       ! distance travelled by luminosity packet
      
      ! variable declarations
      real(kind=rel_kind) :: r_num                          ! random number
      
      ! get random number
      call random_number(r_num)
      
      ! calculate distance travelled through isothermal sphere
      distance=-log(lookup_and_interpolate(r_num,self%mrw_zeta_array,self%mrw_y_array))*(radius/pi)**2*&
         &3._rel_kind*density*self%mrw_inv_planck_ext(abs_rate)
      
      return
   
   end function
   
   ! gives the direction of a ray after emerging from isothermal sphere
   function mrw_random_direction(old_direction) result (new_direction)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: old_direction(n_dim)  ! direction before scatter event
      
      ! result declaration
      real(kind=rel_kind) :: new_direction(n_dim)             ! direction after scattering event
      
      ! variable declarations
      real(kind=rel_kind) :: cos_theta                        ! scattering angles
      
      ! cos theta from lambertian surface
      call random_number(cos_theta)
      cos_theta=sqrt(cos_theta)

      call scattered_direction_3d(new_direction,old_direction,cos_theta)
      
      return
   
   end function

end module