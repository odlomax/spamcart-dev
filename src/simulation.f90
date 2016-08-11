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


! mcrt simulation module
module m_simulation

   use m_kind_parameters
   use m_constants_parameters
   use m_dust
   use m_dust_d03
   use m_source
   use m_source_point_bb
   use m_source_external_bb
   use m_source_external_ps05
   use m_particle
   use m_binary_tree
   use m_kernel
   use m_kernel_m4
   use m_ray
   use m_options
   use m_datacube
   use m_io
   use m_string
   use omp_lib
   
   implicit none
   
   private
   
   public :: simulation
   
   ! define mcrt simulation class
   type :: simulation
   
      type(particle),pointer,contiguous :: particle_array(:)      ! array of sph particles
      type(binary_tree),pointer :: sph_tree                       ! tree of particles
      type(ray),pointer :: lum_packet_array(:)                    ! array of luminosity packets
      type(options),pointer :: sim_params                         ! simulation parameters
      type(datacube),pointer :: intensity_cube                    ! datacube
      class(kernel),pointer :: sph_kernel                         ! smoothing kernel
      class(dust),pointer :: dust_prop                            ! interstellar dust properties
      class(source),pointer :: point_source_array(:)              ! array of point sources (stars)
      class(source),pointer :: isrf_prop                          ! interstellar radiation field properties
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: destroy
      procedure,non_overridable :: perform_iteration
      procedure,non_overridable :: make_datacube
   
   end type
   
   contains
   
   ! initialise simulation object
   subroutine initialise(self,sim_params)
   
      ! argument declarations
      class(simulation),intent(inout) :: self                     ! simulation object
      type(options),intent(in),target :: sim_params               ! simulation parameters
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                               ! counter
      real(kind=rel_kind),allocatable :: position(:,:)            ! particle positions
      real(kind=rel_kind),allocatable :: mass(:)                  ! particle mass
      real(kind=rel_kind),allocatable :: temperature(:)           ! particle temperature
      real(kind=rel_kind),allocatable :: luminosity(:)            ! particle luminosity
      
      ! associate parameters
      self%sim_params=>sim_params
   
      write(*,"(A)") "initialise kernel"
      allocate(kernel_m4::self%sph_kernel)
      call self%sph_kernel%initialise()
      
      write(*,"(A)") "initialise dust"
      allocate(dust_d03::self%dust_prop)
      call self%dust_prop%initialise(self%sim_params%dust_r_v,self%sim_params%dust_t_min,self%sim_params%dust_t_max,&
         &self%sim_params%dust_n_t,iso_scatter=self%sim_params%dust_iso_scatter)
         
      if (self%sim_params%sim_mrw) then
         write(*,"(A)") "initialise modified random walk"
         call self%dust_prop%mrw_initialise(self%sim_params%sim_n_mrw)
      end if
      
      write(*,"(A)") "read in particles"
      call read_in_sph_particles_3d(position,mass,temperature,self%sim_params%sim_cloud_file)
      
      if (.not.self%sim_params%sim_restart) temperature=self%sim_params%sim_initial_t

      allocate(self%particle_array(size(position,2)))
      write(*,"(A)") "initialise particles"
      do i=1,size(self%particle_array)
         if (self%sim_params%sph_scattered_light) then
            call self%particle_array(i)%initialise(position(:,i),(/(0._rel_kind,j=1,n_dim)/),mass(i),&
               &4._rel_kind*pi*self%dust_prop%bol_mass_emissivity(temperature(i)),&
               &self%sim_params%sph_lambda_min,self%sim_params%sph_lambda_max,self%sim_params%sph_n_lambda)
         else
            call self%particle_array(i)%initialise(position(:,i),(/(0._rel_kind,j=1,n_dim)/),mass(i),&
               &4._rel_kind*pi*self%dust_prop%bol_mass_emissivity(temperature(i)))
         end if
      end do
      
      write(*,"(A)") "initialise tree"
      allocate(self%sph_tree)
      call self%sph_tree%initialise(self%particle_array,self%sph_kernel,self%sim_params%sph_eta,self%sim_params%sim_min_d)
      
      if (self%sim_params%point_sources) then
      
         write(*,"(A)") "read in point sources"
         call read_in_point_sources_3d(position,luminosity,temperature,self%sim_params%sim_star_file)
         allocate(source_point_bb::self%point_source_array(size(position,2)))
         
         write(*,"(A)") "initialise point sources"
         do i=1,size(self%point_source_array)
            call self%point_source_array(i)%initialise(position(:,i),(/(0._rel_kind,j=1,n_dim)/),temperature=temperature(i),&
               &luminosity=luminosity(i),wavelength_min=self%sim_params%point_lambda_min,&
               &wavelength_max=self%sim_params%point_lambda_max,n_wavelength=self%sim_params%point_n_lambda)
         end do

      else
         self%point_source_array=>null()
      end if     
      
      write(*,"(A)") "initialise radiation field"
      if (self%sim_params%ext_rf) then
         select case (trim(self%sim_params%ext_rf_type))
      
            case ("ps05")
         
               allocate(source_external_ps05::self%isrf_prop)
               call self%isrf_prop%initialise(self%sph_tree%com,(/(0._rel_kind,i=1,n_dim)/),&
                  &radius=0.5_rel_kind*self%sph_tree%max_length,gal_r=self%sim_params%ext_rf_gal_r,&
                  &gal_z=self%sim_params%ext_rf_gal_r,add_cmb=self%sim_params%ext_rf_cmb)
         
            case ("bb")
            
               allocate(source_external_bb::self%isrf_prop)
               call self%isrf_prop%initialise(self%sph_tree%com,(/(0._rel_kind,i=1,n_dim)/),&
                  &radius=0.5_rel_kind*self%sph_tree%max_length,wavelength_min=self%sim_params%ext_rf_lambda_min,&
                  &wavelength_max=self%sim_params%ext_rf_lambda_max,n_wavelength=self%sim_params%ext_rf_n_lambda,&
                  &temperature=self%sim_params%ext_rf_t,dilution=self%sim_params%ext_rf_d,&
                  &add_cmb=self%sim_params%ext_rf_cmb)
         
            case default
         
               write(*,"(A)") "Unknown radiation field type:",self%sim_params%ext_rf_type
               stop
      
         end select
         
      else
         self%isrf_prop=>null()
      end if

      return
   
   end subroutine
   
   pure subroutine destroy(self)
   
      ! argument declarations
      class(simulation),intent(inout) :: self                     ! simulation object
      
      call self%particle_array%destroy()
      call self%sph_tree%destroy()
      call self%dust_prop%destroy()
      if (associated(self%point_source_array)) call self%point_source_array%destroy()
      if (associated(self%isrf_prop)) call self%isrf_prop%destroy()
      if (associated(self%intensity_cube)) call self%intensity_cube%destroy()
      
      deallocate(self%particle_array)
      deallocate(self%sph_tree)
      deallocate(self%sim_params)
      deallocate(self%dust_prop)
      if (associated(self%point_source_array)) deallocate(self%point_source_array)
      if (associated(self%isrf_prop)) deallocate(self%isrf_prop)
      if (associated(self%intensity_cube)) deallocate(self%intensity_cube)
      
      return
   
   end subroutine
   
   ! perform temperature calculation
   subroutine perform_iteration(self,iteration)
   
      ! argument declarations
      class(simulation),intent(inout) :: self                     ! simulation object
      integer(kind=int_kind),intent(in) :: iteration              ! iteration number
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                               ! counter
      integer(kind=int_kind) :: n_packets_source                  ! number of packets for individual source
      integer(kind=int_kind) :: n_threads                         ! number of available threads
      integer(kind=int_kind) :: thread_num                        ! thread number
      real(kind=rel_kind) :: total_luminosity                     ! total luminosity of all sources
      real(kind=rel_kind) :: luminosity_chunk                     ! luminosity chunk for each source
      real(kind=rel_kind) :: position(n_dim)                      ! position vector
      real(kind=rel_kind),allocatable :: position_array(:,:)      ! array of positions
      character(kind=chr_kind,len=string_length) :: output_file   ! output file        
      
      ! calculate total luminosity
      total_luminosity=0._rel_kind
      if (associated(self%point_source_array)) total_luminosity=total_luminosity+sum(self%point_source_array%luminosity)
      if (associated(self%isrf_prop)) total_luminosity=total_luminosity+self%isrf_prop%luminosity
      
      ! get number of threads
      n_threads=omp_get_num_procs()
      
      ! initialise luminosity packets (one per thread)
      write(*,"(A)") "initialise luminosity packets"
      allocate(self%lum_packet_array(n_threads))
      do i=1,n_threads
         call self%lum_packet_array(i)%initialise(self%sph_kernel,self%sph_tree,self%dust_prop)
      end do
      
      ! reset scattered light bins and sublimation fraction
      call self%particle_array%reset_a_dot(4._rel_kind*pi*self%dust_prop%bol_mass_emissivity(self%sim_params%dust_sub_t_min),&
         &4._rel_kind*pi*self%dust_prop%bol_mass_emissivity(self%sim_params%dust_sub_t_max))
      
      if (associated(self%point_source_array)) then
      
         write(*,"(A)") "follow packets from point sources"
         !$omp parallel num_threads(n_threads) default(shared) private(i,j,thread_num,n_packets_source,luminosity_chunk)
            do i=1,size(self%point_source_array)
         
               thread_num=omp_get_thread_num()+1
               n_packets_source=int(real(self%sim_params%sim_n_packet,rel_kind)*self%point_source_array(i)%luminosity/&
                  &(total_luminosity*real(n_threads,rel_kind)))
               luminosity_chunk=self%point_source_array(i)%luminosity/real(n_packets_source*n_threads,rel_kind)
                  
               do j=1,n_packets_source
                  call self%lum_packet_array(thread_num)%follow(self%point_source_array(i)%position,&
                     &self%point_source_array(i)%random_direction(),self%point_source_array(i)%velocity,&
                     &self%point_source_array(i)%random_wavelength(),luminosity_chunk,&
                     &self%sim_params%sim_mrw,self%sim_params%sim_mrw_gamma)
               
                  if (mod(j,min(1000,n_packets_source))==0.and.thread_num==1) &
                     &write(*,"(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)") &
                        &"iteration ",iteration," of ",self%sim_params%sim_n_it,&
                        &", star ",i," of ",size(self%point_source_array),&
                        &", packet ",j*n_threads," of ",n_packets_source*n_threads
               
               end do
                       
            end do
         !$omp end parallel
      
      end if
      
      if (associated(self%isrf_prop)) then
      
         write(*,"(A)") "follow packets from external radiation field"
         n_packets_source=int(real(self%sim_params%sim_n_packet,rel_kind)*self%isrf_prop%luminosity/&
            &(total_luminosity*real(n_threads,rel_kind)))
         luminosity_chunk=self%isrf_prop%luminosity/real(n_packets_source*n_threads,rel_kind)
         !$omp parallel num_threads(n_threads) default(shared) private(j,thread_num,position)
         
               thread_num=omp_get_thread_num()+1
                  
               do j=1,n_packets_source
                  position=self%isrf_prop%random_position()
                  call self%lum_packet_array(thread_num)%follow(position,&
                     &self%isrf_prop%random_direction(position),self%isrf_prop%velocity,&
                     &self%isrf_prop%random_wavelength(),luminosity_chunk,&
                     &self%sim_params%sim_mrw,self%sim_params%sim_mrw_gamma)
               
                  if (mod(j,min(1000,n_packets_source))==0.and.thread_num==1) &
                     &write(*,"(A,I0,A,I0,A,I0,A,I0)") &
                        &"iteration ",iteration," of ",self%sim_params%sim_n_it,&
                        &", external rad field, packet ",j*n_threads," of ",n_packets_source*n_threads
               
               end do
                       
         !$omp end parallel
      
      end if
      
      ! normalise absorption rate
      call self%particle_array%normalise_a()
      
      ! write out particles
      write(*,"(A)") "write out particles"
      write(output_file,"(A,I2.2,A)") trim(self%sim_params%sim_id)//"_",iteration,".dat"
      allocate(position_array(n_dim,size(self%particle_array)))
      do i=1,size(self%particle_array)
         position_array(:,i)=self%particle_array(i)%r
      end do
      call write_out_sph_particles_3d(position_array,self%particle_array%m,&
         &self%dust_prop%dust_temperature(self%particle_array%a_dot),self%particle_array%rho,self%particle_array%h,output_file)
      
      ! deallocate luminosity packets
      call self%lum_packet_array%destroy()
      deallocate(self%lum_packet_array)
   
      return
   
   end subroutine
   
   subroutine make_datacube(self)
   
      ! argument declarations
      class(simulation),intent(inout) :: self                     ! simulation object
      
      ! allocate datacube
      allocate(self%intensity_cube)
      
      ! make datacube
      write(*,"(A)") "making datacube"
      call self%intensity_cube%initialise(self%sph_kernel,self%sph_tree,self%dust_prop,&
         &self%sim_params%datacube_x_min,self%sim_params%datacube_x_max,self%sim_params%datacube_n_x,&
         &self%sim_params%datacube_y_min,self%sim_params%datacube_y_max,self%sim_params%datacube_n_y,&
         &self%sim_params%datacube_angles,string_to_real(self%sim_params%datacube_lambda_string),&
         &background=self%isrf_prop,point_sources=self%point_source_array)
         
      if (self%sim_params%datacube_convolve) then
         write(*,"(A)") "convolving datacube with psf"
         call self%intensity_cube%psf_convolve&
            &(self%sim_params%datacube_distance,string_to_real(self%sim_params%datacube_fwhm_string))
      end if
      
      write(*,"(A)") "write out datacube"
      call  write_out_datacube_3d(self%intensity_cube%x,self%intensity_cube%y,self%intensity_cube%lambda,&
         &self%intensity_cube%i_lambda,self%sim_params%sim_id)
      
      return
   
   end subroutine

end module