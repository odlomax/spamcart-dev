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
   use m_atomic_update
   use m_dust
   use m_dust_d03
   use m_source
   use m_source_point_bb
   use m_source_point_ms_star
   use m_source_point_custom
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
      procedure,non_overridable :: source_extract
      procedure,non_overridable :: make_datacube
      procedure,non_overridable :: write_particles_bin
      procedure,non_overridable :: read_particles_bin
   
   end type
   
   contains
   
   ! initialise simulation object
   subroutine initialise(self,sim_params)
   
      ! argument declarations
      class(simulation),intent(inout) :: self                     ! simulation object
      type(options),intent(in),target :: sim_params               ! simulation parameters
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                               ! counter
      integer(kind=int_kind) :: string_len                        ! length of string
      logical(kind=log_kind) :: h_present                         ! is smoothing length present in cloud file
      real(kind=rel_kind),allocatable :: position(:,:)            ! particle/star positions
      real(kind=rel_kind),allocatable :: mass(:)                  ! particle/star mass
      real(kind=rel_kind),allocatable :: temperature(:)           ! particle/star temperature
      real(kind=rel_kind),allocatable :: luminosity(:)            ! star luminosity
      real(kind=rel_kind),allocatable :: age(:)                   ! star age
      real(kind=rel_kind),allocatable :: metallicity(:)           ! star metallicity
      character(kind=chr_kind,len=string_length) :: file_name     ! name of file
      
      ! associate parameters
      self%sim_params=>sim_params
   
      write(*,"(A)") "initialise kernel"
      allocate(kernel_m4::self%sph_kernel)
      call self%sph_kernel%initialise()      
      
      write(*,"(A)") "initialise dust"
      allocate(dust_d03::self%dust_prop)
      call self%dust_prop%initialise(self%sim_params%dust_r_v,self%sim_params%dust_t_min,self%sim_params%dust_t_max,&
         &self%sim_params%dust_n_t,iso_scatter=self%sim_params%dust_iso_scatter)
         
      write(*,"(A)") "initialise modified random walk"
      call self%dust_prop%mrw_initialise(self%sim_params%sim_n_mrw)
      
      write(*,"(A)") "read in particles"
      ! check if cloud input is ascii (.dat) or binary (.bin)
      string_len=len_trim(self%sim_params%sim_cloud_file)
      select case (self%sim_params%sim_cloud_file(string_len-2:string_len))
         case ("bin")
            call self%read_particles_bin()
            h_present=.true.
         
         case ("dat")
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
            h_present=.false.
      
      end select
      
      write(*,"(A)") "initialise tree"
      allocate(self%sph_tree)
      call self%sph_tree%initialise(self%particle_array,self%sph_kernel,self%sim_params%sph_eta,self%sim_params%sim_min_d,h_present)
      
      if (self%sim_params%point_sources) then
      
         write(*,"(A)") "read in point sources"
         select case (trim(self%sim_params%point_type))
         
            case ("bb")
            
               call read_in_point_sources_3d(position,luminosity=luminosity,temperature=temperature,&
                  &file_name=self%sim_params%sim_star_file)
               allocate(source_point_bb::self%point_source_array(size(position,2)))
         
               write(*,"(A)") "initialise point sources"
               do i=1,size(self%point_source_array)
                  call self%point_source_array(i)%initialise(position(:,i),(/(0._rel_kind,j=1,n_dim)/),temperature=temperature(i),&
                     &luminosity=luminosity(i),wavelength_min=self%sim_params%point_lambda_min,&
                     &wavelength_max=self%sim_params%point_lambda_max,n_wavelength=self%sim_params%point_n_lambda)
               end do
               
            case ("ms_star")
            
               call read_in_point_sources_3d(position,mass=mass,age=age,metallicity=metallicity,&
                  &file_name=self%sim_params%sim_star_file)
               allocate(source_point_ms_star::self%point_source_array(size(position,2)))
         
               write(*,"(A)") "initialise point sources"
               do i=1,size(self%point_source_array)
                  call self%point_source_array(i)%initialise(position(:,i),(/(0._rel_kind,j=1,n_dim)/),mass=mass(i),&
                     &age=age(i),metallicity=metallicity(i))
               end do
               
            case ("custom")
            
               call read_in_point_sources_3d(position,file_name=self%sim_params%sim_star_file)
               allocate(source_point_custom::self%point_source_array(size(position,2)))
               
               
               write(*,"(A)") "initialise point sources"
               open(2,file=trim(self%sim_params%point_file_list))
               read(2,*) file_name
               do i=1,size(self%point_source_array)
                  call self%point_source_array(i)%initialise(position(:,i),(/(0._rel_kind,j=1,n_dim)/),luminosity_file=file_name)
               end do
               close(2)
                
         end select
         
         ! set sublimation radii for point sources
         call self%point_source_array%set_sublimation_radius(self%dust_prop,self%sim_params%dust_sub_t)
         
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
      integer(kind=int_kind) :: i,j,k                             ! counter
      integer(kind=int_kind) :: n_threads                         ! number of available threads
      integer(kind=int_kind) :: n_packets_source                  ! number of packets for individual source
      real(kind=rel_kind) :: total_luminosity                     ! total luminosity of all sources
      real(kind=rel_kind) :: luminosity_chunk                     ! luminosity chunk for each source
      real(kind=rel_kind) :: position(n_dim)                      ! position vector
      real(kind=rel_kind) :: direction(n_dim)                     ! direction vector
      real(kind=rel_kind) :: wavelength                           ! wavelength
      real(kind=rel_kind) :: tau                                  ! optical depth travelled by luminosity packet
      real(kind=rel_kind) :: bg_mean_tau                          ! mean optical depth travelled by background packets
      real(kind=rel_kind) :: source_luminosity                    ! total luminosity of sources
      real(kind=rel_kind) :: packet_luminosity                    ! packet energy deposition rate
      real(kind=rel_kind) :: particle_luminosity                  ! particle energy emission rate
      real(kind=rel_kind),allocatable :: position_array(:,:)      ! array of positions
      real(kind=rel_kind),allocatable :: ps_mean_tau(:)           ! mean optical depth travelled by point source packets
      character(kind=chr_kind,len=string_length) :: output_file   ! output file        
      
      ! get number of threads
      n_threads=omp_get_num_procs()
      
      ! initialise luminosity packets (one per thread)
      write(*,"(A)") "initialise luminosity packets"
      allocate(self%lum_packet_array(n_threads))
      do i=1,n_threads
         call self%lum_packet_array(i)%initialise(self%sph_kernel,self%sph_tree,self%dust_prop)
      end do
      
      ! reset scattered light bins
      call self%particle_array%reset_a_dot()
      
      write(*,"(A)") "set sublimation fraction"
      ! set sublimation fraction
      ! by default, set to 1
      self%particle_array%f_sub=1._rel_kind
      select case (trim(self%sim_params%dust_sub_type))
      
         case ("variable")
      
            do i=1,size(self%particle_array)
               self%particle_array(i)%f_sub=self%sph_tree%root_node%sph_gather_f_sub_variable(self%sph_kernel,&
                  &self%particle_array(i),4._rel_kind*pi*self%dust_prop%bol_mass_emissivity(self%sim_params%dust_sub_t))
            end do
            
         case ("fixed")
         
            if (associated(self%point_source_array)) then
               do i=1,size(self%particle_array)
                  self%particle_array(i)%f_sub=self%sph_tree%root_node%sph_gather_f_sub_fixed(self%sph_kernel,&
                     &self%particle_array(i),self%point_source_array)
               end do
            end if
         
      end select
      
      if (associated(self%point_source_array)) then
      
         allocate(ps_mean_tau(size(self%point_source_array)))
         ps_mean_tau=0._rel_kind
         total_luminosity=sum(self%point_source_array%luminosity)
         
         write(*,"(A)") "follow packets from point sources"
         
         do i=1,size(self%point_source_array)
         
            if (self%sim_params%sim_equal_packets_per_point) then
               n_packets_source=self%sim_params%sim_n_packet_point/size(self%point_source_array)
            else
               n_packets_source=&
                  &int(real(self%sim_params%sim_n_packet_point,rel_kind)*self%point_source_array(i)%luminosity/total_luminosity)
            end if
            
            n_packets_source=max(1,n_packets_source)
         
            luminosity_chunk=self%point_source_array(i)%luminosity/real(n_packets_source,rel_kind)
            k=0
      
            !$omp parallel do num_threads(n_threads) default(shared) private(j,tau)
               do j=1,n_packets_source
               
                  tau=0._rel_kind
                  call self%lum_packet_array(omp_get_thread_num()+1)%follow(self%point_source_array(i)%position,&
                     &self%point_source_array(i)%random_direction(),self%point_source_array(i)%velocity,&
                     &self%point_source_array(i)%random_wavelength(),luminosity_chunk,&
                     &self%sim_params%sim_mrw,self%sim_params%sim_accurate_mrw,self%sim_params%sim_mrw_gamma,tau)
                     
                  call atomic_integer_add(k,1)
                  call atomic_real_add(ps_mean_tau(i),tau )
                           
                  if (mod(k,min(self%sim_params%sim_n_packet_point/100,n_packets_source))==0) &
                     &write(*,"(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)") &
                     &"iteration ",iteration," of ",self%sim_params%sim_n_it,&
                     &", star ",i," of ",size(self%point_source_array),&
                     &", packet ",k," of ",n_packets_source
               end do
            !$omp end parallel do
            ps_mean_tau(i)=ps_mean_tau(i)/real(n_packets_source,rel_kind)
                   
         end do
      
      end if
      
      if (associated(self%isrf_prop)) then
      
         bg_mean_tau=0._rel_kind
      
         write(*,"(A)") "follow packets from external radiation field"
         luminosity_chunk=self%isrf_prop%luminosity/real(self%sim_params%sim_n_packet_external,rel_kind)
         k=0
         !$omp parallel do num_threads(n_threads) default(shared) private(j,position,direction,wavelength,tau)
                  
            do j=1,self%sim_params%sim_n_packet_external
               position=self%isrf_prop%random_position()
               direction=self%isrf_prop%random_direction(position)
               wavelength=self%isrf_prop%random_wavelength()
               tau=0._rel_kind
               call self%lum_packet_array(omp_get_thread_num()+1)%follow(position,&
                  &direction,self%isrf_prop%velocity,&
                  &wavelength,luminosity_chunk,&
                  &self%sim_params%sim_mrw,self%sim_params%sim_accurate_mrw,self%sim_params%sim_mrw_gamma,tau)
            
               call atomic_integer_add(k,1)
               call atomic_real_add(bg_mean_tau,tau)
            
               if (mod(k,min(self%sim_params%sim_n_packet_external/100,self%sim_params%sim_n_packet_external))==0) &
                  &write(*,"(A,I0,A,I0,A,I0,A,I0)") &
                     &"iteration ",iteration," of ",self%sim_params%sim_n_it,&
                     &", external rad field, packet ",k," of ",self%sim_params%sim_n_packet_external
            
            end do
                       
         !$omp end parallel do
         bg_mean_tau=bg_mean_tau/real(self%sim_params%sim_n_packet_external,rel_kind)
         
      
      end if 
      
      ! normalise absorption rate
      call self%particle_array%normalise_a()
      
      ! write out packet travel stats
      if (associated(self%point_source_array)) then
         do i=1,size(ps_mean_tau)
            write(*,"(A,I0,A,E25.17)") "star ",i," mean escape tau: ",ps_mean_tau(i)
         end do
      end if
      
      if (associated(self%isrf_prop)) write(*,"(A,E25.17)") "external rad field escape tau:   ",bg_mean_tau      
         
      ! calculate packet energy deposition rate
      source_luminosity=0._rel_kind
      packet_luminosity=0._rel_kind
      if (associated(self%point_source_array)) then
         packet_luminosity=packet_luminosity+sum(self%point_source_array%luminosity*ps_mean_tau)
         source_luminosity=source_luminosity+sum(self%point_source_array%luminosity)
      end if
         
      if (associated(self%isrf_prop)) then
         packet_luminosity=packet_luminosity+self%isrf_prop%luminosity*bg_mean_tau
         source_luminosity=source_luminosity+self%isrf_prop%luminosity
      end if
      
      ! calculate particle luminosity
      particle_luminosity=0._rel_kind
      do i=1,size(self%particle_array)
         particle_luminosity=particle_luminosity+&
            &self%particle_array(i)%m*self%particle_array(i)%f_sub*self%particle_array(i)%a_dot
            if (associated(self%particle_array(i)%lambda_array)) particle_luminosity=particle_luminosity+&
               &self%particle_array(i)%m*self%particle_array(i)%f_sub*&
               &sum(self%particle_array(i)%a_dot_scatter_array*(self%particle_array(i)%lambda_array(2:)-&
               &self%particle_array(i)%lambda_array(:size(self%particle_array(i)%lambda_array)-1)))
      end do
      
      write(*,"(A,E25.17)") "source luminosity:   ",source_luminosity
      write(*,"(A,E25.17)") "packet luminosity:   ",packet_luminosity
      write(*,"(A,E25.17)") "particle luminosity: ",particle_luminosity
      
      
      ! write out particles
      write(*,"(A)") "write out particles"
      write(output_file,"(A,I2.2,A)") trim(self%sim_params%sim_id)//"/iteration_",iteration,".dat"
      allocate(position_array(n_dim,size(self%particle_array)))
      do i=1,size(self%particle_array)
         position_array(:,i)=self%particle_array(i)%r
      end do
      call write_out_sph_particles_3d(position_array,self%particle_array%m,&
         &self%dust_prop%dust_temperature(self%particle_array%a_dot),self%particle_array%rho,self%particle_array%h,output_file)
         
      call self%write_particles_bin(iteration)
      
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
      
      
      call self%intensity_cube%initialise(self%sim_params,self%isrf_prop,self%point_source_array,&
         &self%sph_kernel,self%sph_tree,self%dust_prop)
         
      if (self%sim_params%datacube_convolve) then
         write(*,"(A)") "convolving datacube with psf"
         call self%intensity_cube%psf_convolve()
      end if
      
      write(*,"(A)") "write out datacube"

      call self%intensity_cube%write_out()
      
      return
   
   end subroutine
   
   ! write out particles in binary format
   subroutine write_particles_bin(self,iteration)
      
      ! character declarations
      class(simulation),intent(in) :: self                     ! simulation object
      integer(kind=int_kind),intent(in) :: iteration           ! iteration number
      
      ! variable delcarations
      integer(kind=int_kind) :: i                              ! counter
      character(kind=chr_kind,len=string_length) :: file_name  ! file name
      
      ! open file
      write(file_name,"(A,I2.2,A)") trim(self%sim_params%sim_id)//"/iteration_",iteration,".bin"
      open(1,file=trim(file_name),form="unformatted")
      
      ! write number of particles to file
      write(1) size(self%particle_array,kind=int_kind)
      
      do i=1,size(self%particle_array)
      
         ! write out particle properties
         
         write(1) self%particle_array(i)%r
         write(1) self%particle_array(i)%v
         write(1) self%particle_array(i)%m
         write(1) self%particle_array(i)%h
         write(1) self%particle_array(i)%rho
         write(1) self%particle_array(i)%a_dot
         write(1) self%particle_array(i)%f_sub
         if (associated(self%particle_array(i)%lambda_array)) then
            write(1) size(self%particle_array(i)%lambda_array,kind=int_kind)
            write(1) self%particle_array(i)%lambda_array
            write(1) self%particle_array(i)%a_dot_scatter_array
         else
            write(1) 0
         end if
         
      end do
      
      close(1)
      
      return
   
   end subroutine
   
   ! write out particles in binary format
   subroutine read_particles_bin(self)
      
      ! character declarations
      class(simulation),intent(inout) :: self                  ! simulation object
      
      ! variable delcarations
      integer(kind=int_kind) :: i                              ! counter
      integer(kind=int_kind) :: n_particles                    ! number of particles
      integer(kind=int_kind) :: n_array                        ! arbitrary array length
      
      ! open file       
      open(1,file=trim(self%sim_params%sim_cloud_file),form="unformatted")
      
      ! write number of particles to file
      read(1) n_particles
      allocate(self%particle_array(n_particles))
      
      do i=1,size(self%particle_array)
      
         ! write out particle properties
         
         read(1) self%particle_array(i)%r
         read(1) self%particle_array(i)%v
         read(1) self%particle_array(i)%m
         read(1) self%particle_array(i)%h
         read(1) self%particle_array(i)%rho
         read(1) self%particle_array(i)%a_dot
         read(1) self%particle_array(i)%f_sub
         read(1) n_array
         if (n_array>0) then
            allocate(self%particle_array(i)%lambda_array(n_array))
            allocate(self%particle_array(i)%a_dot_scatter_array(n_array-1))
            read(1) self%particle_array(i)%lambda_array
            read(1) self%particle_array(i)%a_dot_scatter_array
         else
            self%particle_array(i)%lambda_array=>null()
            self%particle_array(i)%a_dot_scatter_array=>null()
         end if
         
         ! calculate some derived quantities
         self%particle_array(i)%inv_h=1._rel_kind/self%particle_array(i)%h
         self%particle_array(i)%inv_rho=1._rel_kind/self%particle_array(i)%rho
         self%particle_array(i)%a_dot_new=0._rel_kind
         
      end do
      
      close(1)
      
      return
   
   end subroutine
   
   ! cull particle ensemble to sphere
   subroutine source_extract(self,centre,radius)
   
      ! argument declarations
      class(simulation),intent(inout) :: self            ! simulation object
      real(kind=rel_kind),intent(in) :: centre(n_dim)    ! source centre
      real(kind=rel_kind),intent(in) :: radius           ! radius
      
      ! variable declarations
      type(particle),pointer :: temp_particle_array(:)   ! temporary particle array
      integer(kind=int_kind) :: i                        ! counter
      logical(kind=log_kind) :: array_mask(size(self%particle_array))   ! true if particles are within sphere
      
      ! set mask
      do i=1,size(self%particle_array)
         array_mask(i)=sum((self%particle_array(i)%r-centre)**2)<=radius**2
         self%particle_array(i)%r=self%particle_array(i)%r-centre

      end do
      
      ! copy live particles
      allocate(temp_particle_array(count(array_mask)))
      
      temp_particle_array=pack(self%particle_array,array_mask)
      
      self%particle_array=>temp_particle_array
      
      ! reset particle tree
      call self%sph_tree%destroy()
      call self%sph_tree%initialise(&
         &self%particle_array,self%sph_kernel,self%sim_params%sph_eta,self%sim_params%sim_min_d,.true._log_kind)
         
      ! set luminosities of sources outside sphere to zero
      do i=1,size(self%point_source_array)
         array_mask(i)=sum((self%point_source_array(i)%position-centre)**2)<=radius**2
         self%point_source_array(i)%position=self%point_source_array(i)%position-centre      
         if (.not.array_mask(i)) self%point_source_array(i)%luminosity=0._rel_kind
      end do
      
      
      
      return
   
   end subroutine


end module