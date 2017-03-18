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

! datacube module
module m_datacube

   use m_kind_parameters
   use m_options
   use m_constants_parameters
   use m_maths
   use m_string
   use m_ray
   use m_kernel
   use m_binary_tree
   use m_dust
   use m_source
   use m_image_tree_node
   use omp_lib
   use ieee_arithmetic
   
   implicit none
   
   private
   public :: datacube
   
   ! define 3D space datacube class
   type :: datacube
   
      type(options),pointer :: sim_params
      type(image_tree_node),pointer :: image_tree_root      ! root node of image tree
      type(proj_particle),pointer,contiguous :: image_tree_particle_array(:) ! projected particles for image tree
      real(kind=rel_kind),allocatable :: x(:)               ! x coordinate
      real(kind=rel_kind),allocatable :: y(:)               ! y coordinate
      real(kind=rel_kind),allocatable :: lambda(:)          ! wavelengths
      real(kind=rel_kind),allocatable :: t_bins(:)          ! temperature bins
      real(kind=rel_kind),allocatable :: i_lambda(:,:,:)    ! flux density
      real(kind=rel_kind),allocatable :: sigma(:,:)         ! column density
      real(kind=rel_kind),allocatable :: ps_l_lambda(:,:,:) ! point source luminosity (unattuenuated and attenuated)
      real(kind=rel_kind),allocatable :: bg_i_lambda(:)     ! background intensity
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: psf_convolve
      procedure,non_overridable :: destroy
      procedure,non_overridable :: write_out
   
   end type
   
   contains
   
   subroutine initialise(self,sim_params,background,point_sources,sph_kernel,sph_tree,dust_prop,v_ob_in)
   
      ! argument declarations
      class(datacube),intent(inout) :: self                 ! datacube object
      type(options),intent(in),target :: sim_params         ! simulation parameters object
      class(source),intent(in),optional :: background       ! background radiation field
      class(source),intent(in),optional :: point_sources(:) ! point sources
      class(kernel),intent(inout),target :: sph_kernel      ! SPH kernel object
      type(binary_tree),intent(inout) :: sph_tree           ! SPH tree object
      class(dust),intent(inout),target :: dust_prop         ! dust object
      real(kind=rel_kind),intent(in),optional :: v_ob_in(n_dim)   ! velocity of observer

      
      ! variable declarations
      type(ray),allocatable :: sph_ray(:)                   ! ray object (1 per thread)
      type(image_tree_node),pointer :: temp_node            ! temporary image tree node 
      integer(kind=int_kind) :: i,j,k                       ! counter
      integer(kind=int_kind) :: num_threads                 ! number of OpenMP threads
      integer(kind=int_kind) :: thread_num                  ! thread number
      integer(kind=int_kind) :: node_id                     ! node id variable
      integer(kind=int_kind) :: dc_max_level                ! datacube max level
      real(kind=rel_kind) :: dx                             ! x increment
      real(kind=rel_kind) :: dy                             ! y increment
      real(kind=rel_kind) :: img_x_unit(n_dim)              ! image x unit vector in 3D space (width)
      real(kind=rel_kind) :: img_y_unit(n_dim)              ! image y unit vector in 3D space (height)
      real(kind=rel_kind) :: img_z_unit(n_dim)              ! image z unit vector in 3D space (depth)
      real(kind=rel_kind) :: pix_position(n_dim)            ! pixel position in 3D space
      real(kind=rel_kind) :: v_ob(3)                        ! velocity of observer
      real(kind=rel_kind) :: x_img                          ! x position in image coordinates
      real(kind=rel_kind) :: y_img                          ! y position in image coordinates
      real(kind=rel_kind) :: v_rec                          ! recession velocity of source
      real(kind=rel_kind),allocatable :: aabb_array(:,:,:)  ! array of 2d axis aligned bounding boxes
      real(kind=rel_kind),allocatable :: i_lambda_array(:,:)! array of intensities
      real(kind=rel_kind),allocatable :: j_lambda_array(:,:)! dust emissivity array
      real(kind=rel_kind),allocatable :: t_array(:)        ! array of temperatures
      real(kind=rel_kind),allocatable :: w_array(:)         ! array of weights
      real(kind=rel_kind),allocatable :: sigma_array(:)     ! array of column densities
      real(kind=rel_kind),allocatable :: t_hist_array(:,:)  ! array of temperature histograms
      real(kind=rel_kind),allocatable :: t_moment_array(:,:)! array of temperature moments (element 1 is the mean)
      
      ! associated parameters object
      self%sim_params=>sim_params
      
      ! set number of pixels (nearest power of 2)
      dc_max_level=ceiling(log(real(max(self%sim_params%datacube_n_x,self%sim_params%datacube_n_y),rel_kind))/log(2._rel_kind))
            
      self%lambda=string_to_real(self%sim_params%datacube_lambda_string)
            
      ! allocate arrays
      allocate(self%x(2**dc_max_level))
      allocate(self%y(2**dc_max_level))
      allocate(self%i_lambda(size(self%lambda),size(self%x),size(self%y)))
      allocate(self%sigma(size(self%x),size(self%y)))
      allocate(self%image_tree_root)
      allocate(self%image_tree_particle_array(size(sph_tree%particle_array)))
      allocate(self%bg_i_lambda(size(self%lambda)))
      if (present(point_sources)) allocate(self%ps_l_lambda(2,size(self%lambda),size(point_sources)))
      
      
      ! set up coordinates
      dx=(self%sim_params%datacube_x_max-self%sim_params%datacube_x_min)/real(size(self%x)+1,rel_kind)
      dy=(self%sim_params%datacube_y_max-self%sim_params%datacube_y_min)/real(size(self%y)+1,rel_kind)
      self%x=lin_space(self%sim_params%datacube_x_min,self%sim_params%datacube_x_min+size(self%x)*dx,size(self%x))+0.5_rel_kind*dx
      self%y=lin_space(self%sim_params%datacube_y_min,self%sim_params%datacube_y_min+size(self%y)*dy,size(self%y))+0.5_rel_kind*dy
      if (self%sim_params%datacube_ppmap) self%t_bins=string_to_real(self%sim_params%datacube_t_string)

      if (present(v_ob_in)) then
         v_ob=v_ob_in
      else
         v_ob=0._rel_kind
      end if
      
      ! set up unit vectors
      img_x_unit=(/1._rel_kind,0._rel_kind,0._rel_kind/)
      img_y_unit=(/0._rel_kind,1._rel_kind,0._rel_kind/)
      img_z_unit=(/0._rel_kind,0._rel_kind,1._rel_kind/)
      
      ! transform 1: altitude
      ! move z axis towards y axis
      img_y_unit=rotate_vector(img_y_unit,(/1._rel_kind,0._rel_kind,0._rel_kind/),self%sim_params%datacube_angles(1))
      img_z_unit=rotate_vector(img_z_unit,(/1._rel_kind,0._rel_kind,0._rel_kind/),self%sim_params%datacube_angles(1))
      
      ! transform 2: azimuth
      ! move y axis towards x axis
      img_x_unit=rotate_vector(img_x_unit,(/0._rel_kind,0._rel_kind,1._rel_kind/),self%sim_params%datacube_angles(2))
      img_y_unit=rotate_vector(img_y_unit,(/0._rel_kind,0._rel_kind,1._rel_kind/),self%sim_params%datacube_angles(2))
      img_z_unit=rotate_vector(img_z_unit,(/0._rel_kind,0._rel_kind,1._rel_kind/),self%sim_params%datacube_angles(2))
      
      ! transform 3: rotation
      ! rotate image
      img_x_unit=rotate_vector(img_x_unit,img_z_unit,self%sim_params%datacube_angles(3))
      img_y_unit=rotate_vector(img_y_unit,img_z_unit,self%sim_params%datacube_angles(3))
      
      
      ! set up image tree particles
      do i=1,size(self%image_tree_particle_array)
         call self%image_tree_particle_array(i)%initialise(sph_tree%particle_array(i),&
            &reshape((/img_x_unit,img_y_unit/),(/n_dim,2/)))
      end do
      
      ! build image tree
      node_id=0
      call self%image_tree_root%initialise(self%image_tree_particle_array,reshape((/self%sim_params%datacube_x_min,&
         &self%sim_params%datacube_y_min,self%sim_params%datacube_x_max,self%sim_params%datacube_y_max/),(/2,2/)),node_id)
      allocate(aabb_array(2,2,self%image_tree_root%n_leaf))
      allocate(i_lambda_array(size(self%lambda),self%image_tree_root%n_leaf))
      allocate(j_lambda_array(size(self%lambda),self%image_tree_root%n_leaf))
      allocate(sigma_array(self%image_tree_root%n_leaf))
      aabb_array=self%image_tree_root%get_leaf_aabb()
      
      if (self%sim_params%datacube_ppmap) then
         allocate(t_hist_array(size(self%t_bins)-1,self%image_tree_root%n_leaf))
         allocate(t_moment_array(self%sim_params%datacube_n_t_moments,self%image_tree_root%n_leaf))
      end if

      ! set background intensity
      if (present(background)) then
         v_rec=dot_product(background%velocity-v_ob,-img_z_unit)
         self%bg_i_lambda=background%intensity(self%lambda,v_rec)
      else
         self%bg_i_lambda=0._rel_kind
      end if
      
      ! set up variables for multi-threading
      num_threads=omp_get_num_procs()
      allocate(sph_ray(num_threads))
      do i=1,size(sph_ray)
         call sph_ray(i)%initialise(sph_kernel,sph_tree,dust_prop)
      end do
      
      !$omp parallel num_threads(num_threads) default(shared) private(i,j,k,pix_position,thread_num,t_array,w_array)
         !$omp do schedule(dynamic)
            do j=1,self%image_tree_root%n_leaf
            
               thread_num=omp_get_thread_num()+1
            
               pix_position=sph_tree%com+0.5_rel_kind*sph_tree%max_length*img_z_unit+&
                  &0.5_rel_kind*(aabb_array(1,2,j)+aabb_array(1,1,j))*img_x_unit+&
                  &0.5_rel_kind*(aabb_array(2,2,j)+aabb_array(2,1,j))*img_y_unit
               call sph_ray(thread_num)%ray_trace_initialise(pix_position,-1._rel_kind*img_z_unit)
            
               do i=1,size(self%lambda)
            
                  i_lambda_array(i,j)=sph_ray(thread_num)%ray_trace_i(self%lambda(i),v_ob,self%bg_i_lambda(i))
                  j_lambda_array(i,j)=sph_ray(thread_num)%ray_trace_j(self%lambda(i),v_ob)+self%bg_i_lambda(i)
            
               end do
            
               ! set column density
               sigma_array(j)=sum(sph_ray(thread_num)%item(:sph_ray(thread_num)%n_item)%sigma*&
                  &sph_ray(thread_num)%item(:sph_ray(thread_num)%n_item)%f_sub)
                  
               if (self%sim_params%datacube_ppmap) then
               
                  t_array=dust_prop%dust_temperature(sph_ray(thread_num)%item(:sph_ray(thread_num)%n_item)%a_dot)
                  w_array=sph_ray(thread_num)%item(:sph_ray(thread_num)%n_item)%sigma*&
                     &sph_ray(thread_num)%item(:sph_ray(thread_num)%n_item)%f_sub
               
                  ! calculate temperature differential column density
                  if (sph_ray(thread_num)%n_item>0) then
                     t_hist_array(:,j)=make_hist(t_array,self%t_bins,w_array)
                  else
                     t_hist_array(:,j)=0._rel_kind
                  end if
                  
                  if (self%sim_params%datacube_log_moments) t_array=log10(t_array)
                     
                  ! calculate column weighted temperature mean
                  t_moment_array(1,j)=mean(t_array,w_array)
                  
                  ! calculate column weighted temperature central moments
                  do k=2,size(t_moment_array,1)
                     t_moment_array(k,j)=dist_moment(t_array,t_moment_array(1,j),k,w_array)
                  end do

               end if
                  
            end do
      
         !$omp end do
      !$omp end parallel
      
      ! copy i_lambda_array and sigma_array to tree
      call self%image_tree_root%set_leaf_i_lambda(i_lambda_array,j_lambda_array,sigma_array,t_hist_array,t_moment_array)     
      
      
      ! calculate point source luminosities
      if (present(point_sources)) then
         
         do j=1,size(self%ps_l_lambda,3)
         
            ! attenuate point source luminosity
            call sph_ray(1)%ray_trace_initialise(point_sources(j)%position,img_z_unit)
             v_rec=dot_product(point_sources(j)%velocity-v_ob,-img_z_unit)
         
            do i=1,size(self%ps_l_lambda,2)
            
               ! set point source luminosity
               self%ps_l_lambda(1,i,j)=point_sources(j)%luminosity*point_sources(j)%intensity(self%lambda(i),v_rec)&
                  &/point_sources(j)%bolometric_intensity
               
               ! set attenuated luminosity
               self%ps_l_lambda(2,i,j)=self%ps_l_lambda(1,i,j)*exp(-sph_ray(1)%ray_trace_tau(self%lambda(i),v_ob))
               
               ! add luminosity to intensity map
               x_img=dot_product(point_sources(j)%position-sph_tree%com,img_x_unit)
               y_img=dot_product(point_sources(j)%position-sph_tree%com,img_y_unit)
               
               ! find correct image_tree node
               temp_node=>self%image_tree_root%get_node_pointer((/x_img,y_img/))
               if (associated(temp_node)) then
                  temp_node%i_lambda(i)=temp_node%i_lambda(i)+self%ps_l_lambda(2,i,j)/&
                     &(4._rel_kind*pi*product(temp_node%aabb(:,2)-temp_node%aabb(:,1)))
                  temp_node%j_lambda(i)=temp_node%j_lambda(i)+self%ps_l_lambda(1,i,j)/&
                     &(4._rel_kind*pi*product(temp_node%aabb(:,2)-temp_node%aabb(:,1)))
               end if
               
            end do
         
         end do
         
         call self%image_tree_root%rebuild_tree()
         
      end if
      
      ! build datacube
      do j=1,size(self%y)
         do i=1,size(self%x)
         
            temp_node=>self%image_tree_root%get_node_pointer((/self%x(i),self%y(j)/),dc_max_level)
            self%i_lambda(:,i,j)=temp_node%i_lambda
            self%sigma(i,j)=temp_node%sigma
         
         end do
      end do

      call sph_ray%destroy()

      deallocate(sph_ray)
       
      return
   
   end subroutine
   
   ! convolve datacube images with Gaussian PSF
   subroutine psf_convolve(self)

      include "fftw3.f03"

      ! argument declarations
      class(datacube),intent(inout) :: self                                ! datacube object

      ! variable declarations
      integer(kind=int_kind) :: i,j,k                                      ! counter
      integer(kind=int_kind) :: plan                                       ! fftw plan variable
      real(kind=rel_kind) :: centre(2)                                     ! centre of grid
      real(kind=rel_kind) :: r                                             ! distance from centre of grid
      real(kind=rel_kind) :: sigma                                         ! standard deviation of absolute beam size
      real(kind=rel_kind) :: sigma_conv                                    ! fwhm -> sigma conversion factor
      real(kind=rel_kind) :: psf_norm                                      ! psf normalisation
      real(kind=rel_kind),allocatable :: fwhm_array(:)                     ! full width half maximum array
      complex(kind=cpx_kind) :: intensity_map(size(self%x),size(self%y))   ! intensity map
      complex(kind=cpx_kind) :: psf_map(size(self%x),size(self%y))         ! point spread function

      fwhm_array=string_to_real(self%sim_params%datacube_fwhm_string)

      ! check fwhm array has correct number of elements
      if (size(fwhm_array)/=size(self%lambda)) then
         write(*,"(A)") "Incorrect number of beam sizes. Cannot convolve."
         return      
      end if

      centre=0.5_rel_kind*(/self%x(1)+self%x(size(self%x)),self%y(1)+self%y(size(self%y))/)
      sigma_conv=self%sim_params%datacube_distance*arcsec_rad/fwhm_sigma

      do i=1,size(self%lambda)
 
         ! set complex map
         intensity_map=cmplx(self%i_lambda(i,:,:),0.,rel_kind)
         sigma=max(fwhm_array(i)*sigma_conv,self%x(2)-self%x(1),self%y(2)-self%y(1))
 
         ! set pdf grid
         do concurrent (k=1:size(self%y))
            do concurrent (j=1:size(self%x))
    
               ! calculate distance from centre
               r=norm2((/self%x(j),self%y(k)/)-centre)
       
               ! set psf
               psf_map(j,k)=cmplx(exp(-r**2/(2._rel_kind*sigma**2)),0.,rel_kind)
    
            end do
         end do
         ! normalise psf
         psf_norm=sum(real(psf_map,rel_kind))*real(size(psf_map))
         psf_map=cmplx(real(psf_map,rel_kind)/psf_norm,0.,rel_kind)
 
         ! cshift psf
         psf_map=cshift(psf_map,size(self%x)/2,1)
         psf_map=cshift(psf_map,size(self%y)/2,2)
 
 
         ! Fourier transform intensity map
         call dfftw_plan_dft_2d(plan,size(self%y),size(self%x),intensity_map,intensity_map,fftw_forward,fftw_estimate)
         call dfftw_execute_dft(plan,intensity_map,intensity_map)
 
         ! Fourier transform psf
         call dfftw_plan_dft_2d(plan,size(self%y),size(self%x),psf_map,psf_map,fftw_forward,fftw_estimate)
         call dfftw_execute_dft(plan,psf_map,psf_map)
 
         ! Inverse Fourier transform product
         intensity_map=intensity_map*psf_map
         call dfftw_plan_dft_2d(plan,size(self%y),size(self%x),intensity_map,intensity_map,fftw_backward,fftw_estimate)
         call dfftw_execute_dft(plan,intensity_map,intensity_map)
 
         ! set datacube map to smoothed version
         self%i_lambda(i,:,:)=real(intensity_map,rel_kind)

      end do

      return

   end subroutine
   
   subroutine write_out(self)
   
      ! argument declarations
      class(datacube),intent(in) :: self                                ! datacube object
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                                     ! counter
      integer(kind=int_kind) :: n_ps                                    ! number of point sources
      character(kind=chr_kind,len=string_length) :: file_name           ! output file name
      character(kind=chr_kind,len=string_length) :: format_string       ! format string
      
      ! write out leaf cell mosaic
      write(file_name,"(A)") trim(self%sim_params%sim_id)//"/mosaic.dat"
      open(1,file=trim(file_name))
      
      j=0
      write(1,"(A)") "# 1 x_min (cm)"
      write(1,"(A)") "# 2 y_min (cm)"
      write(1,"(A)") "# 3 x_max (cm)"
      write(1,"(A)") "# 4 y_max (cm)"
      write(1,"(A)") "# 5 sigma (g cm^-2)"
      if (self%sim_params%datacube_ppmap) then
         do i=1,size(self%t_bins)-1
            j=j+1
            write(1,"(A,I0,A,F6.2,A,F6.2,A)") "# ",5+j," sigma_T [",self%t_bins(i),"K < T <= ",self%t_bins(i+1),"] (g cm^-2 K^-1)"
         end do
         j=j+1
         write(1,"(A,I0,A)") "# ",5+j," T mean (K)"
         do i=2,self%sim_params%datacube_n_t_moments
            j=j+1
            write(1,"(A,I0,A,I0,A,I0,A)") "# ",5+j," T central moment ",i," (K)"
         end do
      end if
      do i=1,size(self%lambda)
         write(1,"(A,I0,A,E10.4,A)") "# ",5+i+j," I_lambda=",self%lambda(i)," micron (erg s^-1 sr^-1 cm^-2 micron^-1)"
      end do
      do i=1,size(self%lambda)
         write(1,"(A,I0,A,E10.4,A)") "# ",5+i+j+size(self%lambda),&
            &" J_lambda=",self%lambda(i)," micron (erg s^-1 sr^-1 cm^-2 micron^-1)"
      end do
      write(1,*)
      call self%image_tree_root%write_leaf(1)
      
      close(1)
      
      ! write out datacube
      write(file_name,"(A)") trim(self%sim_params%sim_id)//"/datacube.dat"
      open(1,file=trim(file_name))
      write(format_string,"(A,I0,A)") "(",3+size(self%i_lambda),"(E25.17))"
      
      write(1,"(A)") "# 1 x (cm)"
      write(1,"(A)") "# 2 y (cm)"
      write(1,"(A)") "# 3 sigma (g cm^-3)"
      do i=1,size(self%lambda)
         write(1,"(A,I0,A,E10.4,A)") "# ",i+3," I_lambda=",self%lambda(i)," micron (erg s^-1 sr^-1 cm^-2 micron^-1)"
      end do
      write(1,*)
      do j=1,size(self%y)
         do i=1,size(self%x)
            write(1,trim(format_string)) self%x(i),self%y(j),self%sigma(i,j),self%i_lambda(:,i,j)
         end do
      end do
      
      close(1)
      
      ! write out spectrum
      write(file_name,"(A)") trim(self%sim_params%sim_id)//"/spectrum.dat"
      open(1,file=trim(file_name))
      if (allocated(self%ps_l_lambda)) then
         n_ps=size(self%ps_l_lambda,3)
      else
         n_ps=0
      end if
      write(format_string,"(A,I0,A)") "(",4+2*(n_ps+1),"(E25.17))"
      
      write(1,"(A)") "# 1 lambda (micron)"
      write(1,"(A)") "# 2 F_lambda 4 pi D^2 [full] (erg s^-1 micron^-1)"
      write(1,"(A)") "# 3 F_lambda 4 pi D^2 [full] (erg s^-1 micron^-1) [no extinction]"
      write(1,"(A)") "# 4 F_lambda 4 pi D^2 [background] (erg s^-1 micron^-1)"
      if (n_ps>0) then
         write(1,"(A)") "# 5 F_lambda 4 pi D^2 [total point sources] (erg s^-1 micron^-1) [no extinction]"
         write(1,"(A)") "# 6 F_lamdba 4 pi D^2 [total point sources] (erg s^-1 micron^-1)"
         do i=1,n_ps
            write(1,"(A,I0,A,I0,A)") "# ",5+2*i," F_lambda 4 pi D^2 [point source ",i,"] (erg s^-1 micron^-1) [no extinction]"
            write(1,"(A,I0,A,I0,A)") "# ",6+2*i," F_lambda 4 pi D^2 [point source ",i,"] (erg s^-1 micron^-1)"
         end do
      end if
      write(1,*)
      do i=1,size(self%lambda)
         if (n_ps>0) then
            write(1,trim(format_string)) self%lambda(i),self%image_tree_root%i_lambda(i)*4._rel_kind*pi*&
               &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1)),&
               &self%image_tree_root%j_lambda(i)*4._rel_kind*pi*&
               &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1)),&
               &self%bg_i_lambda(i)*4._rel_kind*pi*&
               &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1)),&
               &sum(self%ps_l_lambda(1,i,:)),sum(self%ps_l_lambda(2,i,:)),&
               &pack(self%ps_l_lambda(:,i,:),.true.)
         else
            write(1,trim(format_string)) self%lambda(i),self%image_tree_root%i_lambda(i)*4._rel_kind*pi*&
               &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1)),&
               &self%image_tree_root%j_lambda(i)*4._rel_kind*pi*&
               &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1)),&
               &self%bg_i_lambda(i)*4._rel_kind*pi*&
               &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1))
         end if
      end do
         
      close(1)
      
      write(*,"(A,E25.17)") "apparent luminosity: ",trapz_intgr(self%lambda,self%image_tree_root%i_lambda)*4._rel_kind*pi*&
         &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1))                 
      write(*,"(A,E25.17)") "emission luminosity: ",trapz_intgr(self%lambda,self%image_tree_root%j_lambda)*4._rel_kind*pi*&
         &product(self%image_tree_root%aabb(:,2)-self%image_tree_root%aabb(:,1))
      
      return
   
   end subroutine
   
   pure subroutine destroy(self)
   
      ! argument declarations
      class(datacube),intent(inout) :: self        ! datacube object
      
      deallocate(self%x)
      deallocate(self%y)
      deallocate(self%lambda)
      deallocate(self%i_lambda)
      deallocate(self%sigma)
      deallocate(self%bg_i_lambda)
      if (allocated(self%ps_l_lambda)) deallocate(self%ps_l_lambda)
      call self%image_tree_root%destroy()
      
      return
   
   end subroutine

end module