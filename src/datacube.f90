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
   use m_constants_parameters
   use m_maths
   use m_ray
   use m_kernel
   use m_binary_tree
   use m_dust
   use m_source
   use omp_lib
   
   implicit none
   
   private
   public :: datacube
   
   ! define 3D space datacube class
   type :: datacube
   
      real(kind=rel_kind),allocatable :: x(:)               ! x coordinate
      real(kind=rel_kind),allocatable :: y(:)               ! y coordinate
      real(kind=rel_kind),allocatable :: lambda(:)          ! wavelengths
      real(kind=rel_kind),allocatable :: i_lambda(:,:,:)    ! intensity map (lambda,x,y)
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: psf_convolve
      procedure,non_overridable :: destroy
   
   end type
   
   contains
   
   subroutine initialise(self,sph_kernel,sph_tree,dust_prop,&
      &x_min,x_max,n_x,y_min,y_max,n_y,&
      &angle,lambda_array,v_ob_in,background,point_sources)
   
      ! argument declarations
      class(datacube),intent(inout) :: self                 ! datacube object
      class(kernel),intent(inout),target :: sph_kernel      ! SPH kernel object
      class(dust),intent(inout),target :: dust_prop         ! dust object
      type(binary_tree),intent(inout) :: sph_tree           ! SPH tree object
      real(kind=rel_kind),intent(in) :: x_min               ! minimum x position (relative to CoM)
      real(kind=rel_kind),intent(in) :: x_max               ! maximum x position (relative to CoM)
      integer(kind=int_kind),intent(in) :: n_x              ! number of x points
      real(kind=rel_kind),intent(in) :: y_min               ! minimum y position (relative to CoM)
      real(kind=rel_kind),intent(in) :: y_max               ! maximum y position (relative to CoM)
      integer(kind=int_kind),intent(in) :: n_y              ! number of y points
      real(kind=rel_kind),intent(in) :: angle(3)            ! (altitude,azimuth,rotation)
      real(kind=rel_kind),intent(in) :: lambda_array(:)     ! custom array of wavelengths
      real(kind=rel_kind),intent(in),optional :: v_ob_in(n_dim)   ! velocity of observer
      class(source),intent(in),optional :: background       ! background radiation field
      class(source),intent(in),optional :: point_sources(:) ! point sources
      
      ! variable declarations
      type(ray),allocatable :: sph_ray(:)                   ! ray object (1 per thread)
      integer(kind=int_kind) :: i,j,k,l                     ! counter
      integer(kind=int_kind) :: num_threads                 ! number of OpenMP threads
      integer(kind=int_kind) :: thread_num                  ! thread number
      real(kind=rel_kind) :: i_bg                           ! background intensity
      real(kind=rel_kind) :: img_x_unit(n_dim)              ! image x unit vector in 3D space (width)
      real(kind=rel_kind) :: img_y_unit(n_dim)              ! image y unit vector in 3D space (height)
      real(kind=rel_kind) :: img_z_unit(n_dim)              ! image z unit vector in 3D space (depth)
      real(kind=rel_kind) :: pix_position(n_dim)            ! pixel position in 3D space
      real(kind=rel_kind) :: v_ob(3)                        ! velocity of observer
      real(kind=rel_kind) :: x_img                          ! x position in image coordinates
      real(kind=rel_kind) :: y_img                          ! y position in image coordinates
      real(kind=rel_kind) :: v_rec                          ! recession velocity of source
      
      ! allocate arrays
      allocate(self%x(n_x))
      allocate(self%y(n_y))
      allocate(self%lambda(size(lambda_array)))
      allocate(self%i_lambda(size(self%lambda),size(self%x),size(self%y)))
      
      ! set up coordinates
      self%x=lin_space(x_min,x_max,n_x)
      self%y=lin_space(y_min,y_max,n_y)
      self%lambda=lambda_array

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
      img_y_unit=rotate_vector(img_y_unit,(/1._rel_kind,0._rel_kind,0._rel_kind/),angle(1))
      img_z_unit=rotate_vector(img_z_unit,(/1._rel_kind,0._rel_kind,0._rel_kind/),angle(1))
      
      ! transform 2: azimuth
      ! move y axis towards x axis
      img_x_unit=rotate_vector(img_x_unit,(/0._rel_kind,0._rel_kind,1._rel_kind/),angle(2))
      img_y_unit=rotate_vector(img_y_unit,(/0._rel_kind,0._rel_kind,1._rel_kind/),angle(2))
      img_z_unit=rotate_vector(img_z_unit,(/0._rel_kind,0._rel_kind,1._rel_kind/),angle(2))
      
      ! transform 3: rotation
      ! rotate image
      img_x_unit=rotate_vector(img_x_unit,img_z_unit,angle(3))
      img_y_unit=rotate_vector(img_y_unit,img_z_unit,angle(3))
      
      ! set recession velocity of source
      if (present(background)) then
         v_rec=dot_product(background%velocity-v_ob,-img_z_unit)
      else
         v_rec=0._rel_kind
      end if
      
      ! set up variables for multi-threading
      num_threads=omp_get_num_procs()
      allocate(sph_ray(num_threads))
      do i=1,size(sph_ray)
         call sph_ray(i)%initialise(sph_kernel,sph_tree,dust_prop)
      end do
      
      ! generate data cube from dust
      !$omp parallel num_threads(num_threads) default(shared) private(i,j,k,pix_position,i_bg,thread_num)
      
         thread_num=omp_get_thread_num()+1
         !$omp do collapse(1)
         
            do k=1,size(self%i_lambda,3)
            
               do j=1,size(self%i_lambda,2)
               
                  pix_position=sph_tree%com+sph_tree%max_length*img_z_unit+self%x(j)*img_x_unit+self%y(k)*img_y_unit
                  call sph_ray(thread_num)%ray_trace_initialise(pix_position,-1._rel_kind*img_z_unit)
               
                  do concurrent (i=1:size(self%i_lambda,1))
                     
                     ! set background intensity
                     if (present(background)) then
                        i_bg=background%intensity(self%lambda(i),v_rec)
                     else
                        i_bg=0._rel_kind
                     end if
                  
                     self%i_lambda(i,j,k)=sph_ray(thread_num)%ray_trace_i(self%lambda(i),v_ob,i_bg)
                  
                  end do
                 
               end do
               
            end do
         
         !$omp end do
      
      !$omp end parallel
      
      ! generate intensities from stars
      if (present(point_sources)) then
      
         !$omp parallel num_threads(num_threads) default(shared) private(i,j,k,l,x_img,y_img,v_rec,thread_num)
      
            thread_num=omp_get_thread_num()+1
            !$omp do
      
               do l=1,size(point_sources)
      
                  ! set up ray from star to screen
                  call sph_ray(thread_num)%ray_trace_initialise(point_sources(l)%position,img_z_unit)
         
                  x_img=dot_product(point_sources(l)%position-sph_tree%com,img_x_unit)
                  y_img=dot_product(point_sources(l)%position-sph_tree%com,img_y_unit)
         
                  if (x_img>self%x(1).and.x_img<self%x(size(self%x)).and.y_img>self%y(1).and.y_img<self%y(size(self%y))) then
         
                     j=binary_search(x_img,self%x)
                     k=binary_search(y_img,self%y)
                     v_rec=dot_product(point_sources(l)%velocity-v_ob,-img_z_unit)
               
                     do concurrent (i=1:size(self%i_lambda,1))
               
                        self%i_lambda(i,j,k)=self%i_lambda(i,j,k)+exp(-sph_ray(thread_num)%ray_trace_tau(self%lambda(i),v_ob))*&
                           &point_sources(l)%intensity(self%lambda(i),v_rec)*point_sources(l)%luminosity/&
                           &(point_sources(l)%bolometric_intensity*4._rel_kind*pi*(self%x(j+1)-self%x(j))*(self%y(k+1)-self%y(k)))

                     end do
         
                  end if
         
               end do
            !$omp end do
            
         !$omp end parallel
      
      end if
      
      call sph_ray%destroy()

      deallocate(sph_ray) 
       
      return
   
   end subroutine
   
   ! convolve datacube images with Gaussian PSF
   subroutine psf_convolve(self,distance,fwhm_array)

      include "fftw3.f03"

      ! argument declarations
      class(datacube),intent(inout) :: self                                ! datacube object
      real(kind=rel_kind) :: distance                                      ! observer source distance
      real(kind=rel_kind) :: fwhm_array(:)                                 ! array of beam fwhm (arcseconds)

      ! variable declarations
      integer(kind=int_kind) :: i,j,k                                      ! counter
      integer(kind=int_kind) :: plan                                       ! fftw plan variable
      real(kind=rel_kind) :: centre(2)                                     ! centre of grid
      real(kind=rel_kind) :: r                                             ! distance from centre of grid
      real(kind=rel_kind) :: sigma                                         ! standard deviation of absolute beam size
      real(kind=rel_kind) :: sigma_conv                                    ! fwhm -> sigma conversion factor
      real(kind=rel_kind) :: psf_norm                                      ! psf normalisation
      complex(kind=cpx_kind) :: intensity_map(size(self%x),size(self%y))   ! intensity map
      complex(kind=cpx_kind) :: psf_map(size(self%x),size(self%y))         ! point spread function

      ! check fwhm array has correct number of elements
      if (size(fwhm_array)/=size(self%lambda)) then
         write(*,"(A)") "Incorrect number of beam sizes. Cannot convolve."
         return      
      end if

      centre=0.5_rel_kind*(/self%x(1)+self%x(size(self%x)),self%y(1)+self%y(size(self%y))/)
      sigma_conv=distance*arcsec_rad/fwhm_sigma

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
   
   pure subroutine destroy(self)
   
      ! argument declarations
      class(datacube),intent(inout) :: self        ! datacube object
      
      deallocate(self%x)
      deallocate(self%y)
      deallocate(self%lambda)
      deallocate(self%i_lambda)
      
      return
   
   end subroutine

end module