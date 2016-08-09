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

! ray module for luminosity packet propagation
module m_ray

   use m_kind_parameters
   use m_sph_parameters
   use m_constants_parameters
   use m_maths
   use m_line
   use m_particle
   use m_kernel
   use m_binary_tree
   use m_binary_tree_node
   use m_dust
   use m_atomic_update
   
   implicit none
   
   private
   public :: ray
   
   ! ray particle type
   type :: ray_particle
      
      real(kind=rel_kind) :: r(n_dim)           ! particle position
      real(kind=rel_kind) :: v(n_dim)           ! paticle velocity
      real(kind=rel_kind) :: m                  ! particle mass
      real(kind=rel_kind) :: inv_h              ! one over smoothing length
      real(kind=rel_kind) :: inv_rho            ! one over density
      real(kind=rel_kind) :: a_dot                  ! temperature of particle
      real(kind=rel_kind) :: t_r                ! particle postion along ray
      real(kind=rel_kind) :: b                  ! impact m_parameter
      real(kind=rel_kind) :: s_0                ! paticle-origin distance
      real(kind=rel_kind) :: s_1                ! particle-terminal distance
      real(kind=rel_kind) :: sigma              ! column density through particle
      real(kind=rel_kind) :: lambda_ob          ! rest frame wavelength
      real(kind=rel_kind) :: kappa_ext          ! dust mass extinction coeff
      real(kind=rel_kind) :: a                  ! dust albedo
      real(kind=rel_kind) :: f_sub              ! dust sublimation fraction
      type(particle),pointer :: particle_ptr    ! particle pointer
      
      contains
      
      procedure,non_overridable :: initialise=>initialise_ray_particle
      
   end type
   
   ! ray type
   type :: ray
   
      integer(kind=int_kind) :: n_item          ! number of items in item array
      real(kind=rel_kind) :: l_chunk            ! chunk of luminosity
      real(kind=rel_kind) :: v_em(n_dim)        ! emission velocity
      real(kind=rel_kind) :: lambda_em          ! wavelength of emission
      real(kind=rel_kind) :: min_length         ! minimum ray length
      real(kind=rel_kind) :: max_length         ! maximum ray length
      type(line) :: path                        ! path of ray
      type(ray_particle),allocatable :: item(:) ! ray particle array
      type(binary_tree_node),pointer :: root_node  ! root node of binary tree
      class(kernel),pointer :: sph_kernel       ! sph kernel object
      class(dust),pointer :: dust_prop          ! dust properties object
      
      contains
      
      procedure,non_overridable :: initialise=>initialise_ray
      procedure,non_overridable :: destroy=>destroy_ray
      procedure,non_overridable :: follow
      procedure,non_overridable :: ray_trace_initialise
      procedure,non_overridable :: ray_trace_i
      procedure,non_overridable :: ray_trace_tau
      procedure,non_overridable,private :: v_sph
      procedure,non_overridable,private :: a_dot_sph
      procedure,non_overridable,private :: a_sph
      procedure,non_overridable,private :: tau_sph
      procedure,non_overridable,private :: inv_mfp
      procedure,non_overridable,private :: grad_inv_mfp
      procedure,non_overridable,private :: find_particles
      procedure,non_overridable,private :: stock_geometric
      procedure,non_overridable,private :: stock_optical
      procedure,non_overridable,private :: stock_sigma
      procedure,non_overridable,private :: update_absorption_rate
      procedure,non_overridable,private :: mrw_check
      procedure,non_overridable,private :: h_0

   end type
   
   contains
   
   ! associate particle from particle_array with ray_particle object
   pure subroutine initialise_ray_particle(self,sph_particle)
   
      ! argument declarations
      class(ray_particle),intent(inout) :: self ! ray particle object
      type(particle),intent(inout),target :: sph_particle       ! sph particle
      
      ! set quantities
      self%r=sph_particle%r
      self%v=sph_particle%v
      self%m=sph_particle%m
      self%inv_h=sph_particle%inv_h
      self%inv_rho=sph_particle%inv_rho
      self%a_dot=sph_particle%a_dot
      self%f_sub=sph_particle%f_sub
      self%particle_ptr=>sph_particle
      
      return
   
   end subroutine
   
   ! initialise the ray
   pure subroutine initialise_ray(self,sph_kernel,sph_tree,dust_prop)
   
      ! argument declarations
      class(ray),intent(inout) :: self                          ! ray object
      class(kernel),intent(inout),target :: sph_kernel          ! sph kernel object
      class(dust),intent(inout),target :: dust_prop             ! dust properties object
      type(binary_tree),intent(inout) :: sph_tree               ! particle tree
      
      ! associate kernel and tree root node
      self%dust_prop=>dust_prop
      self%sph_kernel=>sph_kernel
      self%root_node=>sph_tree%root_node

      self%min_length=sph_tree%min_length
      self%max_length=2._rel_kind*sph_tree%max_length
      
      ! allocate ray_particle array
      allocate(self%item(n_ray))
      
      return
   
   end subroutine

   ! deallocate ray object
   elemental subroutine destroy_ray(self)   
   
      ! argument declarations
      class(ray),intent(inout) :: self          ! ray object
      
      deallocate(self%item)
      
      return
   
   end subroutine
   
   ! follow a luminosity packet
   subroutine follow(self,origin_in,direction_in,v_em_in,lambda_em_in,l_chunk,use_mrw,mrw_gamma)
   
      ! argument declarations
      class(ray),intent(inout) :: self                  ! ray object
      real(kind=rel_kind),intent(in) :: origin_in(n_dim)        ! origin of ray
      real(kind=rel_kind),intent(in) :: direction_in(n_dim)     ! direction of ray
      real(kind=rel_kind),intent(in) :: v_em_in(n_dim)  ! velocity of emitter
      real(kind=rel_kind),intent(in) :: lambda_em_in    ! emission wavelength
      real(kind=rel_kind),intent(in) :: l_chunk         ! chunk of luminosity
      logical(kind=log_kind),intent(in) :: use_mrw      ! use modified random walk
      real(kind=rel_kind),intent(in) :: mrw_gamma       ! mrw gamma parameter
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                     ! counter
      real(kind=rel_kind) :: origin(n_dim)              ! origin of ray
      real(kind=rel_kind) :: direction(n_dim)           ! direction of ray
      real(kind=rel_kind) :: new_direction(n_dim)       ! new direction of ray
      real(kind=rel_kind) :: v_em(n_dim)                ! velocity of emitter
      real(kind=rel_kind) :: v_ob(n_dim)                ! velocity of observer
      real(kind=rel_kind) :: lambda_em                  ! emission wavelength
      real(kind=rel_kind) :: length                     ! length of ray
      real(kind=rel_kind) :: tau                        ! number of optical depths to travel
      real(kind=rel_kind) :: tau_est                    ! estimated value of tau
      real(kind=rel_kind) :: inv_mfp                    ! inverse mean free path (density * kappa_ext)
      real(kind=rel_kind) :: grad_inv_mfp(n_dim)        ! gradient of inverse mean free path
      real(kind=rel_kind) :: l_0                        ! iteration points 
      real(kind=rel_kind) :: l_1
      real(kind=rel_kind) :: l_b                        ! bracketing point
      real(kind=rel_kind) :: tau_0                      ! iteration poitns
      real(kind=rel_kind) :: tau_1
      real(kind=rel_kind) :: tau_b                      ! backeting point
      real(kind=rel_kind) :: dtau_dl                    ! iteration derivative
      real(kind=rel_kind) :: d2tau_dl2                  ! iteration second derivative
      real(kind=rel_kind) :: uni_r                      ! uniform random number
      real(kind=rel_kind) :: lambda_ob                  ! rest frame wavelength
      real(kind=rel_kind) :: h                          ! smoothing length
      logical(kind=log_kind) :: extend_ray              ! are we extending a ray?
      logical(kind=log_kind) :: short_ray               ! is ray shorter than sim resolution?
      real(kind=rel_kind) :: ave_a_dot                  ! ray-averaged a_dot
      real(kind=rel_kind) :: ave_rho                    ! ray-averaged density
      real(kind=rel_kind) :: mrw_distance               ! distance travelled along modified random walk
      real(kind=rel_kind) :: planck_abs                 ! planck mean absorption
      real(kind=rel_kind) :: planck_sca                 ! planck mean scatter
      
      ! set input variables
      origin=origin_in
      direction=direction_in
      v_em=v_em_in
      lambda_em=lambda_em_in
      
      
      self%l_chunk=l_chunk
      self%lambda_em=lambda_em
      self%v_em=v_em
      
      extend_ray=.false.
      
      ! calculate inv_mfp and grad_inv_mfp from tree
      inv_mfp=0._rel_kind
      grad_inv_mfp=0._rel_kind
      h=0._rel_kind
      call self%root_node%sph_inv_mfp(self%sph_kernel,self%dust_prop,&
         &origin,lambda_em,direction,v_em,inv_mfp,grad_inv_mfp,h)
         
      call random_number(tau)
      tau=-log(tau)
         
      ! allocate item list
      if (.not.allocated(self%item)) allocate(self%item(n_ray))
      
      do
      
         ! estimate length from local quantities
         length=min(2._rel_kind*ray_eta*tau/&
            &(inv_mfp+sqrt(max(inv_mfp**2+2._rel_kind*tau*dot_product(grad_inv_mfp,direction),0._rel_kind))),self%max_length)

         ! adjust length if it's too short
         if (length<max(self%min_length,h)) then
            short_ray=.true.
            length=max(length,self%min_length,h)
         else
            short_ray=.false.
         end if 
         
         ! initialise path
         call self%path%initialise(origin,direction,length)
         
         ! set number of items to 0
         self%n_item=0
         
         ! find particles
         call self%find_particles(self%root_node)
         
         ! photon trajectory is complete when n_item=0
         if (self%n_item==0) exit
         
         ! calculate quantities
         call self%stock_geometric()
         call self%stock_sigma()
         call self%stock_optical()
         
         ! check to see if we need to perform modified random walk
         
         if (short_ray.and.use_mrw.and.(.not.extend_ray)) then
         
            ! check if mean free path is less then gamma * h
            if (self%mrw_check(mrw_gamma)) then
            
            ! perform modified random walk
            ave_a_dot=sum(self%item(:self%n_item)%a_dot*self%item(:self%n_item)%sigma)/sum(self%item(:self%n_item)%sigma)
            ave_rho=sum(self%item(:self%n_item)%sigma)/self%path%length
            
            ! calculate mrw variables
            mrw_distance=self%dust_prop%mrw_distance(ave_a_dot,self%path%length,ave_rho)
            planck_abs=self%dust_prop%mrw_planck_abs(ave_a_dot)
            
            ! update particle absorption rates
            do i=1,self%n_item
            
               call atomic_real_add(self%item(i)%particle_ptr%a_dot,&
                  mrw_distance*self%item(i)%sigma*planck_abs/self%path%length)
            
               ! check if scattered light array exists
               if (associated(self%item(i)%particle_ptr%a_dot_scatter_array)) then
               
                  do j=1,size(self%item(i)%particle_ptr%a_dot_scatter_array)
                  
                     planck_sca=self%dust_prop%mrw_planck_sca(ave_a_dot,&
                        &self%item(i)%particle_ptr%lambda_array(j),self%item(i)%particle_ptr%lambda_array(j+1))
                        
                     call atomic_real_add(self%item(i)%particle_ptr%a_dot_scatter_array(j),&
                        &mrw_distance*self%item(i)%sigma*planck_sca/self%path%length)
                  
                  end do
               
               end if
            
            end do
            
            ! set new position, direction, wavelength
            origin=origin+self%path%direction*self%path%length
            direction=self%dust_prop%mrw_random_direction(direction)
            self%lambda_em=self%dust_prop%mrw_random_wavelength(ave_a_dot)
            self%v_em=self%v_sph()
            
            ! restock optical properties
            call self%stock_optical()
            inv_mfp=self%inv_mfp()
            grad_inv_mfp=self%grad_inv_mfp()
            h=self%h_0()
            
            ! skip to next loop iteration
            cycle
         
            end if
         
         end if 
         
         tau_est=self%tau_sph()
            
         ! was initial length estimate too short?
         if (tau_est<tau) then
         
            ! set new properties
            ! new optical depth and origin
            tau=tau-tau_est
            origin=origin+direction*length
            call self%update_absorption_rate()
            extend_ray=.true.
         
         else
         
            extend_ray=.false.
         
            ! perform Newton Raphson iterations until length has converged
            
            ! bracket root (tau=0) between 0 and tau_est
            tau_b=-tau
            tau_0=tau_est-tau
            l_b=0._rel_kind
            l_0=length
            
            do
            
               ! calculate derivatives
               dtau_dl=self%inv_mfp()
               d2tau_dl2=2._rel_kind*((tau_b-tau_0)/(l_b-l_0)**2-dtau_dl/(l_b-l_0))
            
               ! set new value of l
               l_1=l_0-2._rel_kind*tau_0/(dtau_dl+sqrt(dtau_dl**2-2._rel_kind*tau_0*d2tau_dl2))
               
               ! otherwise, set up next iteration
               length=l_1
               call self%path%initialise(origin,direction,length)
               call self%stock_sigma()
               tau_1=self%tau_sph()-tau
!                write(*,*) "iteration"
               
               
               ! exit if l_0 and l_1 have converged
               if (abs(l_0-l_1)/l_1<ray_conv) exit
               
               ! make sure tau_0 and tau_b are either side of solution
               
               if (tau_0*tau_1<0._rel_kind) then
                  ! tau_0 and tau_1 are either side of root (tau=0)
                  l_b=l_0
                  tau_b=tau_0
               end if
               
               l_0=l_1
               tau_0=tau_1
              
            end do
            
            call self%update_absorption_rate()
            
            ! set new properties
            
            ! set local velocity frame wavelength
            v_ob=self%v_sph()
            lambda_ob=self%lambda_em*&
               &(dot_product(v_ob-self%v_em,self%path%direction)/c_light+1._rel_kind)
               
            ! decide whether packet is absorbed or scattered
            call random_number(uni_r)
            if (uni_r<self%dust_prop%albedo(lambda_ob)) then
               
               ! packet is scattered
               new_direction=self%dust_prop%random_scatter_direction(direction,lambda_ob)
               lambda_em=lambda_ob
               
            
            else
            
               ! packet is absorbed and re-emitted
               new_direction=self%dust_prop%random_emission_direction()
               lambda_em=self%dust_prop%random_wavelength(self%a_dot_sph())
            
            end if
            
            ! set new origin
            origin=origin+direction*length
            direction=new_direction
            ! set new emission velocity
            v_em=v_ob
            
            call random_number(tau)
            tau=-log(tau)
            
            self%lambda_em=lambda_em
            self%v_em=v_em
            
            call self%stock_optical()

         end if
         
         ! set new inv_mfp and grad_inv_mfp
         inv_mfp=self%inv_mfp()
         grad_inv_mfp=self%grad_inv_mfp()
         h=self%h_0()
         
!          write(*,*) "wavelength", self%lambda_em
!          write(*,*) "items on ray",self%n_item
!          write(*,*) "live items",count(self%item(:self%n_item)%b<2._rel_kind)
!          write(*,*) "ray origin",self%path%origin
!          write(*,*) "ray terminus",origin
!          write(*,*) "ray length",self%path%length
!          write(*,*) "optical depth",tau
!          write(*,*)
!          
!          write(1,*) self%path%origin
!          write(1,*) origin
!          write(1,*)
     
      end do
!       write(1,*)
!       write(1,*)
     
      return
      
   end subroutine
   
   ! initialise a ray trace (self%initialise has already been called)
   pure subroutine ray_trace_initialise(self,position,direction)
   
      ! argument declarations
      class(ray),intent(inout) :: self                   ! ray object
      real(kind=rel_kind),intent(in) :: position(n_dim)  ! position of ray trace
      real(kind=rel_kind),intent(in) :: direction(n_dim) ! direction of ray trace (away from position)
      
      ! set almost infinitely long path
      call self%path%initialise(position,direction,self%max_length)
      
      ! find particles along path
      self%n_item=0
      call self%find_particles(self%root_node)
      
      ! calculated geometric quantities and col density
      call self%stock_geometric()
      call self%stock_sigma()
      
      ! sort particles along ray
      call ray_particle_sort(self%item(:self%n_item))
           
      return
   
   end subroutine
   
   ! calculate intensity from ray trace (ray_trace_initialise has already been called)
   pure function ray_trace_i(self,lambda_ob,v_ob,i_bg) result (i_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                      ! ray object      
      real(kind=rel_kind),intent(in) :: lambda_ob        ! observer rest frame wavelength
      real(kind=rel_kind),intent(in) :: v_ob(n_dim)      ! observer rest frame velocity
      real(kind=rel_kind),intent(in) :: i_bg             ! background intensity
      
      ! result declaration
      real(kind=rel_kind) :: i_value                     ! monochromatic intensity at position
      
      ! variable declaration
      integer(kind=int_kind) :: i                        ! counter
      real(kind=rel_kind) :: lambda_em                   ! wavelength in emission rest frame
      real(kind=rel_kind) :: j_value                     ! mass emissivity
      real(kind=rel_kind) :: kappa_ext                   ! extinction coeff
      real(kind=rel_kind) :: a                           ! albedo
      
      
      ! loop backwards over ray and calculate j_value
      i_value=i_bg
      do i=self%n_item,1,-1

         ! observed wavelength in emission rest frame
         lambda_em=lambda_ob*(dot_product(v_ob-self%item(i)%v,self%path%direction)/c_light+1._rel_kind)
         
         ! calculate opacity and albedo
         call self%dust_prop%dust_mass_ext_and_albedo(lambda_em,kappa_ext,a)
         
         ! sink term
         i_value=i_value*exp(-kappa_ext*self%item(i)%sigma)
         
         ! emission component source term
         j_value=self%dust_prop%mono_mass_emissivity(lambda_em,self%item(i)%a_dot)
         
         ! scattering source term
         if (associated(self%item(i)%particle_ptr%lambda_array)) then
         
            j_value=j_value+self%item(i)%particle_ptr%a_dot_scatter(lambda_em)*0.25_rel_kind/pi
         
         end if
         
         ! add source term
         i_value=i_value+j_value*(1._rel_kind-exp(-kappa_ext*self%item(i)%sigma))/kappa_ext
         
      end do
      
      return
   
   end function
   
   ! calculate optical depth from ray trace (ray_trace_initialise has already been called)
   pure function ray_trace_tau(self,lambda_ob,v_ob) result (tau_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                      ! ray object      
      real(kind=rel_kind),intent(in) :: lambda_ob        ! observer rest frame wavelength
      real(kind=rel_kind),intent(in) :: v_ob(n_dim)      ! observer rest frame velocity
      
      ! result declaration
      real(kind=rel_kind) :: tau_value                   ! optical depth from infinity to position
      
      ! variable declaration
      integer(kind=int_kind) :: i                        ! counter
      real(kind=rel_kind) :: lambda_em                   ! wavelength in emission rest frame
      
      ! tau value
      tau_value=0._rel_kind
      do i=1,self%n_item
         lambda_em=lambda_ob*(dot_product(v_ob-self%item(i)%v,self%path%direction)/c_light+1._rel_kind)
         tau_value=tau_value+self%item(i)%sigma*self%dust_prop%dust_mass_ext(lambda_em)
      end do
   
      return
   
   end function
   
   ! find all particles along ray
   pure recursive subroutine find_particles(self,tree_node)
   
      ! argument declarations
      class(ray),intent(inout) :: self                  ! ray object
      type(binary_tree_node),intent(inout) :: tree_node ! binary tree node object
      
      ! variable declarations
      type(ray_particle),allocatable :: temp_item(:)    ! temporary ray particle array
      integer(kind=int_kind) :: i                       ! counter
      logical(kind=log_kind),allocatable :: live_particles(:)   ! mask for "live" particles
      
      if (associated(tree_node%left_node)) then
      
         ! walk tree if not at leaf node
         if (tree_node%left_node%ray_aabb_intersection(self%path))&
            &call self%find_particles(tree_node%left_node)
         if (tree_node%right_node%ray_aabb_intersection(self%path))&
            &call self%find_particles(tree_node%right_node)
            
      else
      
         ! at leaf node.
         ! build live_particle array
         allocate(live_particles(size(tree_node%particle_array)))
         do concurrent (i=1:size(tree_node%particle_array))
         
            live_particles(i)=&
               &(sum((self%path%origin+&
               &dot_product(tree_node%particle_array(i)%r-self%path%origin,self%path%direction)*&
               &self%path%direction-&
               &tree_node%particle_array(i)%r)**2)*&
               &tree_node%particle_array(i)%inv_h**2)<self%sph_kernel%r_support**2
         
         end do
         
         ! first check if there's enough storage in item array, resize if necessary
         if (self%n_item+count(live_particles)>size(self%item)) then
         
            allocate(temp_item(2*size(self%item)))
            temp_item(:size(self%item))=self%item
            call move_alloc(temp_item,self%item)
         
         end if
         
         ! add live particles to ray
         do i=1,size(tree_node%particle_array)
         
            if (live_particles(i)) then
            
               self%n_item=self%n_item+1
               call self%item(self%n_item)%initialise(tree_node%particle_array(i))
            
            end if
         
         end do
         deallocate(live_particles)
            
      end if
      
      return
   
   end subroutine
   
   ! calculate geometric properties
   pure subroutine stock_geometric(self)
   
      ! argument declarations
      class(ray),intent(inout) :: self                  ! ray object
      
      ! variable declarations
      integer(kind=int_kind) :: i                       ! counter
      
      do concurrent (i=1:self%n_item)
      
         ! calculate particle position along ray
         self%item(i)%t_r=&
            &dot_product(self%item(i)%r-self%path%origin,self%path%direction)
            
         ! calculate impact parameter
         self%item(i)%b=&
            &norm2(self%path%origin+&
            &self%item(i)%t_r*self%path%direction-&
            &self%item(i)%r)*&
            &self%item(i)%inv_h
            
         ! calculate particle-origin distance
         self%item(i)%s_0=&
            &norm2(self%path%origin-&
            &self%item(i)%r)*&
            &self%item(i)%inv_h
            
         ! calculate particle-terminus distance
         self%item(i)%s_1=&
            &norm2(self%path%terminus-&
            &self%item(i)%r)*&
            &self%item(i)%inv_h
         
      end do
      
      return
   
   end subroutine
   
   ! calculate optical properties
   pure subroutine stock_optical(self)
   
      ! argument declarations
      class(ray),intent(inout) :: self                  ! ray object
      
      ! variable declarations
      integer(kind=int_kind) :: i                       ! counter
      
      do concurrent (i=1:self%n_item)
         
         ! calculate dust mass extinction and albedo
         self%item(i)%lambda_ob=self%lambda_em*&
            &(dot_product(self%item(i)%v-self%v_em,self%path%direction)/c_light+1._rel_kind)
         call self%dust_prop%dust_mass_ext_and_albedo(self%item(i)%lambda_ob,self%item(i)%kappa_ext,self%item(i)%a)
         
         ! adjust extinction by f_sub
         self%item(i)%kappa_ext=self%item(i)%kappa_ext*self%item(i)%f_sub
         
      end do
      
      return
   
   end subroutine
   
   ! calculate quantities along populated ray
   ! similar to stock, but calculates fewer quantities
   pure subroutine stock_sigma(self)
   
      ! argument declarations
      class(ray),intent(inout) :: self                  ! ray object
      
      ! variable declarations
      integer(kind=int_kind) :: i                       ! counter
      real(kind=rel_kind) :: sigma_0                    ! column density 0
      real(kind=rel_kind) :: sigma_1                    ! column density 1
         
      do concurrent (i=1:self%n_item)
            
         ! calculate particle-terminus distance
         self%item(i)%s_1=&
            &norm2(self%path%terminus-&
            &self%item(i)%r)*&
            &self%item(i)%inv_h
      
         ! calculate column density
         sigma_0=self%sph_kernel%sigma(self%item(i)%b,self%item(i)%s_0)
         sigma_1=self%sph_kernel%sigma(self%item(i)%b,self%item(i)%s_1)
         self%item(i)%sigma=self%item(i)%m*self%item(i)%inv_h**(n_dim-1)*&
            &merge(sigma_0+sigma_1,abs(sigma_0-sigma_1),&
            &self%item(i)%t_r>0._rel_kind.and.self%item(i)%t_r<self%path%length)
      
      end do
      
      return
   
   end subroutine
   
   ! update particle energy absorption rates
   subroutine update_absorption_rate(self)
   
      ! argument declarations
      class(ray),intent(inout) :: self                  ! ray object
      
      ! variable declarations
      integer(kind=int_kind) :: i                       ! counter
      integer(kind=int_kind) :: j                       ! index
      
      do i=1,self%n_item
      
         call atomic_real_add(self%item(i)%particle_ptr%a_dot_new,&
            &self%item(i)%sigma*self%item(i)%kappa_ext*(1._rel_kind-self%item(i)%a)*self%l_chunk)
            
         ! update scattered light histogram, if allocated.
         if (associated(self%item(i)%particle_ptr%lambda_array)) then
         
            if (self%item(i)%lambda_ob>=self%item(i)%particle_ptr%lambda_array(1).and.&
               &self%item(i)%lambda_ob<=self%item(i)%particle_ptr%lambda_array(&
               &size(self%item(i)%particle_ptr%lambda_array))) then
               
               j=binary_search(self%item(i)%lambda_ob,self%item(i)%particle_ptr%lambda_array)
               call atomic_real_add(self%item(i)%particle_ptr%a_dot_scatter_array(j),&
                  &self%item(i)%sigma*self%item(i)%kappa_ext*self%item(i)%a*self%l_chunk)
               
            end if   
         
         end if
      
      end do
      
      return
   
   end subroutine
   
   ! get absorption rate at path%length for ray object
   pure function a_dot_sph(self) result (a_dot_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                     ! ray object
      
      ! result declaration
      real(kind=rel_kind) :: a_dot_value                    ! absorption rate path%length      
      
      a_dot_value=sum(self%item(:self%n_item)%m*&
         &self%item(:self%n_item)%inv_rho*&
         &self%item(:self%n_item)%inv_h**(n_dim)*&
         &self%item(:self%n_item)%a_dot*&
         &self%sph_kernel%w(self%item(:self%n_item)%s_1))
      
      return
   
   end function
   
   ! get velocity at path%length for ray object
   pure function v_sph(self) result (v_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                     ! ray object
      
      ! result declaration
      real(kind=rel_kind) :: v_value(n_dim)             ! velocity path%length
      
      ! variable declaration
      integer(kind=int_kind) :: i                       ! counter
      
      
      v_value=0._rel_kind
      do i=1,self%n_item
      
         v_value=v_value+self%item(i)%m*&
            &self%item(i)%inv_rho*&
            &self%item(i)%inv_h**(n_dim)*&
            &self%item(i)%v*&
            &self%sph_kernel%w(self%item(i)%s_1)
         
      end do
      
      return
   
   end function
   
   ! get rest albedo at path%length for ray object
   pure function a_sph(self) result (a_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                     ! ray object
      
      ! result declaration
      real(kind=rel_kind) :: a_value                    ! albedo path%length      

      a_value=sum(self%item(:self%n_item)%m*&
         &self%item(:self%n_item)%inv_rho*&
         &self%item(:self%n_item)%inv_h**(n_dim)*&
         &self%item(:self%n_item)%a*&
         &self%sph_kernel%w(self%item(:self%n_item)%s_1))
      
      return
   
   end function
   
   ! get the optical depth along ray
   pure function tau_sph(self) result (tau_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                     ! ray object
      
      ! result declaration
      real(kind=rel_kind) :: tau_value                  ! optical depth along ray
      
      tau_value=sum(self%item(:self%n_item)%sigma*&
         &self%item(:self%n_item)%kappa_ext)
      
      return
      
   end function
   
   ! get inverse mean free path at path%length for ray object
   pure function inv_mfp(self) result (inv_mfp_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                     ! ray object
      
      ! result declaration
      real(kind=rel_kind) :: inv_mfp_value              ! inv_mfp at path%length
      
      inv_mfp_value=sum(self%item(:self%n_item)%m*&
         &self%item(:self%n_item)%inv_h**(n_dim)*&
         &self%item(:self%n_item)%kappa_ext*&
         &self%sph_kernel%w(self%item(:self%n_item)%s_1))
      
      return
   
   end function
   
   pure function grad_inv_mfp(self) result (grad_inv_mfp_value)
   
      ! argument declarations
      class(ray),intent(in) :: self                     ! ray object
      
      ! result declaration
      real(kind=rel_kind) :: grad_inv_mfp_value(n_dim)  !grad_inv_mfp at path%length
      
      ! variable declarations
      integer(kind=int_kind) :: i                       ! counter
    
      grad_inv_mfp_value=0._rel_kind
      do i=1,self%n_item
      
         grad_inv_mfp_value=grad_inv_mfp_value+self%item(i)%m*&
            &self%item(i)%inv_h**(n_dim+2)*&
            &self%item(i)%kappa_ext*&
            &self%sph_kernel%dw_dr(self%item(i)%s_1)*&
            &(self%path%terminus-self%item(i)%r)/self%item(i)%s_1
      
      end do
      
      return
   
   end function
   
   ! check if we need to perform modified random walk
   pure function mrw_check(self,gamma) result (use_mrw)
   
      ! argument declarations
      class(ray),intent(in) :: self                     ! ray object
      real(kind=rel_kind),intent(in) :: gamma           ! min number of mean free paths for mrw
      
      ! result declaration
      logical(kind=log_kind) :: use_mrw                 ! use modified random walk
      
      ! variable declarations
      real(kind=rel_kind) :: inv_mfp                    ! planck average inv_mfp at origin
      
      ! calculate inverse mean free path
      inv_mfp=sum(self%item(:self%n_item)%m*self%item(:self%n_item)%inv_h**n_dim*self%sph_kernel%w(self%item(:self%n_item)%s_0)*&
         &self%item(:self%n_item)%f_sub*self%dust_prop%mrw_inv_planck_ext(self%item(:self%n_item)%a_dot))
         
      use_mrw=(self%path%length>gamma*inv_mfp)
         
      return
   
   end function 
   
   ! get smoothing length at origin
   pure function h_0(self) result (h)
   
      ! argument declarations
      class(ray),intent(in) :: self                      ! ray object
      
      ! result declaration
      real(kind=rel_kind) :: h                           ! smoothing length
      
      h=sum(self%item(:self%n_item)%m*self%item(:self%n_item)%inv_h**(n_dim-1)*self%item(:self%n_item)%inv_rho*&
         &self%sph_kernel%w(self%item(:self%n_item)%s_0))
      
      return
         
   end function
   
   ! sort ray items in order of position along ray
   pure recursive subroutine ray_particle_sort(array)
   
      ! argument declarations
      type(ray_particle),intent(inout) :: array(:) ! array to be sorted
      
      ! variable declarations
      integer(kind=int_kind) :: i               ! counter
      integer(kind=int_kind) :: pivot           ! pivot element
      
      ! sort array
      if (size(array)>1) then

         ! set pivot to middle element
         pivot=1+(size(array)-1)/2
         
         ! move pivot to end of array
         call ray_particle_swap(array(pivot),array(size(array)))
         
         ! loop over array
         pivot=1
         do i=1,size(array)-1
         
            if (array(i)%t_r<=array(size(array))%t_r) then
            
               call ray_particle_swap(array(pivot),array(i))
               pivot=pivot+1
            
            end if
         
         end do
         
         ! move pivot back
         call ray_particle_swap(array(pivot),array(size(array)))
         
         ! sort rest of array
         call ray_particle_sort(array(:pivot-1))
         call ray_particle_sort(array(pivot+1:))
      
      end if      
      
      return
   
   end subroutine
   
   ! swaps ray item a with ray item b
   pure subroutine ray_particle_swap(a,b)
   
      ! argument declarations
      type(ray_particle),intent(inout) :: a
      type(ray_particle),intent(inout) :: b
      
      ! variable declarations
      type(ray_particle) :: temp
      
      temp=a
      a=b
      b=temp
      
      return
   
   end subroutine
   

end module