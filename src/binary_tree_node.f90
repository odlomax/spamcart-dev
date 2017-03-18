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

! binary tree node module
module m_binary_tree_node

   use m_kind_parameters
   use m_sph_parameters
   use m_constants_parameters
   use m_particle
   use m_kernel
   use m_line
   use m_dust
   use m_source

   implicit none
   
   private
   public :: binary_tree_node
   
   ! class definition for a binary tree node
   type :: binary_tree_node
      
      real(kind=rel_kind) :: aabb(n_dim,2)              ! axis aligned bounding box
   
      type(particle),contiguous,pointer :: particle_array(:)   ! array of particles
      type(binary_tree_node),pointer :: left_node       ! node for left hand side or particle array
      type(binary_tree_node),pointer :: right_node      ! node for right hand side of particle array
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: destroy
      procedure,non_overridable :: stock
      procedure,non_overridable :: sph_gather_density
      procedure,non_overridable :: sph_gather_f_sub_variable
      procedure,non_overridable :: sph_gather_f_sub_fixed
      procedure,non_overridable :: sph_scatter
      procedure,non_overridable :: sph_inv_mfp
      procedure,non_overridable :: ray_aabb_intersection
      procedure,non_overridable,private :: particle_variance
      procedure,non_overridable,private :: swap_particle
      procedure,non_overridable,private :: partition
   
   end type
   
   contains
   
   ! build all nodes of binary tree beneath this one
   pure recursive subroutine initialise(self,particle_array)
   
      ! argument declarations
      class(binary_tree_node),intent(inout) :: self               ! binary tree node object
      type(particle),intent(inout),contiguous,target :: particle_array(:)  ! array of particles
      
      ! variable declarations
      integer(kind=int_kind) :: k_dim                           ! splitting dimension
      integer(kind=int_kind) :: k_part                          ! partition index
      
      ! associate particle array
      self%particle_array=>particle_array
      
      ! build more tree if this is not a leaf node
      if (size(self%particle_array)>n_leaf) then
         
         ! get splitting dimension
         k_dim=maxloc(self%particle_variance(),1)
         
         ! partition array (split about median value)
         k_part=size(self%particle_array)/2
         call self%partition(k_part,k_dim)
         
         ! build left branch of tree
         allocate(self%left_node)
         call self%left_node%initialise(particle_array(:k_part))
         
         ! build right branch of tree
         allocate(self%right_node)
         call self%right_node%initialise(particle_array(k_part+1:))
         
      else
      
         ! point branches to null
         self%left_node=>null()
         self%right_node=>null()
      
      end if
        
      return
   
   end subroutine
   
   ! destroy all nodes of binary tree beneath this one
   pure recursive subroutine destroy(self)
   
      ! argument declarations
      class(binary_tree_node),intent(inout) :: self             ! binary tree node object
      
      ! check if node has children
      if (associated(self%left_node)) then
      
         call self%left_node%destroy()
         call self%right_node%destroy()
         
         deallocate(self%left_node,self%right_node)
         
      end if
   
      return
   
   end subroutine
   
   ! stock all nodes of binary tree
   pure recursive subroutine stock(self,sph_kernel)
   
      ! argument declarations
      class(binary_tree_node),intent(inout) :: self             ! binary tree node object
      class(kernel),intent(in) :: sph_kernel                    ! sph kernel object
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      ! check if node has children
      if (associated(self%left_node)) then
      
         call self%left_node%stock(sph_kernel)
         call self%right_node%stock(sph_kernel)
         
         self%aabb(:,1)=min(self%left_node%aabb(:,1),self%right_node%aabb(:,1))
         self%aabb(:,2)=max(self%left_node%aabb(:,2),self%right_node%aabb(:,2))
         
      else
      
         do i=1,n_dim
            self%aabb(i,1)=minval(self%particle_array%r(i)-&
               &self%particle_array%h*sph_kernel%r_support)
            self%aabb(i,2)=maxval(self%particle_array%r(i)+&
               &self%particle_array%h*sph_kernel%r_support)
         end do
      
      end if
   
      return
   
   end subroutine
   
   ! peform sph gather to calculate density
   pure recursive function sph_gather_density(self,sph_kernel,r,h_gath) result (rho)
   
      ! argument declarations
      class(binary_tree_node),intent(in) :: self                ! binary tree node object
      class(kernel),intent(in) :: sph_kernel                    ! sph kernel object
      real(kind=rel_kind),intent(in) :: r(n_dim)                ! position
      real(kind=rel_kind),intent(in) :: h_gath                  ! gather smoothing length
      
      ! result declaration
      real(kind=rel_kind) :: rho                                ! density
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      rho=0._rel_kind
      
      ! check if virtual particle overlaps aabb of node
      if (all(r+h_gath*sph_kernel%r_support>self%aabb(:,1).and.&
         r-h_gath*sph_kernel%r_support<self%aabb(:,2))) then
      
         if (associated(self%left_node)) then
         
            ! recurse down tree
            rho=rho+self%left_node%sph_gather_density(sph_kernel,r,h_gath)
            rho=rho+self%right_node%sph_gather_density(sph_kernel,r,h_gath)
            
         else
         
            ! calculate quantities from leaf
            do i=1,size(self%particle_array)
            
               rho=rho+self%particle_array(i)%m*&
                  &sph_kernel%w(norm2(r-self%particle_array(i)%r)/h_gath)/h_gath**n_dim
                  
            end do
         
         end if
         
      end if
      
      return
   
   end function
   
   ! peform sph gather to calculate smoothed sublimation fraction
   pure recursive function sph_gather_f_sub_variable(self,sph_kernel,sph_particle,a_dot_sub) result (f_sub)
   
      ! argument declarations
      class(binary_tree_node),intent(in) :: self                ! binary tree node object
      class(kernel),intent(in) :: sph_kernel                    ! sph kernel object
      type(particle),intent(in) :: sph_particle                 ! sph particle
      real(kind=rel_kind),intent(in) :: a_dot_sub               ! sublimation energy absorption rate 
      
      ! result declaration
      real(kind=rel_kind) :: f_sub                              ! sublimation fraction
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      f_sub=0._rel_kind
      
      ! check if virtual particle overlaps aabb of node
      if (all(sph_particle%r+sph_particle%h*sph_kernel%r_support>self%aabb(:,1).and.&
         sph_particle%r-sph_particle%h*sph_kernel%r_support<self%aabb(:,2))) then
      
         if (associated(self%left_node)) then
         
            ! recurse down tree
            f_sub=f_sub+self%left_node%sph_gather_f_sub_variable(sph_kernel,sph_particle,a_dot_sub)
            f_sub=f_sub+self%right_node%sph_gather_f_sub_variable(sph_kernel,sph_particle,a_dot_sub)
            
         else
         
            ! calculate quantities from leaf
            do i=1,size(self%particle_array)
               f_sub=f_sub+merge(1._rel_kind,epsilon(0._rel_kind),self%particle_array(i)%a_dot<a_dot_sub)*&
                  &self%particle_array(i)%m*sph_kernel%w(norm2(sph_particle%r-self%particle_array(i)%r)/sph_particle%h)/&
                  &(sph_particle%h**n_dim*sph_particle%rho)
            end do
         end if
         
      end if
      
      return
   
   end function
   
   ! peform sph gather to calculate smoothed sublimation fraction
   pure recursive function sph_gather_f_sub_fixed(self,sph_kernel,sph_particle,point_source_array) result (f_sub)
   
      ! argument declarations
      class(binary_tree_node),intent(in) :: self                ! binary tree node object
      class(kernel),intent(in) :: sph_kernel                    ! sph kernel object
      type(particle),intent(in) :: sph_particle                 ! sph particle
      class(source),intent(in) :: point_source_array(:)         ! array of point sources
      
      ! result declaration
      real(kind=rel_kind) :: f_sub                              ! sublimation fraction
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                             ! counter
      logical(kind=log_kind) :: in_radius                       ! is particle within a sublimation radius?
      
      f_sub=0._rel_kind
      
      ! check if virtual particle overlaps aabb of node
      if (all(sph_particle%r+sph_particle%h*sph_kernel%r_support>self%aabb(:,1).and.&
         sph_particle%r-sph_particle%h*sph_kernel%r_support<self%aabb(:,2))) then
      
         if (associated(self%left_node)) then
         
            ! recurse down tree
            f_sub=f_sub+self%left_node%sph_gather_f_sub_fixed(sph_kernel,sph_particle,point_source_array)
            f_sub=f_sub+self%right_node%sph_gather_f_sub_fixed(sph_kernel,sph_particle,point_source_array)
            
         else
         
            ! calculate quantities from leaf
            do i=1,size(self%particle_array)
            
               in_radius=.false.
               do j=1,size(point_source_array)
                  if (sum((self%particle_array(i)%r-point_source_array(j)%position)**2)<=point_source_array(j)%radius**2) then
                     in_radius=.true.
                     exit
                  end if
               end do   
            
            
               f_sub=f_sub+merge(epsilon(0._rel_kind),1._rel_kind,in_radius)*&
                  &self%particle_array(i)%m*sph_kernel%w(norm2(sph_particle%r-self%particle_array(i)%r)/sph_particle%h)/&
                  &(sph_particle%h**n_dim*sph_particle%rho)
            end do
         end if
         
      end if
      
      return
   
   end function
   
   ! peform sph scatter to calculate a quantity
   pure recursive subroutine sph_scatter(self,sph_kernel,r,v,a_dot,rho)
   
      ! argument declarations
      class(binary_tree_node),intent(in) :: self                ! binary tree node object
      class(kernel),intent(in) :: sph_kernel                    ! sph kernel object
      real(kind=rel_kind),intent(in) :: r(n_dim)                ! position
      real(kind=rel_kind),intent(inout),optional :: v(n_dim)    ! velocity
      real(kind=rel_kind),intent(inout),optional :: a_dot       ! energy absorption rate
      real(kind=rel_kind),intent(inout),optional :: rho         ! density
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      real(kind=rel_kind) :: w                                  ! kernel density
      real(kind=rel_kind) :: w_dv                               ! kernel mass
      
      ! check if position is within aabb of node
      if (any(r<self%aabb(:,1)).or.any(r>self%aabb(:,2))) return
      
      if (associated(self%left_node)) then
      
         ! recurse down tree
         call self%left_node%sph_scatter(sph_kernel,r,v,a_dot,rho)
         call self%right_node%sph_scatter(sph_kernel,r,v,a_dot,rho)
         
      else
      
         ! calculate quantities from leaf
         do i=1,size(self%particle_array)
         
            w=sph_kernel%w&
               (sqrt(sum((r-self%particle_array(i)%r)**2))*self%particle_array(i)%inv_h)*&
               &self%particle_array(i)%inv_h**(n_dim)
            w_dv=w*self%particle_array(i)%m*self%particle_array(i)%inv_rho
            
            if (present(v)) v=v+self%particle_array(i)%v*w_dv
            if (present(a_dot)) a_dot=a_dot+self%particle_array(i)%a_dot*w_dv
            if (present(rho)) rho=rho+self%particle_array(i)%m*w
               
         end do
      
      end if
      
      return
   
   end subroutine
   
   ! calculate the inverse mean free path, and its gradient
   pure recursive subroutine sph_inv_mfp(self,sph_kernel,dust_prop,&
      &r,lambda_em,n_em,v_em,inv_mfp,grad_inv_mfp,inv_h,ave_inv_mfp,grad_ave_inv_mfp)
   
      ! argument declarations
      class(binary_tree_node),intent(in) :: self                ! binary tree node object
      class(kernel),intent(in) :: sph_kernel                    ! sph kernel object
      class(dust),intent(in) :: dust_prop                       ! dust properties object
      real(kind=rel_kind),intent(in) :: r(n_dim)                ! position
      real(kind=rel_kind),intent(in) :: lambda_em               ! wavelength of emission
      real(kind=rel_kind),intent(in) :: n_em(n_dim)             ! direction of emission
      real(kind=rel_kind),intent(in) :: v_em(n_dim)             ! velocity of emission
      real(kind=rel_kind),intent(inout) :: inv_mfp              ! inverse mean free path
      real(kind=rel_kind),intent(inout) :: grad_inv_mfp(n_dim)  ! gradient of inverse mean free path
      real(kind=rel_kind),intent(inout) :: inv_h                ! inverse smoothing length to the power n_dim
      real(kind=rel_kind),intent(inout) :: ave_inv_mfp          ! planck averaged mean free path
      real(kind=rel_kind),intent(inout) :: grad_ave_inv_mfp(n_dim)   ! gradient of planck average mean free path
       
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      real(kind=rel_kind) :: s                                  ! normalised radius
      real(kind=rel_kind) :: w                                  ! kernel density
      real(kind=rel_kind) :: grad_w(n_dim)                      ! kernel density gradient
      real(kind=rel_kind) :: lambda_ob                          ! wavelength in rest frame of particle
      real(kind=rel_kind) :: dust_mass_ext                      ! mass extinction coefficient
      real(kind=rel_kind) :: inv_planck_ext                     ! inverse mean planck extinction
      
      ! check if position is within aabb of node
      if (any(r<self%aabb(:,1)).or.any(r>self%aabb(:,2))) return
      
      if (associated(self%left_node)) then
      
         ! recurse down tree
         call self%left_node%sph_inv_mfp(sph_kernel,dust_prop,&
            &r,lambda_em,n_em,v_em,inv_mfp,grad_inv_mfp,inv_h,ave_inv_mfp,grad_ave_inv_mfp)
         call self%right_node%sph_inv_mfp(sph_kernel,dust_prop,&
            &r,lambda_em,n_em,v_em,inv_mfp,grad_inv_mfp,inv_h,ave_inv_mfp,grad_ave_inv_mfp)
         
      else
      
         ! calculate quantities from leaf
         do i=1,size(self%particle_array)
         
            s=norm2(r-self%particle_array(i)%r)*self%particle_array(i)%inv_h
            w=sph_kernel%w(s)*self%particle_array(i)%inv_h**(n_dim)
            grad_w=sph_kernel%dw_dr(s)*self%particle_array(i)%inv_h**(n_dim+2)*&
               &(r-self%particle_array(i)%r)/s
            
            lambda_ob=lambda_em*&
               &(dot_product(self%particle_array(i)%v-v_em,n_em)/c_light+1._rel_kind)
            dust_mass_ext=dust_prop%dust_mass_ext(lambda_ob)*self%particle_array(i)%f_sub
            inv_planck_ext=dust_prop%mrw_inv_planck_ext(self%particle_array(i)%a_dot)*self%particle_array(i)%f_sub
            
            inv_mfp=inv_mfp+self%particle_array(i)%m*dust_mass_ext*w
            grad_inv_mfp=grad_inv_mfp+self%particle_array(i)%m*dust_mass_ext*grad_w
            
            inv_h=inv_h+w
            
            ave_inv_mfp=ave_inv_mfp+self%particle_array(i)%m*inv_planck_ext*w
            grad_ave_inv_mfp=grad_ave_inv_mfp+self%particle_array(i)%m*inv_planck_ext*grad_w
            
         end do
      
      end if
      
      return
   
   end subroutine 
   
   ! particle position variance along all cardinal dimensions
   pure function particle_variance(self) result (sigma_sqrd)
   
      ! argument declarations
      class(binary_tree_node),intent(in) :: self                ! binary tree node object
      
      ! result declaration
      real(kind=rel_kind) :: sigma_sqrd(n_dim)                  ! particle position variance
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      real(kind=rel_kind) :: mu(n_dim)                          ! particle position mean
      
      ! calculate mean
      mu=0._rel_kind
      do i=1,size(self%particle_array)
         mu=mu+self%particle_array(i)%r
      end do
      mu=mu/real(size(self%particle_array),rel_kind)
      
      ! calculate variance
      sigma_sqrd=0._rel_kind
      do i=1,size(self%particle_array)
         sigma_sqrd=sigma_sqrd+(self%particle_array(i)%r-mu)**2
      end do
      sigma_sqrd=sigma_sqrd/real(size(self%particle_array),rel_kind)
      
      return
   
   end function
   
   ! swap two particles in array
   pure subroutine swap_particle(self,i,j)
   
      ! argument declarations
      class(binary_tree_node),intent(inout) :: self            ! binary tree node object
      integer(kind=int_kind),intent(in) :: i                   ! particle index
      integer(kind=int_kind),intent(in) :: j                   ! particle index
      
      ! variable declarations
      type(particle) :: temp_particle                          ! temporary particle
      
      temp_particle=self%particle_array(i)
      self%particle_array(i)=self%particle_array(j)
      self%particle_array(j)=temp_particle
      
      return
   
   end subroutine
   
   ! partition particle array about k_part along dimension k_dim
   pure subroutine partition(self,k_part,k_dim)
   
      ! argument declarations
      class(binary_tree_node),intent(inout) :: self             ! binary tree node object
      integer(kind=int_kind),intent(inout) :: k_part            ! partition element
      integer(kind=int_kind),intent(in) ::    k_dim             ! splitting dimension
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      integer(kind=int_kind) :: i_min,i_max                     ! bounds on i
      integer(kind=int_kind) :: k_part_guess                    ! guess value of k_part
      
      i_min=1
      i_max=size(self%particle_array)
      
      ! perform quickselect partition
      do
      
         ! make a guess of k_part
         k_part_guess=i_min+(i_max-i_min)/2
         
         ! move k_part_guess to end of array
         call self%swap_particle(k_part_guess,i_max)
         
         ! loop over array
         k_part_guess=i_min
         do i=i_min,i_max-1
            if (self%particle_array(i)%r(k_dim)<=self%particle_array(i_max)%r(k_dim)) then
               call self%swap_particle(k_part_guess,i)
               k_part_guess=k_part_guess+1
            end if
         end do
         
         ! move k_part_guess back
         call self%swap_particle(k_part_guess,i_max)
      
         if (k_part_guess==k_part) then
            exit                        ! guess was correct
         elseif (k_part_guess<k_part) then
            i_min=k_part_guess+1        ! guess was too low
         else
            i_max=k_part_guess-1        ! guess was too high
         end if
      
      end do
  
      return
   
   end subroutine
   
   ! check if ray intersects aabb
   elemental function ray_aabb_intersection(self,ray) result (intersection)
   
      ! argument declations
      class(binary_tree_node),intent(in) :: self                ! binary tree node object
      type(line),intent(in) :: ray                              ! ray
      
      ! result declaration
      logical(kind=log_kind) :: intersection                    ! true if ray intersects aabb
      
      ! variable declarations
      real(kind=rel_kind) :: t_0(n_dim)                         ! slab intersection points
      real(kind=rel_kind) :: t_1(n_dim)                         ! slab intersection points
      real(kind=rel_kind) :: t_min(n_dim)                       ! lower intersection points
      real(kind=rel_kind) :: t_max(n_dim)                       ! upper intersection points
      real(kind=rel_kind) :: t_min_scalar                       ! scalar t_min
      real(kind=rel_kind) :: t_max_scalar                       ! scalar t_max
      
      ! calc slab intersection points
      t_0=(self%aabb(:,1)-ray%origin)*ray%inv_direction
      t_1=(self%aabb(:,2)-ray%origin)*ray%inv_direction
      
      ! calculation minima and maxima of interection points
      t_min=min(t_0,t_1)
      t_max=max(t_0,t_1)
      
      ! line intersects aabb if maxval(t_min) < minval(t_max)
      t_min_scalar=maxval(t_min)
      t_max_scalar=minval(t_max)
      
      ! also perform checks for line end points
      intersection=(t_min_scalar<t_max_scalar).and.&
         &(t_max_scalar>0._rel_kind).and.&
         &(t_min_scalar<ray%length)
   
      return
   
   end function
   
end module