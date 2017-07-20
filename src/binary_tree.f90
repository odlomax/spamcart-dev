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

! binary tree module
module m_binary_tree

   use m_kind_parameters
   use m_sph_parameters
   use m_binary_tree_node
   use m_particle
   use m_kernel
   use omp_lib
   
   implicit none
   
   private
   public :: binary_tree
   
   ! define binary_tree class
   type :: binary_tree
   
      type(binary_tree_node),pointer :: root_node               ! root node of tree
      type(particle),contiguous,pointer :: particle_array(:)    ! pointer to particle array
      class(kernel),pointer :: sph_kernel                       ! pointer to kernel object
      real(kind=rel_kind) :: max_length                         ! diameter of circumsphere containing all particles
      real(kind=rel_kind) :: min_length                         ! shortest smoothing length
      real(kind=rel_kind) :: com(n_dim)                         ! centre of mass
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: destroy
      procedure,non_overridable :: calc_h
   
   end type
   
   contains
   
   ! initialise binary tree
   subroutine initialise(self,particle_array,sph_kernel,eta,h_present)
   
      ! argument declarations
      class(binary_tree),intent(inout) :: self                  ! binary tree object
      type(particle),intent(inout),target :: particle_array(:)  ! array of particles
      class(kernel),intent(inout),target :: sph_kernel          ! sph kernel
      real(kind=rel_kind),intent(in) :: eta                     ! smoothing length scale factor
      logical(kind=log_kind),intent(in) :: h_present            ! has smoothing length been calculated?
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      ! associate particle array
      self%particle_array=>particle_array
      
      ! associate kernel
      self%sph_kernel=>sph_kernel
      
      ! build tree from root node
      write(*,"(A)") "build tree"
      allocate(self%root_node)
      call self%root_node%initialise(particle_array)
      
      ! stock tree from root node (round 1)
      write(*,"(A)") "stock tree (round 1)"
      call self%root_node%stock(sph_kernel)      
      
      if (.not.h_present) then
      
         ! calculate smoothing lengths
         write(*,"(A)") "calculate h"
         call self%calc_h(eta)
      
         ! stock tree from root node (round 2)
         write(*,"(A)") "stock tree (round 2)"
         call self%root_node%stock(self%sph_kernel)
         
      end if
      
      
      ! calculate com, max_length and min_length
      self%com=0._rel_kind
      do i=1,size(self%particle_array)
         self%com=self%particle_array(i)%r*self%particle_array(i)%m
      end do
      self%com=self%com/sum(self%particle_array%m)
      self%max_length=0._rel_kind
      do i=1,size(self%particle_array)
         self%max_length=max((norm2(self%particle_array(i)%r-self%com)+self%particle_array(i)%h),self%max_length)
      end do
      self%max_length=self%max_length*2._rel_kind
      self%min_length=minval(self%particle_array%h)
      
      return
   
   end subroutine
   
   ! destroy binary tree
   pure subroutine destroy(self)
   
      ! argument declarations
      class(binary_tree),intent(inout) :: self                  ! binary tree object
      
      ! destroy the root node
      call self%root_node%destroy()
      deallocate(self%root_node)
      
      return
   
   end subroutine
   
   subroutine calc_h(self,eta)
   
      ! argument declarations
      class(binary_tree),intent(inout) :: self                  ! binary tree object
      real(kind=rel_kind),intent(in) :: eta                     ! h scale factor
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      integer(kind=int_kind) :: j                               ! counter
      integer(kind=int_kind) :: n_threads                       ! number of threads
      real(kind=rel_kind) :: delta_r(n_dim)                     ! size of box
      real(kind=rel_kind) :: delta_r_max                        ! largest dimension
      real(kind=rel_kind) :: h_init                             ! starting smoothing length
      real(kind=rel_kind) :: h_low                              ! lower bound for bisector
      real(kind=rel_kind) :: h_high                             ! uppber bound for bisector
      real(kind=rel_kind) :: h                                  ! smoothing length
      real(kind=rel_kind) :: h_new                              ! new smoothing length
      real(kind=rel_kind) :: rho                                ! density
      logical(kind=log_kind) :: found_h                         ! did we find h?
      
      ! make smoothing length guess from volume of system
      do i=1,n_dim
         delta_r(i)=maxval(self%particle_array%r(i))-minval(self%particle_array%r(i))
      end do
      ! divide out max dimension to avoid float overflow
      delta_r_max=maxval(delta_r)
      delta_r=delta_r/delta_r_max
      
      ! set initial guess of smoothing length
      h_init=eta*delta_r_max*(product(delta_r)/&
         &real(size(self%particle_array),rel_kind))**(1._rel_kind/real(n_dim,rel_kind))
      
      ! loop over all particles
      n_threads=omp_get_num_procs()
      !$omp parallel do schedule(dynamic) num_threads(n_threads) default(shared)&
      !$omp& private(i,j,h_low,h_high,h,h_new,rho,found_h)
      
         do i=1,size(self%particle_array)
         
            found_h=.false.
            if (self%particle_array(i)%inv_h>0._rel_kind) then
               h_new=1._rel_kind/self%particle_array(i)%inv_h
            else
               h_new=h_init
            end if
            h_low=huge(0._rel_kind)
            h_high=0._rel_kind
         
            ! find h using fixed point iteration
            do j=1,h_iterations
            
               
               h=h_new
               rho=self%root_node%sph_gather_density(self%sph_kernel,self%particle_array(i)%r,h)
               
               h_new=eta*(self%particle_array(i)%m/rho)**(1._rel_kind/real(n_dim,rel_kind))
               
               if (abs(h_new-h)/h<delta_h) then
               
                  found_h=.true.
                  exit
               
               end if
            
            end do
            
            ! use bisector
            if (.not.found_h) then
               
               ! get lower bound
               h_low=h
               do
               
                  rho=self%root_node%sph_gather_density(self%sph_kernel,self%particle_array(i)%r,h_low)
                  if (rho<eta**n_dim*self%particle_array(i)%m/h_low**n_dim) exit
                  h_low=0.5_rel_kind*h_low
               
               end do
               
               ! get upper bound
               h_high=h
               do
               
                  rho=self%root_node%sph_gather_density(self%sph_kernel,self%particle_array(i)%r,h_high)
                  if (rho>eta**n_dim*self%particle_array(i)%m/h_high**n_dim) exit
                  h_high=2._rel_kind*h_high
               
               end do
               
               do
                  
                  ! perform bisections
                  h=0.5_rel_kind*(h_high+h_low)
                  rho=self%root_node%sph_gather_density(self%sph_kernel,self%particle_array(i)%r,h)
                  
                  h_new=eta*(self%particle_array(i)%m/rho)**(1._rel_kind/real(n_dim,rel_kind))
                  if (abs(h_new-h)/h<delta_h) exit 
                  
                  if (rho>eta**n_dim*self%particle_array(i)%m/h**n_dim) then
                  
                     h_high=h
                     
                  else
                  
                     h_low=h
                  
                  end if               
                  
               end do
               
            end if
            
            self%particle_array(i)%h=h
            self%particle_array(i)%inv_h=1._rel_kind/h
            self%particle_array(i)%rho=rho
            self%particle_array(i)%inv_rho=1._rel_kind/rho
         
         end do
      
      !$omp end parallel do
      
      return
   
   end subroutine

end module