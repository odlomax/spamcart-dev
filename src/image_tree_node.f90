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

! node for image quadtree
module m_image_tree_node

   use m_kind_parameters
   use m_particle
   
   implicit none
      
   private
   public :: image_tree_node
   public :: proj_particle
   
   ! module parameters
   integer(kind=int_kind),parameter :: n_dim_im=2                  ! dimensions of image
   
   ! projected particle type
   type :: proj_particle
   
      real(kind=rel_kind) :: r(n_dim_im)                        ! projected particle position
      real(kind=rel_kind) :: h                                  ! particle smoothing length
      
      contains
      
      procedure,non_overridable :: initialise=>initialise_proj_particle
   
   end type
   
   ! define quadtree node type
   type :: image_tree_node
   
      real(kind=rel_kind) :: aabb(n_dim_im,2)                  ! axis aligned bounding box (grid volume)
      integer(kind=int_kind) :: n_leaf                         ! number of leaves on tree
      integer(kind=int_kind) :: level                          ! level of tree
      integer(kind=int_kind) :: max_level                      ! maximum level of tree
      real(kind=rel_kind) :: sigma                             ! column density
      real(kind=rel_kind),allocatable :: i_lambda(:)           ! intensity array
      real(kind=rel_kind),allocatable :: j_lambda(:)           ! intensity array
      type(proj_particle),pointer,contiguous :: particle_array(:) ! particle array
      type(image_tree_node),pointer,contiguous :: children(:)  ! child node array
      
      contains
      
      procedure,non_overridable :: initialise=>initialise_image_tree_node
      procedure,non_overridable :: destroy
      procedure,non_overridable :: get_leaf_aabb
      procedure,non_overridable :: set_leaf_i_lambda
      procedure,non_overridable :: get_node_pointer
      procedure,non_overridable :: write_leaf
      procedure,non_overridable :: rebuild_tree
   
   end type
   
   contains
   
   pure subroutine initialise_proj_particle(self,sph_particle,unit_vector)
   
      ! argument declarations
      class(proj_particle),intent(inout) :: self               ! projected particle class
      type(particle),intent(in) :: sph_particle                ! sph particle
      real(kind=rel_kind),intent(in) :: unit_vector(n_dim,n_dim_im)  ! unit vectors of image plane
      
      ! variable declarations
      integer(kind=int_kind) :: i                              ! counter
      
      ! set position
      do i=1,n_dim_im
         self%r(i)=dot_product(sph_particle%r,unit_vector(:,i))
      end do
      
      ! set smoothing length
      self%h=sph_particle%h
   
      return
   
   end subroutine
   
   pure recursive subroutine initialise_image_tree_node(self,particle_array,aabb,level)
   
      ! argument declarations
      class(image_tree_node),intent(inout) :: self             ! quadtree node
      type(proj_particle),intent(inout),target :: particle_array(:)! array of particles
      real(kind=rel_kind),intent(in) :: aabb(n_dim_im,2)       ! axis aligned bounding box
      integer(kind=int_kind),intent(in),optional :: level      ! level on tree
      
      ! set level
      if (present(level)) then
         self%level=level
      else
         self%level=0
      end if
      
      
      ! associate particle array
      self%particle_array=>particle_array
      
      ! set aabb
      self%aabb=aabb
      
      ! build more sells if minimum smoothing length is smaller than cell   
      if (size(self%particle_array)>0) then
      
         if (minval(self%particle_array%h)<maxval(self%aabb(:,2)-self%aabb(:,1))) then
   
            ! recursively build children
            allocate(self%children(2**n_dim_im))
            call build_children(self%particle_array,self%children,self%aabb,self%level)
            self%n_leaf=sum(self%children%n_leaf)
            self%max_level=maxval(self%children%max_level)
            return
                  
         end if
         
      end if
      
      self%children=>null()
      self%n_leaf=1
      self%max_level=self%level
   
      return
   
   end subroutine
   
   pure recursive subroutine destroy(self)
   
      ! argument declarations
      class(image_tree_node),intent(inout) :: self              ! quadtree node
      
      ! variable declarations
      integer(kind=int_kind) :: i                               ! counter
      
      if (allocated(self%i_lambda)) deallocate(self%i_lambda)
      if (allocated(self%j_lambda)) deallocate(self%j_lambda)
      
      if (associated(self%children)) then
      
         do i=1,size(self%children)
            call self%children(i)%destroy()
         end do
         deallocate(self%children)
      
      end if    
      
      return
   
   end subroutine
   
   ! get array of aabb of all leaves in tree
   pure recursive function get_leaf_aabb(self) result(aabb)
   
      ! argument declarations
      class(image_tree_node),intent(in) :: self                ! quadtree node
      
      ! result declaration
      real(kind=rel_kind) :: aabb(n_dim_im,2,self%n_leaf)      ! array of axis-aligned bounding boxes
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                            ! counter
      
      if (associated(self%children)) then
      
         j=1
         do i=1,size(self%children)
            aabb(:,:,j:j+self%children(i)%n_leaf-1)=self%children(i)%get_leaf_aabb()
            j=j+self%children(i)%n_leaf
         end do
            
      else
      
         aabb(:,:,1)=self%aabb
      
      end if
      
      return
   
   end function
   
   ! set assign i_lambda array to all leaves in tree
   pure recursive subroutine set_leaf_i_lambda(self,i_lambda,j_lambda,sigma)
   
      ! argument declarations
      class(image_tree_node),intent(inout) :: self          ! quadtree node
      real(kind=rel_kind),intent(in),contiguous :: i_lambda(:,:)  ! i_lambda array. second dimension must equal self%n_leaf
      real(kind=rel_kind),intent(in),contiguous :: j_lambda(:,:)  ! j_lambda_array
      real(kind=rel_kind),intent(in),contiguous :: sigma(:)       ! column density array
      
      ! variable declarations
      integer(kind=int_kind) :: i,j                         ! counter
      
      allocate(self%i_lambda(size(i_lambda,1)))
      allocate(self%j_lambda(size(j_lambda,1)))
      
      if (associated(self%children)) then
      
         ! recurse down tree and average i_lambda over children cells
         j=1
         self%i_lambda=0._rel_kind
         self%j_lambda=0._rel_kind
         self%sigma=0._rel_kind
         do i=1,size(self%children)
            call self%children(i)%set_leaf_i_lambda(i_lambda(:,j:j+self%children(i)%n_leaf-1),&
               &j_lambda(:,j:j+self%children(i)%n_leaf-1),sigma(j:j+self%children(i)%n_leaf-1))
            j=j+self%children(i)%n_leaf
            self%i_lambda=self%i_lambda+self%children(i)%i_lambda
            self%j_lambda=self%j_lambda+self%children(i)%j_lambda
            self%sigma=self%sigma+self%children(i)%sigma
         end do
         self%i_lambda=self%i_lambda/real(size(self%children),rel_kind)
         self%j_lambda=self%j_lambda/real(size(self%children),rel_kind)
         self%sigma=self%sigma/real(size(self%children),rel_kind)
      
      else
      
         self%i_lambda=i_lambda(:,1)
         self%j_lambda=j_lambda(:,1)
         self%sigma=sigma(1)
      
      end if
   
      return
   
   end subroutine
   
   ! returns pointer to the node containing r
   recursive function get_node_pointer(self,r,max_level_in) result(node_pointer)
   
      ! argument declarations
      class(image_tree_node),intent(in),target :: self   ! quadtree node
      real(kind=rel_kind),intent(in) :: r(n_dim_im)      ! position
      integer(kind=int_kind),intent(in),optional :: max_level_in  ! maximum recursion level
      
      ! result declaration
      class(image_tree_node),pointer :: node_pointer     ! i_lambda
      
      ! variable declarations
      integer(kind=int_kind) :: j                        ! counter
      integer(kind=int_kind) :: max_level                ! maximum recursion level
      
      ! set max level to max_level_in or really big number
      if (present(max_level_in)) then
         max_level=max_level_in
      else
         max_level=huge(1)
      end if
      
      if (associated(self%children).and.self%level<max_level) then
      
         ! find child cell which contains r
         do j=1,size(self%children)
         
            if (all(r>=self%children(j)%aabb(:,1).and.r<self%children(j)%aabb(:,2))) then
               node_pointer=>self%children(j)%get_node_pointer(r,max_level_in)
               return
            end if
         
         end do
         
         ! if r is not in child cells, pointer set to null
         node_pointer=>null()
         
      else
         
         ! found correct cell
         node_pointer=>self
         
      end if

      return
   
   end function
   
   ! rebuild branch nodes from modified leaves
   pure recursive subroutine rebuild_tree(self)
   
      ! argument declarations
      class(image_tree_node),intent(inout) :: self                ! quadtree node
      
      ! variable declarations
      integer(kind=int_kind) :: i                                 ! counter
      
      if (associated(self%children)) then
      
         self%i_lambda=0._rel_kind
         self%j_lambda=0._rel_kind
         do i=1,size(self%children)
         
            call self%children(i)%rebuild_tree()
            self%i_lambda=self%i_lambda+self%children(i)%i_lambda
            self%j_lambda=self%j_lambda+self%children(i)%j_lambda
         
         end do
         self%i_lambda=self%i_lambda/real(size(self%children),rel_kind)
         self%j_lambda=self%j_lambda/real(size(self%children),rel_kind)
      
      end if
      
      return
   
   end subroutine
   
   recursive subroutine write_leaf(self,file_id)
   
      ! argument declarations
      class(image_tree_node),intent(in) :: self                   ! quadtree node
      integer(kind=int_kind),intent(in) :: file_id                ! file identifier
      
      ! variable declarations
      integer(kind=int_kind) :: i                                 ! counter
      character(kind=chr_kind,len=string_length) :: format_string ! character format string
      
      if (associated(self%children)) then
      
         do i=1,size(self%children)
            call self%children(i)%write_leaf(file_id)
         end do
      
      else
      
         write(format_string,"(A,I0,A)") "(",size(self%aabb)+1+2*size(self%i_lambda),"(E25.17))"
         write(file_id,trim(format_string)) self%aabb,self%sigma,self%i_lambda,self%j_lambda
      
      end if 
      
      return
   
   end subroutine

   ! partition particle array about boundary along dimension k_dim
   pure subroutine partition_particles(particle_array,k_part,k_dim,boundary)
   
      ! argument declarations
      type(proj_particle),intent(inout),contiguous :: particle_array(:) ! array of particles
      integer(kind=int_kind),intent(out) :: k_part              ! partition partilce
      integer(kind=int_kind),intent(in) :: k_dim                ! splitting dimension
      real(kind=rel_kind),intent(in) :: boundary                ! split array about this boundary
      
      ! variable declarations
      type(proj_particle) :: temp_particle                      ! temporary particle
      integer(kind=int_kind) :: i                               ! counter
      
      k_part=0
      
      ! k_part = number of particles to the left of boundary
      do i=1,size(particle_array)
      
         if (particle_array(i)%r(k_dim)<boundary) then
         
            k_part=k_part+1
            temp_particle=particle_array(k_part)
            particle_array(k_part)=particle_array(i)
            particle_array(i)=temp_particle
         
         end if
      
      end do
      
      return
      
   end subroutine
   
   ! build children nodes
   pure recursive subroutine build_children(particle_array,children,aabb,tree_level,centre_in,level_in,lhs_in)
   
      ! argument declarations
      type(proj_particle),intent(inout),contiguous,target :: particle_array(:)   ! array of particles
      type(image_tree_node),intent(inout),contiguous :: children(:)         ! array of child nodes
      real(kind=rel_kind),intent(in) :: aabb(n_dim_im,2)                    ! parent axis aligned bounding box
      integer(kind=int_kind),intent(in) :: tree_level                       ! level on tree
      real(kind=rel_kind),intent(in),optional :: centre_in(n_dim_im)        ! centre of aabb
      integer(kind=int_kind),intent(in),optional :: level_in                ! level of recursion (1 to n_dim_im)
      logical(kind=log_kind),intent(inout),optional :: lhs_in(n_dim_im)     ! (left hand side) trace which cell we're in
      
      ! variable declarations
      integer(kind=int_kind) :: k_part                                      ! partition element
      integer(kind=int_kind) :: level                                       ! level of recursion
      real(kind=rel_kind) :: new_aabb(n_dim_im,2)                           ! new axis aligned bounding box
      real(kind=rel_kind) :: centre(n_dim_im)                               ! centre of aabb
      logical(kind=log_kind) :: lhs(n_dim_im)                              
      
      ! optional arguments aren't set on first call
      ! set aabb centre
      if (present(centre_in)) then
         centre=centre_in
      else
         centre=0.5_rel_kind*(aabb(:,1)+aabb(:,2))
      end if
      ! set level
      if (present(level_in)) then
         level=level_in
      else
         level=1
      end if
      ! set lhs
      if (present(lhs_in)) lhs=lhs_in
      
      
      ! if level is less than n_dim_im, split up the volume
      ! when level exceeds n_dim_im, assign particles to child nodes
      if (level<=n_dim_im) then
      
         ! partition particle array
         call partition_particles(particle_array,k_part,level,centre(level))
         
         ! split up particle array and children until children array has only one element   
         lhs(level)=.true.
         call build_children(particle_array(:k_part),&
            &children(:size(children)/2),aabb,tree_level,centre,level+1,lhs)
            
         lhs(level)=.false.
         call build_children(particle_array(k_part+1:),&
            &children(size(children)/2+1:),aabb,tree_level,centre,level+1,lhs)
         
      else
      
         ! children array now only has a single element
         new_aabb(:,1)=merge(aabb(:,1),centre,lhs)
         new_aabb(:,2)=merge(centre,aabb(:,2),lhs)
         call children(1)%initialise(particle_array,new_aabb,tree_level+1)
      
      end if
      
      return
   
   end subroutine

end module