module m_line

   use m_kind_parameters
   
   implicit none
   
   private
   public :: line
   
   type :: line 
   
      real(kind=rel_kind) :: origin(n_dim)              ! origin of line
      real(kind=rel_kind) :: terminus(n_dim)            ! terminus of line
      real(kind=rel_kind) :: direction(n_dim)           ! direction unit vector
      real(kind=rel_kind) :: inv_direction(n_dim)       ! 1/direction
      real(kind=rel_kind) :: length                     ! length of line
      
      contains
      
      procedure,non_overridable :: initialise
   
   end type
   
   contains
   
   pure subroutine initialise(self,origin,direction,length)
   
      ! argument declarations
      class(line),intent(inout) :: self                 ! line object
      real(kind=rel_kind),intent(in) :: origin(n_dim)   ! origin of line
      real(kind=rel_kind),intent(in) :: direction(n_dim)! direction of line
      real(kind=rel_kind),intent(in) :: length          ! length of line
      
      ! set variables
      self%origin=origin
      self%direction=direction/norm2(direction)
      self%length=length
      self%inv_direction=1._rel_kind/self%direction
      self%terminus=self%origin+self%length*self%direction
      
      
      return
   
   end subroutine

end module