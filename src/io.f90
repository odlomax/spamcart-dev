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

! file input/output module
module m_io
   
   use m_kind_parameters
   use m_maths

   implicit none
   
   private
   public :: read_in_sph_particles_3d
   public :: read_in_point_sources_3d
   public :: write_out_sph_particles_3d
   public :: write_out_datacube_3d

   integer(kind=int_kind),parameter :: h_str=1000                             ! max length of header line
   
   real(kind=rel_kind),parameter :: m_in_cm=1.e+2                             ! m in cm
   real(kind=rel_kind),parameter :: pc_in_cm=3.0856776e+18                    ! parsec in cm
   real(kind=rel_kind),parameter :: au_in_cm=1.4959787e+13                    ! au in cm
   real(kind=rel_kind),parameter :: kg_in_g=1.e+3                             ! kg in g
   real(kind=rel_kind),parameter :: msun_in_g=1.9891e+33                      ! solar mass in g
   real(kind=rel_kind),parameter :: w_in_ergs=1.e+7                           ! W in erg/s
   real(kind=rel_kind),parameter :: z_sun=1.34e-2_rel_kind                    ! solar metalicity
   real(kind=rel_kind),parameter :: s_in_year=3.1556926e+7                    ! seconds in a year      
   
   contains
   
   ! read in sph particles from file
   subroutine read_in_sph_particles_3d(position,mass,temperature,file_name)
   
      ! argument declarations
      real(kind=rel_kind),intent(out),allocatable :: position(:,:)            ! particle positions
      real(kind=rel_kind),intent(out),allocatable :: mass(:)                  ! particles masses
      real(kind=rel_kind),intent(out),allocatable :: temperature(:)           ! particle temperatures
      character(kind=chr_kind,len=string_length),intent(in) :: file_name                        ! name of input file
      
      ! variable declarations
      integer(kind=int_kind) :: i                                             ! counter
      integer(kind=int_kind) :: h_space                                       ! index of space in header
      integer(kind=int_kind) :: u_space                                       ! index of space in unit line
      integer(kind=int_kind) :: n_col                                         ! number of columns
      integer(kind=int_kind) :: n_row                                         ! number of rows
      integer(kind=int_kind) :: c_x1                                          ! position columns
      integer(kind=int_kind) :: c_x2
      integer(kind=int_kind) :: c_x3
      integer(kind=int_kind) :: c_m                                           ! mass column
      integer(kind=int_kind) :: c_t                                           ! temperature column
      integer(kind=int_kind) :: read_status                                   ! read status of file
      real(kind=rel_kind) :: unit_factor(5)                                   ! unit conversion factor
      real(kind=rel_kind),allocatable :: data(:,:)                            ! data array from file
      character(len=h_str) :: header                                    ! header line from file
      character(len=h_str) :: unit_line                                 ! units line from file
      character(kind=chr_kind,len=string_length) :: format_string                               ! character format string
      
      ! input routine converts columns into cgs units
      ! unit_factor(1) = position(1) conversion factor
      ! unit_factor(2) = position(2) conversion factor
      ! unit_factor(3) = position(3) conversion factor
      ! unit_factor(4) = mass conversion factor
      ! unit_factor(5) = temperature conversion factor
      
      unit_factor=1._rel_kind
      
      ! open parameters file
      open(1,file=trim(file_name),status="old",form="formatted")
      
      ! set format string
      write(format_string,"(A2,I0,A1)") "(A",h_str,")"
      
      ! read header line
      read(1,trim(format_string)) header
      
      ! read unit line
      read(1,trim(format_string)) unit_line
      
      ! use header line to figure out column numbers
      header=adjustl(header)
      unit_line=adjustl(unit_line)
      n_col=0
      c_x1=0
      c_x2=0
      c_x3=0
      c_m=0
      c_t=0
      
      do
      
         ! find whitespace index
         h_space=index(header," ")
         u_space=index(unit_line," ")
         
         ! exit if first character is whitespace
         if (h_space==1) exit
         
         ! count number of columns
         n_col=n_col+1
         
         ! set column numbers
         select case (header(:h_space-1))
         
            case ("x_1")
               c_x1=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("cm")
                     unit_factor(1)=1._rel_kind
                     
                  case ("m")
                     unit_factor(1)=m_in_cm
               
                  case ("au")
                     unit_factor(1)=au_in_cm
                     
                  case ("pc")
                     unit_factor(1)=pc_in_cm
                     
                  case default
                    write(6,"(A)") "Error: unknown length unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
                     
               
            case ("x_2")
               c_x2=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("cm")
                     unit_factor(2)=1._rel_kind
               
                  case ("m")
                     unit_factor(2)=m_in_cm
               
                  case ("au")
                     unit_factor(2)=au_in_cm
                     
                  case ("pc")
                     unit_factor(2)=pc_in_cm
                     
                  case default
                    write(6,"(A)") "Error: unknown length unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
               
            case ("x_3")
               c_x3=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("cm")
                     unit_factor(3)=1._rel_kind
               
                  case ("m")
                     unit_factor(3)=m_in_cm
               
                  case ("au")
                     unit_factor(3)=au_in_cm
                     
                  case ("pc")
                     unit_factor(3)=pc_in_cm
                     
                  case default
                    write(6,"(A)") "Error: unknown length unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
               
            case ("m")
               c_m=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("g")
                     unit_factor(4)=1._rel_kind
               
                  case ("kg")
                     unit_factor(4)=kg_in_g
                     
                  case ("M_sun")
                     unit_factor(4)=msun_in_g
                     
                  case default
                    write(6,"(A)") "Error: unknown mass unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
               
            case ("T")
               c_t=n_col
               select case (unit_line(:u_space-1))
               
                  case ("K")
                     unit_factor(5)=1._rel_kind
                     
                  case default
                    write(6,"(A)") "Error: unknown temperature unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
               
         end select
         
         ! remove column label from header line
         header=adjustl(header(h_space+1:))
         unit_line=adjustl(unit_line(u_space+1:))
                  
      end do
      
      ! read through rest of file to get number of rows
      n_row=0
      do
      
         read(1,*,iostat=read_status)
         
         if (read_status<0) exit
         
         n_row=n_row+1
         
      end do
      
      
      ! allocate data arrays
      allocate(data(n_col,n_row))
      
      ! make sure arrays haven't already been allocated
      if (allocated(position)) deallocate(position)
      if (allocated(mass)) deallocate(mass)
      if (allocated(temperature)) deallocate(temperature)
      
      ! allocate arrays
      allocate(position(3,n_row))
      allocate(mass(n_row))
      allocate(temperature(n_row))
      
      ! rewind file
      rewind(1)
      ! skip header
      read(1,*)
      read(1,*)
      
      ! read in data block
      do i=1,n_row
      
         read(1,*) data(:,i)
         
      end do
      
      ! close file
      close(1)
      
      ! assign to arrays
      position(1,:)=data(c_x1,:)*unit_factor(1)
      position(2,:)=data(c_x2,:)*unit_factor(2)
      position(3,:)=data(c_x3,:)*unit_factor(3)
      mass=data(c_m,:)*unit_factor(4)
      temperature=data(c_t,:)*unit_factor(5)
      
      ! deallocate data block
      deallocate(data)
      
      return
   
   end subroutine
   
   ! read in point_sources from file
   subroutine read_in_point_sources_3d(position,temperature,luminosity,mass,age,metallicity,file_name)
   
      ! argument declarations
      real(kind=rel_kind),intent(out),allocatable :: position(:,:)            ! particle positions
      real(kind=rel_kind),intent(out),allocatable,optional :: temperature(:)  ! particle temperatures
      real(kind=rel_kind),intent(out),allocatable,optional :: luminosity(:)   ! particles luminosities
      real(kind=rel_kind),intent(out),allocatable,optional :: mass(:)         ! stellar mass
      real(kind=rel_kind),intent(out),allocatable,optional :: age(:)          ! stellar age
      real(kind=rel_kind),intent(out),allocatable,optional :: metallicity(:)  ! stellar metallicity
      character(kind=chr_kind,len=string_length),intent(in) :: file_name      ! name of input file
      
      ! variable declarations
      integer(kind=int_kind) :: i                                             ! counter
      integer(kind=int_kind) :: h_space                                       ! index of space in header
      integer(kind=int_kind) :: u_space                                       ! index of space in unit line
      integer(kind=int_kind) :: n_col                                         ! number of columns
      integer(kind=int_kind) :: n_row                                         ! number of rows
      integer(kind=int_kind) :: c_x1                                          ! position columns
      integer(kind=int_kind) :: c_x2
      integer(kind=int_kind) :: c_x3
      integer(kind=int_kind) :: c_l                                           ! luminosity column
      integer(kind=int_kind) :: c_t                                           ! temperature column
      integer(kind=int_kind) :: c_m                                           ! mass column
      integer(kind=int_kind) :: c_a                                           ! age column
      integer(kind=int_kind) :: c_z                                           ! metallicity column
      integer(kind=int_kind) :: read_status                                   ! read status of file
      real(kind=rel_kind) :: unit_factor(8)                                   ! unit conversion factor
      real(kind=rel_kind),allocatable :: data(:,:)                            ! data array from file
      character(len=h_str) :: header                                    ! header line from file
      character(len=h_str) :: unit_line                                 ! units line from file
      character(kind=chr_kind,len=string_length) :: format_string                               ! character format string
      
      ! input routine converts columns into cgs units
      ! unit_factor(1) = position(1) conversion factor
      ! unit_factor(2) = position(2) conversion factor
      ! unit_factor(3) = position(3) conversion factor
      ! unit_factor(4) = luminosity conversion factor
      ! unit_factor(5) = temperature conversion factor
      ! unit_factor(6) = mass conversion factor
      ! unit_factor(7) = age conversion factor
      ! unit_factor(8) = metallicity conversion factor
      
      unit_factor=1._rel_kind
      
      
      ! open parameters file
      open(1,file=trim(file_name),status="old",form="formatted")
      
      ! set format string
      write(format_string,"(A2,I0,A1)") "(A",h_str,")"
      
      ! read header line
      read(1,trim(format_string)) header
      
      ! read unit line
      read(1,trim(format_string)) unit_line
      
      ! use header line to figure out column numbers
      header=adjustl(header)
      unit_line=adjustl(unit_line)
      n_col=0
      c_x1=0
      c_x2=0
      c_x3=0
      c_l=0
      c_t=0
      c_m=0
      c_a=0
      c_z=0
      
      do
      
         ! find whitespace index
         h_space=index(header," ")
         u_space=index(unit_line," ")
         
         ! exit if first character is whitespace
         if (h_space==1) exit
         
         ! count number of columns
         n_col=n_col+1
         
         ! set column numbers
         
         select case (header(:h_space-1))
         
            case ("x_1")
               c_x1=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("cm")
                     unit_factor(1)=1._rel_kind
                     
                  case ("m")
                     unit_factor(1)=m_in_cm
               
                  case ("au")
                     unit_factor(1)=au_in_cm
                     
                  case ("pc")
                     unit_factor(1)=pc_in_cm
                     
                  case default
                    write(6,"(A)") "Error: unknown length unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
                     
               
            case ("x_2")
               c_x2=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("cm")
                     unit_factor(2)=1._rel_kind
               
                  case ("m")
                     unit_factor(2)=m_in_cm
               
                  case ("au")
                     unit_factor(2)=au_in_cm
                     
                  case ("pc")
                     unit_factor(2)=pc_in_cm
                     
                  case default
                    write(6,"(A)") "Error: unknown length unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
               
            case ("x_3")
               c_x3=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("cm")
                     unit_factor(3)=1._rel_kind
               
                  case ("m")
                     unit_factor(3)=m_in_cm
               
                  case ("au")
                     unit_factor(3)=au_in_cm
                     
                  case ("pc")
                     unit_factor(3)=pc_in_cm
                     
                  case default
                    write(6,"(A)") "Error: unknown length unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
               
            case ("L")
               c_l=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("erg_s^-1")
                     unit_factor(4)=1._rel_kind
               
                  case ("W")
                     unit_factor(4)=w_in_ergs
                     
                  case default
                    write(6,"(A)") "Error: unknown luminosity unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
               
            case ("T")
               c_t=n_col
               select case (unit_line(:u_space-1))
               
                  case ("K")
                     unit_factor(5)=1._rel_kind
                     
                  case default
                    write(6,"(A)") "Error: unknown temperature unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
                
            case ("M")
               c_m=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("g")
                     unit_factor(6)=1._rel_kind/msun_in_g
               
                  case ("kg")
                     unit_factor(6)=kg_in_g/msun_in_g
                     
                  case ("M_sun")
                     unit_factor(6)=1._rel_kind
                     
                  case default
                    write(6,"(A)") "Error: unknown mass unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
                
            case ("Age")
               c_a=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("s")
                     unit_factor(7)=1._rel_kind/s_in_year
               
                  case ("yr")
                     unit_factor(7)=1._rel_kind
                     
                  case default
                    write(6,"(A)") "Error: unknown time unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
                
            case ("Z")
               c_z=n_col
               ! set unit
               select case (unit_line(:u_space-1))
               
                  case ("Z_sun")
                     unit_factor(8)=z_sun
               
                  case ("Z/(X+Y+Z)")
                     unit_factor(8)=1._rel_kind
                     
                  case default
                    write(6,"(A)") "Error: unknown time unit "//trim(unit_line(:u_space-1))
                    stop
                     
                end select
                
         end select
         
         ! remove column label from header line
         header=adjustl(header(h_space+1:))
         unit_line=adjustl(unit_line(u_space+1:))
                  
      end do
      
      ! read through rest of file to get number of rows
      n_row=0
      do
      
         read(1,*,iostat=read_status)
         
         if (read_status<0) exit
         
         n_row=n_row+1
         
      end do
      
      
      ! allocate data arrays
      allocate(data(n_col,n_row))
      
      
      ! allocate arrays
      allocate(position(3,n_row))
      
      if (present(luminosity)) allocate(luminosity(n_row))
      if (present(temperature)) allocate(temperature(n_row))
      if (present(mass)) allocate(mass(n_row))
      if (present(age)) allocate(age(n_row))
      if (present(metallicity)) allocate(metallicity(n_row))
      
      ! rewind file
      rewind(1)
      ! skip header
      read(1,*)
      read(1,*)
      
      ! read in data block
      do i=1,n_row
      
         read(1,*) data(:,i)
         
      end do
      
      ! close file
      close(1)
      
      ! assign to arrays
      position(1,:)=data(c_x1,:)*unit_factor(1)
      position(2,:)=data(c_x2,:)*unit_factor(2)
      position(3,:)=data(c_x3,:)*unit_factor(3)
      if (present(luminosity)) luminosity=data(c_l,:)*unit_factor(4)
      if (present(temperature)) temperature=data(c_t,:)*unit_factor(5)
      if (present(mass)) mass=data(c_m,:)*unit_factor(6)
      if (present(age)) age=data(c_a,:)*unit_factor(7)
      if (present(metallicity)) metallicity=data(c_z,:)*unit_factor(8)
      
      
      ! deallocate data block
      deallocate(data)
      
      return
   
   end subroutine

   ! write out sph particles to file
   subroutine write_out_sph_particles_3d&
      &(position,mass,temperature,density,smoothing_length,file_name)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: position(:,:)                         ! particle positions
      real(kind=rel_kind),intent(in) :: mass(:)                               ! particles masses
      real(kind=rel_kind),intent(in) :: temperature(:)                        ! particle temperatures
      real(kind=rel_kind),intent(in) :: density(:)                            ! density
      real(kind=rel_kind),intent(in) :: smoothing_length(:)                   ! smoothing length
      character(kind=chr_kind,len=string_length),intent(in) :: file_name      ! name of input file
      
      ! units(1) = position(1) unit
      ! units(2) = position(2) unit
      ! units(3) = position(3) unit
      ! units(4) = mass unit
      ! units(5) = temperature unit
      ! units(6) = density unit
      ! units(7) = smoothing length unit
      
      ! variable declarations
      integer(kind=int_kind) :: i                                             ! counter
      
      ! open file
      open(1,file=trim(file_name),status="unknown",form="formatted")
      
      ! write header
      write(1,"(7(A27))") "x_1","x_2","x_3","m","T","rho","h"
      
      ! internal cgs units
      write(1,"(7(A25))")&
         &"cm","cm","cm","g","K","g_cm^-3","cm"
            
      ! write out data
      do i=1,size(position,2)
      
         write(1,"(7(E27.17E3))") position(:,i),mass(i),temperature(i),density(i),smoothing_length(i)
      
      end do
      
      ! close file
      close(1)
      
      return
   
   end subroutine
   
   subroutine write_out_datacube_3d(x,y,lambda,i_lambda,sigma,sim_id)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x(:)                         ! x positions
      real(kind=rel_kind),intent(in) :: y(:)                         ! y positions
      real(kind=rel_kind),intent(in) :: lambda(:)                    ! wavelengths
      real(kind=rel_kind),intent(in) :: i_lambda(:,:,:)              ! intensity grid (lambda,x,y)
      real(kind=rel_kind),intent(in) :: sigma(:,:)                   ! column density
      character(kind=chr_kind,len=string_length),intent(in) :: sim_id! run id
      
      ! variable declarations
      real(kind=rel_kind) :: d_a                                     ! pixel size
      character(kind=chr_kind,len=string_length) :: file_name        ! file name
      
      ! variable declarations
      integer(kind=int_kind) :: i,j,k                                ! counter
      
      ! write out image for each wavelength
      do i=1,size(lambda)
      
         write(file_name,"(A,E10.4,A)") trim(sim_id)//"/image_lambda=",lambda(i),"_micron.dat" 
         open(1,file=trim(file_name))
         ! write header 
         write(1,"(3(A27))") "x_1","x_2","i_lambda"
         write(1,"(2(A27),A31)") "cm", "cm","erg s^-1 sr^-1 cm^-2 micron^-1"
         
            do k=1,size(y)
               do j=1,size(x)
                  write(1,"(3(E27.17E3))") x(j),y(k),i_lambda(i,j,k)
               end do
            end do
         
         close(1)
      
      end do
      
      ! write out spectrum
      d_a=((x(size(x))-x(1))*(y(size(y))-y(1)))/real(size(x)*size(y),rel_kind)
      open(1,file=trim(sim_id)//"/spectrum.dat")
      
      write(1,"(2(A27))") "lambda","I_lambda"
      write(1,"(A27,A27)") "micron", "erg s^-1 sr^-1 micron^-1"
      
      do i=1,size(lambda)
      
         write(1,"(2(E27.17E3))") lambda(i),sum(i_lambda(i,:,:))*d_a
      end do
      
      close(1)
      
      open(1,file=trim(sim_id)//"/column_density.dat")
      write(1,"(3(A27))") "x_1","x_2","sigma"
      write(1,"(3(A27))") "cm", "cm","g cm^-2"
         do j=1,size(y)
            do i=1,size(x)
               write(1,"(3(E27.17E3))") x(i),y(j),sigma(i,j)
            end do
         end do
      close(1)
      
      return
   
   end subroutine

end module