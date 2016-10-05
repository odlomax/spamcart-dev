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

! options module for simulation parameters
! see Rotitaille 2010, A&A, 520, 70 for Modified Random Walk details
module m_options

   use m_kind_parameters
   
   implicit none
   
   private
   public :: options
   
   ! define simulation options class
   type :: options

      ! simulation params
      character(kind=chr_kind,len=string_length) :: sim_cloud_file            ! cloud input file
      character(kind=chr_kind,len=string_length) :: sim_star_file             ! stars input file
      character(kind=chr_kind,len=string_length) :: sim_id                    ! simulation run id
      integer(kind=int_kind) :: sim_n_packet_point                            ! number of point source luminosity packets per iteration
      integer(kind=int_kind) :: sim_n_packet_external                         ! number of external radiation field luminosity packets per iteration
      logical(kind=log_kind) :: sim_equal_packets_per_point                   ! equal number packets per point source
      integer(kind=int_kind) :: sim_n_it                                      ! number of iterations
      logical(kind=log_kind) :: sim_restart                                   ! restart sim? (don't use sim_initial_t)
      real(kind=rel_kind) :: sim_initial_t                                    ! initial background temperature
      real(kind=rel_kind) :: sim_min_d                                        ! simulation minimum resolvable distance
      logical(kind=log_kind) :: sim_mrw                                       ! use modified random walk for optically thick regions
      real(kind=rel_kind) :: sim_mrw_gamma                                    ! mrw gamma variable
      integer(kind=int_kind) :: sim_n_mrw                                     ! number of mrw lookup table points
      
      ! dust params
      real(kind=rel_kind) :: dust_t_min                                       ! minimum dust temperature
      real(kind=rel_kind) :: dust_t_max                                       ! maximum dust temperature
      integer(kind=int_kind) :: dust_n_t                                      ! number of dust temperature samples
      real(kind=rel_kind) :: dust_r_v                                         ! dust extinction cure coeff
      logical(kind=log_kind) :: dust_iso_scatter                              ! isotropic scattering
      real(kind=rel_kind) :: dust_sub_t_min                                   ! minimum dust sublimation temperature
      real(kind=rel_kind) :: dust_sub_t_max                                   ! maximum dust sublimation temperature
      
      ! external radiation field params
      logical(kind=log_kind) :: ext_rf                                        ! include an external radiation field
      character(kind=chr_kind,len=string_length) :: ext_rf_type               ! external radiation field type (bb or ps05)
      logical(kind=log_kind) :: ext_rf_cmb                                    ! add a CMB ro radiation field
      real(kind=rel_kind) :: ext_rf_gal_r                                     ! galactic radius (kiloparsecs)
      real(kind=rel_kind) :: ext_rf_gal_z                                     ! galactic height (kiloparsecs)
      real(kind=rel_kind) :: ext_rf_t                                         ! blackbody temperature
      real(kind=rel_kind) :: ext_rf_d                                         ! blackbody dilution factor
      real(kind=rel_kind) :: ext_rf_lambda_min                                ! blackbody min wavelength
      real(kind=rel_kind) :: ext_rf_lambda_max                                ! blackbody max wavelength
      integer(kind=int_kind) :: ext_rf_n_lambda                               ! number of wavelength samples
      
      ! point source params
      logical(kind=log_kind) :: point_sources                                 ! include point sources
      real(kind=rel_kind) :: point_lambda_min                                 ! minimum wavelength
      real(kind=rel_kind) :: point_lambda_max                                 ! maximum wavelength
      integer(kind=int_kind) :: point_n_lambda                                ! number of wavelengths
      
      ! sph params
      character(kind=chr_kind,len=string_length) :: sph_kernel                ! sph kernel type
      real(kind=rel_kind) :: sph_eta                                          ! smoothing length scale factor
      logical(kind=log_kind) :: sph_scattered_light                           ! capture scattered light for each particle
      real(kind=rel_kind) :: sph_lambda_min                                   ! minimum scattered light wavelength
      real(kind=rel_kind) :: sph_lambda_max                                   ! maximum scattered light wavelength
      integer(kind=int_kind) :: sph_n_lambda                                  ! number of wavelengths
      
      ! datacube params
      logical(kind=log_kind) :: datacube_make                                 ! make a datacube
      logical(kind=log_kind) :: datacube_convolve                             ! convolve datacube with psf
      character(kind=chr_kind,len=string_length) :: datacube_lambda_string    ! string of wavelenghts for data cube
      character(kind=chr_kind,len=string_length) :: datacube_fwhm_string      ! string of beam fwhm
      real(kind=rel_kind) :: datacube_distance                                ! distance from source to observer
      real(kind=rel_kind) :: datacube_x_min                                   ! minimum x coordinate     
      real(kind=rel_kind) :: datacube_x_max                                   ! maximum x coordinate
      real(kind=rel_kind) :: datacube_y_min                                   ! minimum y coordinate
      real(kind=rel_kind) :: datacube_y_max                                   ! maximum y coordinate
      real(kind=rel_kind) :: datacube_angles(3)                               ! altitude, azimuth and image rotation angles
      integer(kind=int_kind) :: datacube_n_x                                  ! number of pixels along x axis
      integer(kind=int_kind) :: datacube_n_y                                  ! number of pixels along y axis
      
      contains
      
      procedure,non_overridable :: initialise

   end type
   
   contains
   
   ! set defaults and read in params file
   subroutine initialise(self,params_file)
   
      ! argument declarations 
      class(options),intent(inout) :: self                                          ! options class
      character(kind=chr_kind,len=string_length),intent(in) :: params_file          ! parameters file name
      
      ! variable declaration
      integer(kind=int_kind) :: read_status                                         ! read status variable
      integer(kind=int_kind) :: i_hash                                              ! index of hash character
      integer(kind=int_kind) :: i_equal                                             ! inde of equal sign character
      character(kind=chr_kind,len=string_length) :: file_name                       ! file name
      character(kind=chr_kind,len=string_length) :: params_file_line                ! line from parameters file
      character(kind=chr_kind,len=string_length) :: format_string                   ! format string for params_file_line
      character(kind=chr_kind,len=string_length) :: param_name                      ! name of parameter
      character(kind=chr_kind,len=string_length) :: param_value                     ! value of parameter
      
      ! set defaults
      
      ! simulation defaults
      self%sim_cloud_file=""
      self%sim_star_file=""
      self%sim_id=""
      self%sim_n_packet_point=1000000
      self%sim_n_packet_external=1000000
      self%sim_equal_packets_per_point=.false.
      self%sim_n_it=5
      self%sim_restart=.false.
      self%sim_initial_t=10._rel_kind
      self%sim_min_d=0._rel_kind
      self%sim_mrw=.true.
      self%sim_mrw_gamma=2._rel_kind
      self%sim_n_mrw=1001
      
      
      ! dust table parameters
      self%dust_t_min=2.e+0_rel_kind
      self%dust_t_max=2.e+4_rel_kind
      self%dust_n_t=801
      self%dust_r_v=5.5_rel_kind
      self%dust_iso_scatter=.true.
      self%dust_sub_t_min=900._rel_kind
      self%dust_sub_t_max=1100._rel_kind
      
      ! external radiation field parameters
      self%ext_rf=.true.
      self%ext_rf_type="ps05"
      self%ext_rf_cmb=.true.
      self%ext_rf_gal_r=8.5_rel_kind
      self%ext_rf_gal_z=0._rel_kind
      self%ext_rf_t=1.e4_rel_kind
      self%ext_rf_d=5.e-16
      self%ext_rf_lambda_min=1.e-1_rel_kind
      self%ext_rf_lambda_max=1.e+4_rel_kind
      self%ext_rf_n_lambda=1001
      
      ! point source parameters
      self%point_sources=.true.
      self%point_lambda_min=1.e-1_rel_kind
      self%point_lambda_max=1.e+4_rel_kind
      self%point_n_lambda=1001
      
      ! sph parameteres
      self%sph_kernel="m4"
      self%sph_eta=1.2_rel_kind
      self%sph_scattered_light=.false.
      self%sph_lambda_min=1.e-1_rel_kind
      self%sph_lambda_max=1.e+1_rel_kind
      self%sph_n_lambda=11
      
      ! set datacube  params
      self%datacube_make=.true.
      self%datacube_convolve=.true.
      self%datacube_lambda_string="3.6 4.5 5.8 8.0 24. 70. 100. 160. 250. 350. 500."    ! IRAC/MIPS/PACS/SPIRE (micron)
      self%datacube_fwhm_string="1.7 1.7 1.7 1.9 6.0 5.6 6.8 12.0 17.6 23.9 35.2"       ! PSF FWHM (arcsec)
      self%datacube_distance=3.0856776e+18                                          ! distance from source (1 pc)
      self%datacube_x_min=-4.4879361e+16_rel_kind                                   ! 3000 au
      self%datacube_x_max=4.4879361e+16_rel_kind                                    ! 3000 au
      self%datacube_y_min=-4.4879361e+16_rel_kind                                   ! 3000 au
      self%datacube_y_max=4.4879361e+16_rel_kind                                    ! 3000 au
      self%datacube_angles=(/0._rel_kind,0._rel_kind,0._rel_kind/)                  ! x-y projection
      self%datacube_n_x=256
      self%datacube_n_y=256
      
      
      ! read in parameters from file
      if (len(trim(params_file))>0) then
         file_name=params_file
      else
         file_name="params.dat"
      end if
      
      ! open parameters file
      open(1,file=trim(file_name),status="old",form="formatted")
      
      ! set format string
      write(format_string,"(A2,I0,A1)") "(A",string_length,")"
      
      ! read file line by line
      do
      
         ! read line
         read(1,trim(format_string),iostat=read_status) params_file_line
         
         ! have we reached end of file?
         if (read_status<0) exit
      
         ! was there an error
         if (read_status>0) then
            write(6,"(A)") "Error reading params file."
            stop
         end if
         
         ! remove leading white space
         params_file_line=adjustl(params_file_line)
         
         ! cycle if line is comment line
         if (params_file_line(1:1)=="#") cycle
         
         ! look for "=" character in line
         i_equal=index(params_file_line,"=")
         
         if (i_equal>0) then
            
            ! strip out comment if present
            i_hash=index(params_file_line,"#")
            if (i_hash>0) then
            
               params_file_line=params_file_line(:i_hash-1)
            
            end if
            
            param_name=adjustl(params_file_line(:i_equal-1))
            param_value=adjustl(params_file_line(i_equal+1:))
            
            select case (trim(param_name))
            
               case ("sim_cloud_file")
                  self%sim_cloud_file=param_value
               
               case ("sim_star_file")
                  self%sim_star_file=param_value
                  
               case ("sim_id")
                  self%sim_id=param_value
               
               case ("sim_n_packet_point")
                  read(param_value,*) self%sim_n_packet_point
                  
               case ("sim_n_packet_external")
                  read(param_value,*) self%sim_n_packet_external
                  
               case ("sim_equal_packets_per_point")
                  read(param_value,*) self%sim_equal_packets_per_point
                  
               case ("sim_n_it")
                  read(param_value,*) self%sim_n_it
               
               case ("sim_restart")
                  read(param_value,*) self%sim_restart
                  
               case ("sim_initial_t")
                  read(param_value,*) self%sim_initial_t
                  
               case ("sim_min_d")
                  read(param_value,*) self%sim_min_d
                  
               case ("sim_mrw")
                  read(param_value,*) self%sim_mrw
                  
               case ("sim_mrw_gamma")
                  read(param_value,*) self%sim_mrw_gamma
                  
               case ("sim_n_mrw")
                  read(param_value,*) self%sim_n_mrw
               
               case ("dust_t_min")
                  read(param_value,*) self%dust_t_min
                  
               case ("dust_t_max")
                  read(param_value,*) self%dust_t_max
                  
               case ("dust_r_v")
                  read(param_value,*) self%dust_r_v
               
               case ("dust_iso_scatter")
                  read(param_value,*) self%dust_iso_scatter
                  
               case ("dust_sub_t_min")
                  read(param_value,*) self%dust_sub_t_min
                  
               case ("dust_sub_t_max")
                  read(param_value,*) self%dust_sub_t_max
                  
               case ("ext_rf")
                  read(param_value,*) self%ext_rf
                  
               case ("ext_rf_type")
                  self%ext_rf_type=param_value
                  
               case ("ext_rf_cmb")
                  read(param_value,*) self%ext_rf_cmb
               
               case ("ext_rf_gal_r")
                  read(param_value,*) self%ext_rf_gal_r
                  
               case ("ext_rf_gal_z")
                  read(param_value,*) self%ext_rf_gal_z
                  
               case ("ext_rf_t")
                  read(param_value,*) self%ext_rf_t
                  
               case ("ext_rf_d")
                  read(param_value,*) self%ext_rf_d
                  
               case ("ext_rf_lambda_min")
                  read(param_value,*) self%ext_rf_lambda_min
                  
               case ("ext_rf_lambda_max")
                  read(param_value,*) self%ext_rf_lambda_max
                  
               case ("ext_rf_n_lambda")
                  read(param_value,*) self%ext_rf_n_lambda
                                    
               case ("point_sources")
                  read(param_value,*) self%point_sources
               
               case ("point_lambda_min")
                  read(param_value,*) self%point_lambda_min
                  
               case ("point_lambda_max")
                  read(param_value,*) self%point_lambda_max
                  
               case ("point_n_lambda")
                  read(param_value,*) self%point_n_lambda
                  
               case ("sph_kernel")
                  self%sph_kernel=param_value
                  
               case ("sph_eta")
                  read(param_value,*) self%sph_eta
                  
               case ("sph_scattered_light")
                  read(param_value,*) self%sph_scattered_light
                  
               case ("sph_lambda_min")
                  read(param_value,*) self%sph_lambda_min
                  
               case ("sph_lambda_max")
                  read(param_value,*) self%sph_lambda_max
                  
               case ("sph_n_lambda")
                  read(param_value,*) self%sph_n_lambda
               
               case ("datacube_make")
                  read(param_value,*) self%datacube_make
                  
               case ("datacube_convolve")
                  read(param_value,*) self%datacube_convolve
                                  
               case ("datacube_lambda_string")
                  self%datacube_lambda_string=param_value
                  
               case ("datacube_fwhm_string")
                  self%datacube_fwhm_string=param_value
                  
               case ("datacube_distance")
                  read(param_value,*) self%datacube_distance
                  
               case ("datacube_x_min")
                  read(param_value,*) self%datacube_x_min
                  
               case ("datacube_x_max")
                  read(param_value,*) self%datacube_x_max
                  
               case ("datacube_y_min")
                  read(param_value,*) self%datacube_y_min
                  
               case ("datacube_y_max")
                  read(param_value,*) self%datacube_y_max
                  
               case ("datacube_angles")
                  read(param_value,*) self%datacube_angles
                  
               case ("datacube_n_x")
                  read(param_value,*) self%datacube_n_x
                  
               case ("datacube_n_y")
                  read(param_value,*) self%datacube_n_y
                  
               case default
               
                  write(6,"(A)") 'Warning: unknown parameter "'&
                  &//trim(adjustl(params_file_line(:i_equal-1)))//'"'
               
            end select
         
         end if
      
      end do
      
      ! sort out conflicts
      if (self%sim_n_it==0) self%sph_scattered_light=.false.
      
      ! make sure directory exists for file output
      call execute_command_line("mkdir "//trim(self%sim_id))
      
      ! close params file
      close(1)
   
   end subroutine
   
end module