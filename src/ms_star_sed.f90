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

! main sequence sed module
module m_ms_star_sed

   use m_kind_parameters
   use m_maths
   
   implicit none
   
   private
   public :: ms_star_sed
   
   ! main sequence star SED class
   type :: ms_star_sed
   
      real(kind=rel_kind) :: bol_luminosity                          ! bolometric luminosity [erg s^-1]
      real(kind=rel_kind),allocatable :: wavelength_array(:)         ! array of wavelengths [micron]
      real(kind=rel_kind),allocatable :: intensity_array(:)          ! stellar surface intensity [erg s^-1 sr^-1 cm^-2 micron^-1]
      real(kind=rel_kind),allocatable :: mono_luminosity(:)          ! monochromatic luminosity [erg s^-1 micron^-1]
      
      contains
      
      procedure,non_overridable :: initialise
      procedure,non_overridable :: destroy
   
   end type
   
   contains
   
   ! read in files and set SED
   subroutine initialise(self,mass,age,metallicity)
   
      ! argument declarations
      class(ms_star_sed),intent(inout) :: self                       ! SED object
      real(kind=rel_kind),intent(in) :: mass                         ! mass of star [M_sun]
      real(kind=rel_kind),intent(in) :: age                          ! age of star [years]
      real(kind=rel_kind),intent(in) :: metallicity                  ! metallicity of star [Z value]
      
      ! variable declarations
      integer(kind=int_kind) :: i                                    ! counter
      integer(kind=int_kind) :: i_0,i_1                              ! indicies
      integer(kind=int_kind) :: j_0,j_1                              ! indicies
      real(kind=rel_kind) :: log_age                                 ! log10 age
      real(kind=rel_kind) :: log_metallicity                         ! log10 metallicity
      real(kind=rel_kind) :: luminosity                              ! interpolated luminosity
      real(kind=rel_kind) :: temperature                             ! interpolated temperature
      real(kind=rel_kind) :: surface_gravity                         ! interpolated surface gravity
      real(kind=rel_kind) :: luminosity_00,luminosity_01             ! luminosity interpolation points
      real(kind=rel_kind) :: luminosity_10,luminosity_11             ! luminosity interpolation points
      real(kind=rel_kind) :: temperature_00,temperature_01           ! temperature interpolation points
      real(kind=rel_kind) :: temperature_10,temperature_11           ! temperature interpolation points
      real(kind=rel_kind) :: surface_gravity_00,surface_gravity_01   ! surface gravity interpolation points
      real(kind=rel_kind) :: surface_gravity_10,surface_gravity_11   ! surface gravity interpolation points
      real(kind=rel_kind) :: metallicity_value(7)                    ! metallicity values
      real(kind=rel_kind) :: age_value(71)                           ! age values
      real(kind=rel_kind) :: surface_gravity_value(11)               ! surface gravity values
      real(kind=rel_kind) :: temperature_value(75)                   ! temperature values
      real(kind=rel_kind),allocatable :: lambda_0(:),lambda_1(:)        ! wavelength array
      real(kind=rel_kind),allocatable :: lambda_00(:),lambda_01(:)      ! wavelength array
      real(kind=rel_kind),allocatable :: lambda_10(:),lambda_11(:)      ! wavelength array
      real(kind=rel_kind),allocatable :: f_lambda_0(:),f_lambda_1(:)    ! flux interpolation points
      real(kind=rel_kind),allocatable :: f_lambda_00(:),f_lambda_01(:)  ! flux interpolation points
      real(kind=rel_kind),allocatable :: f_lambda_10(:),f_lambda_11(:)  ! flux interpolation points
      character(kind=chr_kind,len=string_length) :: metallicity_label_1(7)    ! metallicity labels (luminosity files)
      character(kind=chr_kind,len=string_length) :: metallicity_label_2(7)    ! metallicity labels (sed files)
      character(kind=chr_kind,len=string_length) :: age_label(71)             ! age labels
      character(kind=chr_kind,len=string_length) :: surface_gravity_label(11) ! surface gravity labels
      character(kind=chr_kind,len=string_length) :: temperature_label(75)     ! temperature labels
      character(kind=chr_kind,len=string_length) :: file_name        ! name of file
      
      ! set metallicity labels
      metallicity_label_1=(/"-2.00","-1.50","-1.00","-0.50","+0.00","+0.20","+0.50"/)
      metallicity_label_2=(/"m20","m15","m10","m05","p00","p02","p05"/)
      
      ! read values from labels
      do i=1,size(metallicity_label_1)
         read(metallicity_label_1(i),*) metallicity_value(i)
      end do
      metallicity_value=log10(0.0134_rel_kind)+metallicity_value

      ! set age labels
      age_label=(/"3.9811E+06","4.4668E+06","5.0119E+06","5.6234E+06","6.3096E+06","7.0795E+06","7.9433E+06","8.9125E+06",&
         &"1.0000E+07","1.1220E+07","1.2589E+07","1.4125E+07","1.5849E+07","1.7783E+07","1.9953E+07","2.2387E+07","2.5119E+07",&
         &"2.8184E+07","3.1623E+07","3.5482E+07","3.9811E+07","4.4669E+07","5.0119E+07","5.6235E+07","6.3096E+07","7.0795E+07",&
         &"7.9434E+07","8.9126E+07","1.0000E+08","1.1220E+08","1.2589E+08","1.4126E+08","1.5849E+08","1.7783E+08","1.9953E+08",&
         &"2.2388E+08","2.5119E+08","2.8184E+08","3.1623E+08","3.5482E+08","3.9811E+08","4.4669E+08","5.0120E+08","5.6235E+08",&
         &"6.3097E+08","7.0796E+08","7.9434E+08","8.9127E+08","1.0000E+09","1.1220E+09","1.2590E+09","1.4126E+09","1.5849E+09",&
         &"1.7783E+09","1.9953E+09","2.2388E+09","2.5119E+09","2.8184E+09","3.1624E+09","3.5482E+09","3.9812E+09","4.4669E+09",&
         &"5.0120E+09","5.6236E+09","6.3097E+09","7.0797E+09","7.9435E+09","8.9128E+09","1.0000E+10","1.1221E+10","1.2590E+10"/)
         
      ! read values from labels
      do i=1,size(age_label)
         read(age_label(i),*) age_value(i)
      end do
      age_value=log10(age_value)
      
      ! set surface gravity labels
      surface_gravity_label=(/"0.00000","0.50000","1.00000","1.50000",&
         &"2.00000","2.50000","3.00000","3.50000","4.00000","4.50000","5.00000"/)
         
      do i=1,size(surface_gravity_label)
         read(surface_gravity_label(i),*) surface_gravity_value(i)
      end do
      
      temperature_label=(/"3500. ","3750. ","4000. ","4250. ","4500. ","4750. ","5000. ","5250. ","5500. ","5750. ","6000. ",&
         &"6250. ","6500. ","6750. ","7000. ","7250. ","7500. ","7750. ","8000. ","8250. ","8500. ","8750. ","9000. ","9250. ",&
         &"9500. ","9750. ","10000.","10250.","10500.","10750.","11000.","11250.","11500.","11750.","12000.","12250.","12500.",&
         &"12750.","13000.","14000.","15000.","16000.","17000.","18000.","19000.","20000.","21000.","22000.","23000.","24000.",&
         &"25000.","26000.","27000.","28000.","29000.","30000.","31000.","32000.","33000.","34000.","35000.","36000.","37000.",&
         &"38000.","39000.","40000.","41000.","42000.","43000.","44000.","45000.","46000.","47000.","48000.","49000."/)      
      
      do i=1,size(temperature_label)
         read(temperature_label(i),*) temperature_value(i)
      end do
      temperature_value=log10(temperature_value)
      
      
      ! Interpolate luminosity, temperature and surface gravity from lookup tables
      
      ! get log values (make sure they are in range if necessary)
      log_metallicity=max(min(log10(metallicity),metallicity_value(size(metallicity_value))),metallicity_value(1))
      log_age=max(min(log10(age),age_value(size(age_value))),age_value(1))
      
      ! get indicies
      i_0=binary_search(log_metallicity,metallicity_value)
      i_1=i_0+1
      j_0=binary_search(log_age,age_value)
      j_1=j_0+1
      
      
      ! read in first set of interpolation points
      file_name="msstardata/luminosity/logZ="//trim(metallicity_label_1(i_0))//"_t="//trim(age_label(j_0))//".dat"
      call get_luminosity_interpolation_points(mass,luminosity_00,temperature_00,surface_gravity_00,file_name)
            
      ! read in second set of interpolation points
      file_name="msstardata/luminosity/logZ="//trim(metallicity_label_1(i_1))//"_t="//trim(age_label(j_0))//".dat"
      call get_luminosity_interpolation_points(mass,luminosity_01,temperature_01,surface_gravity_01,file_name)

      ! read in third set of interpolation points
      file_name="msstardata/luminosity/logZ="//trim(metallicity_label_1(i_0))//"_t="//trim(age_label(j_1))//".dat"
      call get_luminosity_interpolation_points(mass,luminosity_10,temperature_10,surface_gravity_10,file_name)
      
      ! read in fourth set of interpolation points
      file_name="msstardata/luminosity/logZ="//trim(metallicity_label_1(i_1))//"_t="//trim(age_label(j_1))//".dat"
      call get_luminosity_interpolation_points(mass,luminosity_11,temperature_11,surface_gravity_11,file_name)
      
      
      ! interpolate out age and metallicity     
      luminosity=lerp(log_age,log_metallicity,age_value(j_0),age_value(j_1),metallicity_value(i_0),metallicity_value(i_1),&
         &luminosity_00,luminosity_01,luminosity_10,luminosity_11)
      temperature=lerp(log_age,log_metallicity,age_value(j_0),age_value(j_1),metallicity_value(i_0),metallicity_value(i_1),&
         &temperature_00,temperature_01,temperature_10,temperature_11)
      surface_gravity=lerp(log_age,log_metallicity,age_value(j_0),age_value(j_1),metallicity_value(i_0),metallicity_value(i_1),&
         &surface_gravity_00,surface_gravity_01,surface_gravity_10,surface_gravity_11)
            
      self%bol_luminosity=10._rel_kind**luminosity
      
      ! Interpolate SED from lookup table
      
      ! get indicies
      j_0=binary_search(temperature,temperature_value)
      j_1=j_0+1

      ! get first SED
      call get_sed_interpolation_points(surface_gravity,surface_gravity_value,&
         &metallicity_label_2(i_0),temperature_label(j_0),surface_gravity_label,lambda_00,f_lambda_00)
      
      ! get second SED
      call get_sed_interpolation_points(surface_gravity,surface_gravity_value,&
         &metallicity_label_2(i_0),temperature_label(j_1),surface_gravity_label,lambda_01,f_lambda_01)
      
      ! get third SED
      call get_sed_interpolation_points(surface_gravity,surface_gravity_value,&
         &metallicity_label_2(i_1),temperature_label(j_0),surface_gravity_label,lambda_10,f_lambda_10)
      
      ! get fourth SED
      call get_sed_interpolation_points(surface_gravity,surface_gravity_value,&
         &metallicity_label_2(i_1),temperature_label(j_1),surface_gravity_label,lambda_11,f_lambda_11)
      
      ! interpolate out temperature dependence
      call interpolate_sed(temperature,temperature_value(j_0),temperature_value(j_1),&
         &lambda_00,lambda_01,f_lambda_00,f_lambda_01,lambda_0,f_lambda_0)
      
      call interpolate_sed(temperature,temperature_value(j_0),temperature_value(j_1),&
         &lambda_10,lambda_11,f_lambda_10,f_lambda_11,lambda_1,f_lambda_1)
         
      ! interpolate out metallicity dependence
      call interpolate_sed(log_metallicity,metallicity_value(i_0),metallicity_value(i_1),&
         &lambda_0,lambda_1,f_lambda_0,f_lambda_1,self%wavelength_array,self%intensity_array)
      
      ! convert wavelength [angstrom] to [micron]
      self%wavelength_array=self%wavelength_array*1.e-4_rel_kind
      
      ! convert flux [erg s^-1 cm^-2 angstrom^-1] to intensity [erg s^-1 sr^-1 cm^-2 micron^-1]
      self%intensity_array=self%intensity_array*1.e+4_rel_kind/pi
      
      ! convert luminosity [L_sun] to [erg s^-1]
      self%bol_luminosity=self%bol_luminosity*3.846e+33_rel_kind
      
      ! calculate monochromatic luminosity
      self%mono_luminosity=self%bol_luminosity*self%intensity_array/&
         &trapz_intgr(self%wavelength_array,self%intensity_array)
 
      return
   
   end subroutine
   
   ! deallocate array in ms_star_sed object
   pure subroutine destroy(self)
   
      ! argument declarations
      class(ms_star_sed),intent(inout) :: self                                ! SED object
      
      deallocate(self%wavelength_array)
      deallocate(self%intensity_array)
      deallocate(self%mono_luminosity)
      
      return
   
   end subroutine
   
   ! interpolated luminosity temperature and surface gravity from file
   subroutine get_luminosity_interpolation_points(mass,luminosity,temperature,surface_gravity,file_name)
      
      ! variable declarations
      real(kind=rel_kind),intent(in) :: mass                               ! mass
      real(kind=rel_kind),intent(out) :: luminosity                        ! luminosity
      real(kind=rel_kind),intent(out) :: temperature                       ! temperature
      real(kind=rel_kind),intent(out) :: surface_gravity                   ! surface gravity
      character(kind=chr_kind,len=string_length),intent(in) :: file_name   ! name of file
           
      ! variable declarations
      integer(kind=int_kind) :: i                                       ! counter
      integer(kind=int_kind) :: n_lines                                 ! number of lines in file
      integer(kind=int_kind) :: read_status                             ! file iostat variable
      real(kind=rel_kind) :: mass_0,mass_1                              ! mass interpolation points
      real(kind=rel_kind) :: luminosity_0,luminosity_1                  ! luminosity interpolation points
      real(kind=rel_kind) :: temperature_0,temperature_1                ! temperature interpolation points
      real(kind=rel_kind) :: surface_gravity_0,surface_gravity_1        ! surface gravity interpolation points
      real(kind=rel_kind) :: dummy_real                                 ! dummy real variable
                           
      open(1,file=trim(file_name))
      n_lines=0
      do                                        ! get number of lines in file
         read(1,*,iostat=read_status)
         if (read_status<0) exit
         n_lines=n_lines+1
      end do
      rewind(1)
      read(1,*)
      read(1,*)
      read(1,*) dummy_real,dummy_real,mass_0,dummy_real,luminosity_0,temperature_0,surface_gravity_0
      read(1,*) dummy_real,dummy_real,mass_1,dummy_real,luminosity_1,temperature_1,surface_gravity_1
      do i=1,n_lines-4
         if (mass_1>mass) exit                  ! found mass interpolation points
         mass_0=mass_1                          ! set lower bounds and read in new upper bounds
         luminosity_0=luminosity_1
         temperature_0=temperature_1
         surface_gravity_0=surface_gravity_1
         read(1,*) dummy_real,dummy_real,mass_1,dummy_real,luminosity_1,temperature_1,surface_gravity_1
      end do
      close(1)
      
      luminosity=lerp_forced(log(mass),log(mass_0),log(mass_1),luminosity_0,luminosity_1)
      temperature=lerp_forced(log(mass),log(mass_0),log(mass_1),temperature_0,temperature_1)
      surface_gravity=lerp_forced(log(mass),log(mass_0),log(mass_1),surface_gravity_0,surface_gravity_1)
      
      return
      
   end subroutine
   
   ! interpolated sed from file
   subroutine get_sed_interpolation_points(surface_gravity,surface_gravity_value,&
      &metallicity_label,temperature_label,surface_gravity_label,lambda,f_lambda)
      
      ! argument declarations
      real(kind=rel_kind),intent(in) :: surface_gravity                                      ! surface gravity
      real(kind=rel_kind),intent(in) :: surface_gravity_value(:)                             ! available surface gravity values
      character(kind=chr_kind,len=string_length),intent(in) :: metallicity_label             ! metallicity label
      character(kind=chr_kind,len=string_length),intent(in) :: temperature_label             ! temperature label
      character(kind=chr_kind,len=string_length),intent(in) :: surface_gravity_label(:)      ! surface gravity label
      real(kind=rel_kind),intent(out),allocatable :: lambda(:)                               ! wavelength array
      real(kind=rel_kind),intent(out),allocatable :: f_lambda(:)                             ! flux array
      
      ! variable declarations
      integer(kind=int_kind) :: i                                       ! counter
      integer(kind=int_kind) :: i_0,i_1                                 ! indicies
      integer(kind=int_kind) :: n_lines                                 ! number of lines in file
      integer(kind=int_kind) :: read_status                             ! file iostat variable
      real(kind=rel_kind),allocatable :: lambda_0(:)                    ! wavelength array
      real(kind=rel_kind),allocatable :: lambda_1(:)                    ! wavelength array
      real(kind=rel_kind),allocatable :: f_lambda_0(:)                  ! flux interpolation points
      real(kind=rel_kind),allocatable :: f_lambda_1(:)                  ! flux interpolation points
      real(kind=rel_kind),allocatable :: avail_surface_gravity_value(:) ! available surface gravity values
      logical(kind=log_kind) :: surface_gravity_mask(size(surface_gravity_label))   ! available surface gravity files
      character(kind=chr_kind,len=string_length),allocatable :: avail_surface_gravity_label(:)  ! available surface gravity labels
      character(kind=chr_kind,len=string_length) :: file_name           ! name of file
      
      ! get available surface gravity values
      do i=1,size(surface_gravity_label)      
         file_name="msstardata/sed/f"//trim(metallicity_label)//"k2odfnew.pck.teff="//trim(temperature_label)//".logg="//&
         &trim(surface_gravity_label(i))//".dat.txt"
         inquire(file=trim(file_name),exist=surface_gravity_mask(i))
      end do
      
      avail_surface_gravity_value=pack(surface_gravity_value,surface_gravity_mask)
      avail_surface_gravity_label=pack(surface_gravity_label,surface_gravity_mask)
      
      i_0=binary_search(surface_gravity,surface_gravity_value)
      i_1=i_0+1
      
      ! read in values from first file
      file_name="msstardata/sed/f"//trim(metallicity_label)//"k2odfnew.pck.teff="//trim(temperature_label)//".logg="//&
         &trim(surface_gravity_label(i_0))//".dat.txt"
      open(1,file=trim(file_name))
      n_lines=0
      do                                        ! get number of lines in file
         read(1,*,iostat=read_status)
         if (read_status<0) exit
         n_lines=n_lines+1
      end do
      rewind(1)
      allocate(lambda_0(n_lines-3))
      allocate(f_lambda_0(n_lines-3))
      read(1,*)                                 ! ignore first three lines
      read(1,*)
      read(1,*)
      do i=1,n_lines-3
         read(1,*) lambda_0(i),f_lambda_0(i)
      end do
      close(1)
      
      ! read in values from second file
      file_name="msstardata/sed/f"//trim(metallicity_label)//"k2odfnew.pck.teff="//trim(temperature_label)//".logg="//&
         &trim(surface_gravity_label(i_1))//".dat.txt"
      open(1,file=trim(file_name))
      n_lines=0
      do                                        ! get number of lines in file
         read(1,*,iostat=read_status)
         if (read_status<0) exit
         n_lines=n_lines+1
      end do
      rewind(1)
      allocate(lambda_1(n_lines-3))
      allocate(f_lambda_1(n_lines-3))
      read(1,*)                                 ! ignore first three lines
      read(1,*)
      read(1,*)
      do i=1,n_lines-3
         read(1,*) lambda_1(i),f_lambda_1(i)
      end do
      close(1)
      
      ! interpolate two SEDs
      call interpolate_sed(surface_gravity,surface_gravity_value(i_0),surface_gravity_value(i_1),&
         &lambda_0,lambda_1,f_lambda_0,f_lambda_1,lambda,f_lambda)
      
      return
      
   end subroutine
   
   ! interpolate together two SEDs
   pure subroutine interpolate_sed(x,x_0,x_1,lambda_0,lambda_1,f_lambda_0,f_lambda_1,lambda,f_lambda)
   
      ! argument declarations
      real(kind=rel_kind),intent(in) :: x                            ! x value
      real(kind=rel_kind),intent(in) :: x_0,x_1                      ! x interpolation points
      real(kind=rel_kind),intent(in) :: lambda_0(:),lambda_1(:)      ! wavelength arrays
      real(kind=rel_kind),intent(in) :: f_lambda_0(:),f_lambda_1(:)  ! flux interpolation points
      real(kind=rel_kind),intent(out),allocatable :: lambda(:)       ! output wavelength array
      real(kind=rel_kind),intent(out),allocatable :: f_lambda(:)     ! output flux array
      
      ! variable declarations
      real(kind=rel_kind) :: f_lambda_0_temp(max(size(f_lambda_0),size(f_lambda_1)))   ! temp flux arrays
      real(kind=rel_kind) :: f_lambda_1_temp(max(size(f_lambda_0),size(f_lambda_1)))   ! temp flux arrays
      
      ! allocate arrays
      allocate(lambda(max(size(lambda_0),size(lambda_1)))) 
      allocate(f_lambda(max(size(f_lambda_0),size(f_lambda_1))))
      
      ! set f_lambda arrays to (practically) zero
      f_lambda_0_temp=minval(f_lambda_0)
      f_lambda_1_temp=minval(f_lambda_1)
      
      ! copy arrays      
      f_lambda_0_temp(max(size(f_lambda_1)-size(f_lambda_0),0)+1:)=f_lambda_0
      f_lambda_1_temp(max(size(f_lambda_0)-size(f_lambda_1),0)+1:)=f_lambda_1
      if (size(lambda_0)<size(lambda_1)) then
         lambda=lambda_1
      else
         lambda=lambda_0
      end if
      
      ! interpolate arrays
      f_lambda=exp(lerp_forced(x,x_0,x_1,log(f_lambda_0_temp),log(f_lambda_1_temp)))
      
      return
   
   end subroutine

end module