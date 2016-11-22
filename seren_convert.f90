module sink_mod

   implicit none
   
   ! set parameters
   integer,parameter :: dp=selected_real_kind(p=15)     ! double precision real
   integer,parameter :: sp=selected_real_kind(p=6)      ! single precision real
   integer,parameter :: pr=selected_real_kind(p=15)     ! double precision real
   integer,parameter :: ilp=selected_int_kind(r=15)     ! long integer
   integer,parameter :: DMDT_RANGE=32
   integer,parameter :: NDIM=3
   integer,parameter :: VDIM=3
   integer,parameter :: BDIM=3
   
   contains
   
   SUBROUTINE read_data_seren_unform&
      &(out_file,sink_position,sink_luminosity,sink_temperature,gas_position,gas_mass,gas_temperature,time)
   
!    use particle_module
!    use hydro_module
!    use time_module
!    use type_module
!    use sink_module
!    use mhd_module

   character(len=*), intent(in) :: out_file    ! Name of output file
   real(kind=dp),intent(out),allocatable :: sink_position(:,:)       ! sink position
   real(kind=dp),intent(out),allocatable :: sink_luminosity(:)       ! luminosity
   real(kind=dp),intent(out),allocatable :: sink_temperature(:)      ! temperature
   real(kind=dp),intent(out),allocatable :: gas_position(:,:)
   real(kind=dp),intent(out),allocatable :: gas_mass(:)
   real(kind=dp),intent(out),allocatable :: gas_temperature(:)
   real(kind=dp),intent(out) :: time                            ! time
   
   
   integer :: stot                              ! number of sinks
   integer :: ptot                              ! number of sph particles
   integer :: pgas
   
   

   character(len=20) :: unit_data(1:50)        ! Char ids of arrays written
   character(len=20) :: data_id(1:50)          ! Char ids of arrays written
   character(len=20) :: format_id              ! File format (for verification)
   logical, allocatable :: ldummy(:)           ! Logical dummy array
   integer :: dim_check                        ! Dimension check
   integer :: dmdt_range_aux                   ! Accretion history array size
   integer :: i                                ! Aux. loop counter
   integer :: idata(1:50)                      ! Integer data
   integer :: j                                ! Aux. particle id
   integer :: jtot                             ! No. of particles in array
   integer :: ndata                            ! Number of arrays written
   integer :: nunit                            ! Number of units
   integer :: nvartype(1:6)                    ! ..
   integer :: p                                ! Particle counter
   integer :: pr_check                         ! Precision check
   integer :: pfirst                           ! ..
   integer :: plast                            ! ..
   integer :: typedata(1:5,1:50)               ! type data header array
   integer, allocatable :: idummy(:)           ! ..
   integer(kind=ILP) :: ilpdata(1:50)          ! Long integer data
   integer, allocatable :: ilpdummy(:)         ! ..
   real(kind=DP) :: dpdata(1:50)               ! Double precision data array
   real(kind=PR) :: rdata(1:50)                ! Real data array
   real(kind=PR), allocatable :: rdummy1(:)    ! real dummy array
   real(kind=PR), allocatable :: rdummy2(:,:)  ! real vector dummy array
   real(kind=DP), allocatable :: dpdummy1(:)   ! real dummy array
   real(kind=DP), allocatable :: dpdummy2(:,:) ! real vector dummy array
   integer :: sink_data_length                 ! ..
   integer :: idummy2(1:2)                     ! Integer dummy array 2
   integer :: idummy3(1:6)                     ! ..
   !integer :: nl,ni,nli,npr,ndp,nchar          ! ..
   integer :: s                                ! Sink counter
   real(kind=PR), allocatable :: raux(:)       ! Aux. variable

   open(1, file=out_file, status="old", form="unformatted")
!    write(6,*) "Snapshot file : ",trim(out_file),"   (SEREN binary snapshot)"

   ! Read information identifying format and precision of file
   ! Then check if each value corresponds to current Seren values
   read(1) format_id
   format_id = trim(adjustl(format_id))
!    write(6,*) "Snapshot file :",trim(out_file),"  ",trim(format_id)
   if (trim(adjustl(format_id)) /= "SERENBINARYDUMPV2") then
      stop 'Incorrect format of IC file'
   end if 
   read(1) pr_check
   if (PR == SP .and. pr_check /= 1 .and. pr_check /=4) then 
      stop 'Incorrect PR of IC file'
   else if (PR == DP .and. pr_check /= 2 .and. pr_check /=8) then 
      stop 'Incorrect PR of IC file'
   end if
   read(1) dim_check
   if (dim_check /= NDIM) then
      stop 'Incorrect NDIM of IC file'
   end if
   read(1) dim_check
   if (dim_check /= VDIM) then
      stop 'Incorrect VDIM of IC file'
   end if
   read(1) dim_check
   if (dim_check /= BDIM) then
      stop 'Incorrect BDIM of IC file'
   end if

   ! Read in all important header information
   read(1) idata
   read(1) ilpdata
   read(1) rdata
   read(1) dpdata
   ptot           = idata(1)
   stot           = idata(2)
!    pboundary      = idata(3)
!    picm           = idata(4)
     pgas           = idata(5)
!    pcdm           = idata(6)
!    pdust          = idata(7)
!    pion           = idata(8)
   nunit          = idata(20)
   ndata          = idata(21)
   dmdt_range_aux = idata(30)
!    pgas_orig      = idata(31)


!    snapshot       = ilpdata(1)
!    nsteps         = ilpdata(2)
!    ntempnext      = ilpdata(3)
!    ndiagnext      = ilpdata(4)
!    nsnapnext      = ilpdata(5)
!    nsinknext      = ilpdata(6)
   time           = dpdata(1)
!    lastsnap       = dpdata(2)
!    mgas_orig      = dpdata(3)

!    write(6,*) "SPH Particles  = ", ptot,"    Sink Particles = ", stot
!    write(6,*) "Gas            = ", pgas
!    write(6,*) "Boundary       = ", pboundary
!    write(6,*) "Intercloud     = ", picm
!    write(6,*) "Dark matter    = ", pcdm
!    write(6,*) "Dust           = ", pdust
!    write(6,*) "Ions           = ", pion
   
   
   
   ! allocate arrays
   
   allocate(sink_position(NDIM,stot))
   allocate(sink_temperature(stot))
   allocate(sink_luminosity(stot))

   allocate(gas_position(NDIM,pgas))
   allocate(gas_mass(pgas))
   allocate(gas_temperature(pgas))
   
   
   
!    if (ptot /= (pgas + picm + pboundary + pcdm + pdust + pion)) &
!          &stop "Fatal error: particles do not add up"

   if (nunit > 0) read(1) unit_data(1:nunit)
   if (ndata > 0) read(1) data_id(1:ndata)
   if (ndata > 0) read(1) typedata(1:5,1:ndata)

   sink_data_length = 19+NDIM+VDIM+2*dmdt_range_aux
   allocate(raux(1:sink_data_length))

   ! Now loop through array ids and read each array in turn
   ! ============================================================================
   do i=1,ndata

      ! Find pfirst, plast from typedata for this data set
      pfirst = typedata(2,i); plast = typedata(3,i)
      jtot = plast - pfirst + 1

      ! porig
      ! -----------------------------------------------------------------------
      if (data_id(i)=="porig") then
         allocate(idummy(1:jtot))
         read(1) idummy
         deallocate(idummy)

      ! Positions
      ! -----------------------------------------------------------------------
      else if (data_id(i)=="r") then
         allocate(rdummy2(1:NDIM,1:jtot))
         read(1) rdummy2
         gas_position=rdummy2
         deallocate(rdummy2)

      ! Masses
      ! -----------------------------------------------------------------------
      else if (data_id(i)=="m") then
         allocate(rdummy1(1:jtot))
         read(1) rdummy1
         gas_mass=rdummy1
         deallocate(rdummy1)

      ! Smoothing length
      ! -----------------------------------------------------------------------
      else if (data_id(i)=="h") then
         allocate(rdummy1(1:jtot))
         read(1) rdummy1
         deallocate(rdummy1)

      ! Velocities
      ! -----------------------------------------------------------------------
      else if (data_id(i)=="v") then
         allocate(rdummy2(1:VDIM,1:jtot))
         read(1) rdummy2
         deallocate(rdummy2)

      ! Density
      ! -----------------------------------------------------------------------
      else if (data_id(i)=="rho") then
         allocate(rdummy1(1:jtot))
         read(1) rdummy1
         deallocate(rdummy1)

      ! Temperature
      ! -----------------------------------------------------------------------
      else if (data_id(i)=='temp') then
         allocate(rdummy1(1:jtot))
         read(1) rdummy1
         gas_temperature=rdummy1
         deallocate(rdummy1)

      ! Internal energy
      ! -----------------------------------------------------------------------
      else if (data_id(i)=='u') then
         allocate(rdummy1(1:jtot))
         read(1) rdummy1
         deallocate(rdummy1)

      ! B-field
      ! -----------------------------------------------------------------------
      else if (data_id(i)=='B') then
         allocate(rdummy2(1:BDIM,pfirst:plast))
         read(1) rdummy2
         deallocate(rdummy2)

      ! Sinks
      ! -----------------------------------------------------------------------
      else if (data_id(i)=='sink_v1') then
         !read(1) nl,ni,nli,npr,ndp,nchar
         read(1) idummy3(1:6)
         allocate(ldummy(1:3))
         if (stot > 0) then
            do j=pfirst,plast
               s = j
               read(1) ldummy
               read(1) idummy2
               read(1) raux
               
               
               sink_position(:,s)=raux(2:NDIM+1)
               sink_luminosity(s)=raux(NDIM+VDIM+18+2*DMDT_RANGE)
               sink_temperature(s)=raux(NDIM+VDIM+19+2*DMDT_RANGE)
               
               
               
               
!                sink(s)%id          = idummy2(1)
!                sink(s)%ncreate     = idummy2(2)
!                sink(s)%accrete     = ldummy(1)
!                sink(s)%static      = ldummy(2)
!                sink(s)%tcreate     = real(raux(1),DP)
!                sink(s)%r(1:NDIM)   = raux(2:NDIM+1)
!                sink(s)%v(1:NDIM)   = raux(NDIM+2:NDIM+VDIM+1)
!                sink(s)%m           = raux(NDIM+VDIM+2)
!                sink(s)%h           = raux(NDIM+VDIM+3)
!                sink(s)%radius      = raux(NDIM+VDIM+4)
!                sink(s)%angmom(1:3) = raux(NDIM+VDIM+5:NDIM+VDIM+7)
!                sink(s)%dmdt        = raux(NDIM+VDIM+8)
!                sink(s)%star_radius = raux(NDIM+VDIM+9)
!                sink(s)%luminosity  = raux(NDIM+VDIM+10)
!                sink(s)%temperature = raux(NDIM+VDIM+11)
!                if (dmdt_range_aux > 0) then
!                   dmdt_range_min = min(dmdt_range_aux,DMDT_RANGE)
!                   sink(s)%macc = 0.0_PR
!                   sink(s)%macc(1:dmdt_range_min) = &
!                         &real(raux(NDIM+VDIM+12:&
!                         &NDIM+VDIM+11+dmdt_range_min),DP)
!                   sink(s)%tacc = 0.0_PR
!                   sink(s)%tacc(1:dmdt_range_min) = &
!                         &real(raux(NDIM+VDIM+12+dmdt_range_aux:&
!                         &NDIM+VDIM+11+dmdt_range_min+dmdt_range_aux),DP)
!                end if
!                sink(s)%mmax        = raux(NDIM+VDIM+12+2*dmdt_range_aux)
!                sink(s)%m_star = raux(NDIM+VDIM+13+2*DMDT_RANGE)
!                sink(s)%m_star_acc = raux(NDIM+VDIM+14+2*DMDT_RANGE)
!                sink(s)%t0_acc = raux(NDIM+VDIM+15+2*DMDT_RANGE)
!                sink(s)%delta_t = raux(NDIM+VDIM+16+2*DMDT_RANGE)
!                sink(s)%prev_time = raux(NDIM+VDIM+17+2*DMDT_RANGE)
!                sink(s)%star_lum = raux(NDIM+VDIM+18+2*DMDT_RANGE)
!                sink(s)%star_temp = raux(NDIM+VDIM+19+2*DMDT_RANGE)
!                sink(s)%rapid_acc = ldummy(3)
            end do
         end if
         deallocate(ldummy)

      ! Skip through arbitrary 1-D or 2-D array
      ! -----------------------------------------------------------------------
      else if (typedata(1,i) >= 1) then
         if (typedata(4,i)==1) then
            allocate(ldummy(1:plast-pfirst+1))
            read(1) ldummy
            deallocate(ldummy)
         else if (typedata(4,i)==2) then
            allocate(idummy(1:plast-pfirst+1))
            read(1) idummy
            deallocate(idummy)
         else if (typedata(4,i)==3) then
            allocate(ilpdummy(1:plast-pfirst+1))
            read(1) ilpdummy
            deallocate(ilpdummy)
         else if (typedata(4,i)==4 .and. typedata(1,i)==1) then
            allocate(rdummy1(1:plast-pfirst+1))
            read(1) rdummy1
            deallocate(rdummy1)
         else if (typedata(4,i)==4 .and. typedata(1,i)>1) then
            allocate(rdummy2(1:typedata(1,i),1:plast-pfirst+1))
            read(1) rdummy2
            deallocate(rdummy2)
         else if (typedata(4,i)==5 .and. typedata(1,i)==1) then
            allocate(dpdummy1(1:plast-pfirst+1))
            read(1) dpdummy1
            deallocate(dpdummy1)
         else if (typedata(4,i)==5 .and. typedata(1,i)>1) then
            allocate(dpdummy2(1:typedata(1,i),1:plast-pfirst+1))
            read(1) dpdummy2
            deallocate(dpdummy2)
         end if

      ! Skip through arbitrary data structure
      ! -----------------------------------------------------------------------
      else if (typedata(1,i) == 0 .and. typedata(4,i) == 7) then
         read(1) nvartype
         if (nvartype(1) > 0) allocate(ldummy(1:nvartype(1)))
         if (nvartype(2) > 0) allocate(idummy(1:nvartype(2)))
         if (nvartype(3) > 0) allocate(ilpdummy(1:nvartype(3)))
         if (nvartype(4) > 0) allocate(rdummy1(1:nvartype(4)))
         if (nvartype(5) > 0) allocate(dpdummy1(1:nvartype(5)))
         do p=pfirst,plast
            if (nvartype(1) > 0) read(1) ldummy
            if (nvartype(2) > 0) read(1) idummy
            if (nvartype(3) > 0) read(1) ilpdummy
            if (nvartype(4) > 0) read(1) rdummy1
            if (nvartype(5) > 0) read(1) dpdummy1
         end do
         if (nvartype(5) > 0) deallocate(dpdummy1)
         if (nvartype(4) > 0) deallocate(rdummy1)
         if (nvartype(3) > 0) deallocate(ilpdummy)
         if (nvartype(2) > 0) deallocate(idummy)
         if (nvartype(1) > 0) deallocate(ldummy)

      end if
   end do
   ! ============================================================================

   ! Close file once finished
   ! ----------------------------------------------------------------------------
   close(1)

   return
   END SUBROUTINE read_data_seren_unform



end module


program sink_conv

   use sink_mod

   implicit none
   
   ! declare variables
   character(len=1024) :: run_id
   character(len=1024) :: seren_file_name
   character(len=1024) :: snap_file_name
   character(len=1024) :: sink_file_name
   character(len=1024) :: params_file_name
   
   integer :: i,j,k,l
   integer :: n_vel,n_sim,n_snap
   integer,allocatable :: id(:)  
   
   logical :: file_exists
   
   real(kind=dp),allocatable :: position(:,:)
   real(kind=dp),allocatable :: mass(:)
   real(kind=dp),allocatable :: temperature(:)
   real(kind=dp),allocatable :: star_position(:,:)
   real(kind=dp),allocatable :: star_luminosity(:)
   real(kind=dp),allocatable :: star_temperature(:)
   real(kind=dp) :: time
   
   logical :: hot
   
   
   
   do j=1,400
   
   
   write(seren_file_name,"(A,I5.5)") "/Users/oliverlomax/Documents/turb_study_ii/m3_r3000_sig0.44_0_2_002.su.",j
   
   write(*,*) trim(seren_file_name)
   
      call read_data_seren_unform&
         &(seren_file_name,star_position,star_luminosity,star_temperature,position,mass,temperature,time)
      
      write(snap_file_name,"(A,I3.3,A)") "face_on_r2_",j,"/iteration_02.bin"
       
!       open(1,file=trim(snap_file_name))
!       write(1,"(5(A25))") "x_1", "x_2", "x_3", "m", "T"
!       write(1,"(5(A25))") "pc", "pc", "pc", "M_sun", "K"
!       do i=1,size(position,2)
!          write(1,"(5(E25.17))") position(:,i),mass(i),temperature(i)
!       end do
!       close(1)
      
      write(sink_file_name,"(A,I3.3,A)") "stars_",j,".dat"
      
!       open(1,file=trim(sink_file_name))
!       write(1,"(5(A25))") "x_1", "x_2", "x_3", "L", "T"
!       write(1,"(5(A25))") "pc", "pc", "pc", "W", "K"
!       do i=1,size(star_position,2)
!          write(1,"(5(E25.17))") star_position(:,i),star_luminosity(i),star_temperature(i)
!       end do
!       close(1)
      hot=any(star_luminosity>1.e27)
      if (hot) write(*,*) "HOT"
      
      write(params_file_name,"(A,I3.3,A)") "params_edge_on_r2_",j,".dat"
      open(1,file=trim(params_file_name))
      
      write(run_id,"(A,I3.3)") "edge_on_r2_",j
      
      write(1,"(A)") "sim_cloud_file = "//trim(snap_file_name)
      if (hot) then
         write(1,"(A)") "sim_n_packet_point=3000000"
      else
         write(1,"(A)") "sim_n_packet_point=3000000"
      end if
      write(1,"(A)") "sim_n_packet_external=3000000"
      write(1,"(A)") "sim_star_file = "//trim(sink_file_name)
      write(1,"(A)") "sim_id = "//trim(run_id)
      write(1,"(A)") "sim_n_it= 0"
      write(1,"(A)") "sim_restart = .false."
      write(1,"(A)") "sim_min_d = 1.6256865e+12"
      write(1,"(A)") "sph_scattered_light=.true."
      write(1,"(A)") "point_sources=.true."
      write(1,"(A)") "datacube_angles = 1.5707963 0. 0."
      write(1,"(A)") "datacube_convolve=.true."
      write(1,"(A)") "datacube_n_x = 512"
      write(1,"(A)") "datacube_n_y = 512"
      write(1,"(A)") "datacube_x_min = -1.234271e+16"
      write(1,"(A)") "datacube_x_max = 1.234271e+16"
      write(1,"(A)") "datacube_y_min = -1.234271e+16"
      write(1,"(A)") "datacube_y_max = 1.234271e+16"
      write(1,"(A)") "datacube_make=.true."
      write(1,"(A)") "datacube_lambda_string = 1.000E-01 1.122E-01 1.259E-01 1.413E-01 1.585E-01 1.778E-01 1.995E-01 2.239E-01&
         & 2.512E-01 2.818E-01 3.162E-01 3.548E-01 3.981E-01 4.467E-01 5.012E-01 5.623E-01 6.310E-01 7.079E-01 7.943E-01 8.913E-01&
         & 1.000E+00 1.122E+00 1.259E+00 1.413E+00 1.585E+00 1.778E+00 1.995E+00 2.239E+00 2.512E+00 2.818E+00 3.162E+00 3.548E+00&
         & 3.981E+00 4.467E+00 5.012E+00 5.623E+00 6.310E+00 7.079E+00 7.943E+00 8.913E+00 1.000E+01 1.122E+01 1.259E+01 1.413E+01&
         & 1.585E+01 1.778E+01 1.995E+01 2.239E+01 2.512E+01 2.818E+01 3.162E+01 3.548E+01 3.981E+01 4.467E+01 5.012E+01 5.623E+01&
         & 6.310E+01 7.079E+01 7.943E+01 8.913E+01 1.000E+02 1.122E+02 1.259E+02 1.413E+02 1.585E+02 1.778E+02 1.995E+02 2.239E+02&
         & 2.512E+02 2.818E+02 3.162E+02 3.548E+02 3.981E+02 4.467E+02 5.012E+02 5.623E+02 6.310E+02 7.079E+02 7.943E+02 8.913E+02&
         & 1.000E+03 1.122E+03 1.259E+03 1.413E+03 1.585E+03 1.778E+03 1.995E+03 2.239E+03 2.512E+03 2.818E+03 3.162E+03 3.548E+03&
         & 3.981E+03 4.467E+03 5.012E+03 5.623E+03 6.310E+03 7.079E+03 7.943E+03 8.913E+03 1.000E+04"
      write(1,"(A)") "datacube_fwhm_string =   1.678E-04 1.882E-04 2.112E-04 2.370E-04 2.659E-04 2.983E-04 3.347E-04 3.756E-04&
         & 4.214E-04 4.728E-04 5.305E-04 5.952E-04 6.679E-04 7.494E-04 8.408E-04 9.434E-04 1.059E-03 1.188E-03 1.333E-03 1.495E-03&
         & 1.678E-03 1.882E-03 2.112E-03 2.370E-03 2.659E-03 2.983E-03 3.347E-03 3.756E-03 4.214E-03 4.728E-03 5.305E-03 5.952E-03&
         & 6.679E-03 7.494E-03 8.408E-03 9.434E-03 1.059E-02 1.188E-02 1.333E-02 1.495E-02 1.678E-02 1.882E-02 2.112E-02 2.370E-02&
         & 2.659E-02 2.983E-02 3.347E-02 3.756E-02 4.214E-02 4.728E-02 5.305E-02 5.952E-02 6.679E-02 7.494E-02 8.408E-02 9.434E-02&
         & 1.059E-01 1.188E-01 1.333E-01 1.495E-01 1.678E-01 1.882E-01 2.112E-01 2.370E-01 2.659E-01 2.983E-01 3.347E-01 3.756E-01&
         & 4.214E-01 4.728E-01 5.305E-01 5.952E-01 6.679E-01 7.494E-01 8.408E-01 9.434E-01 1.059E+00 1.188E+00 1.333E+00 1.495E+00&
         & 1.678E+00 1.882E+00 2.112E+00 2.370E+00 2.659E+00 2.983E+00 3.347E+00 3.756E+00 4.214E+00 4.728E+00 5.305E+00 5.952E+00&
         & 6.679E+00 7.494E+00 8.408E+00 9.434E+00 1.059E+01 1.188E+01 1.333E+01 1.495E+01 1.678E+01"
      write(1,"(A)") "datacube_distance = 3.0856776e+20"
      
      close(1)
      
      
      
   
   end do
   
   
   stop
   
   

end program