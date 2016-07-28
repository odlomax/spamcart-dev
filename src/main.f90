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

! SPAMCART main program
program p_main

   use m_kind_parameters
   use m_simulation
   use m_options
   
   implicit none
   
   type(simulation) :: mcrt_sim                                ! simulation object
   type(options) :: sim_params                                 ! parameters object
   
   integer(kind=int_kind) :: i                                 ! counter
   character(kind=chr_kind,len=string_length) :: params_file   ! parameters file
   
   ! get parameters file, if present
   call get_command_argument(1,params_file)
   
   ! initialise parameters
   call sim_params%initialise(params_file)
   
   ! initialise simulation
   call mcrt_sim%initialise(sim_params)
   
   ! run simulation
   do i=1,sim_params%sim_n_it
      call mcrt_sim%perform_iteration(i)
   end do
   
   ! make datacube
   if (sim_params%datacube_make) call mcrt_sim%make_datacube()

end program