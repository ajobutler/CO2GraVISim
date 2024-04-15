program CO2GraVISim_single_run
   ! This program performs a single run of CO2GraVISim.
   ! Various input parameters and arrays are read into memory, the main solver program
   ! is called, and then the corresponding output files are written.

   use CO2GraVISim_global !wp, iter_max, tmax, dt_init, dt_min, errmax_t_step, errmax_P_iter, P_iter_omega, t0
   use CO2GraVISim_input_parameters  !H0, B0, D0, perm_h, poro_h, n_inj_locs, n_flux_times, Q_inj_locs, Q_flux_vals, t_flux_vals
   use CO2GraVISim_solver

   implicit none
   integer :: j_y, j_p, io, i_arg, arg_count
   real(wp), dimension(:,:,:), allocatable :: h_array, h_res_array, P_array
   real(wp), dimension(:),     allocatable :: plot_times

   !Main Input and Output folders
   character(1000)  :: input_folder       !Default is "./Input"
   character(1000)  :: output_folder      !Default is "./Output"
   character(1000)  :: arg

   !Output files
   character(1000)  ::  h_output_file          !Default is "./Output/Current_Thickness/h_n.txt"       for n=0,1,...
   character(1000)  ::  h_res_output_file      !Default is "./Output/Current_Thickness/h_res_n.txt"   for n=0,1,...
   character(1000)  ::  P_output_file          !Default is "./Output/Current_Pressure/P_n.txt"        for n=0,1,...   
   character(1000)  ::  output_times_file      !Default is "./Output/Other/plot_times.txt"
   character(1000)  ::  param_output_file      !Default is "./Output/Other/parameter_values.txt"
   character(1000)  ::  inj_locs_output_file   !Default is "./Output/Other/injection_locations.txt"


   ! -- Read in Command Line inputs --------------------------------------------------------------------------------

   !Default values
   input_folder  = './Input'
   output_folder = './Output'


   !Parse Command Line inputs - based on code by Jon D'Souza-Eva
   arg_count = command_argument_count()
   i_arg = 1
   do while( i_arg < arg_count )
      call get_command_argument(i_arg, arg)

      select case(arg)

      case('-input')
         if (i_arg < arg_count) then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, input_folder)
         endif

      case('-output')
         if (i_arg < arg_count) then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, output_folder)
         endif

      case default
         write(*,*) 'Unknown argument: ', trim(arg)
         stop

      end select
      
      i_arg = i_arg + 1

   end do
   
   !Remove any trailing '/' or '\' that might have been included
   call trim_slashes(input_folder)
   call trim_slashes(output_folder)

   write(*,*) 'Input Folder : [' , trim(input_folder)  , ']'
   write(*,*) 'Output Folder: [' , trim(output_folder) , ']'
   

   ! -- Read in Inputs ----------------------------------------------------------------------------------------------

   call read_reservoir_inputs(input_folder) !from the -confined_current_input_parameters- module

   !Set some temporary values in the array of actual plot_times.
   !The final values will be close to the target_plot_times values, but 
   !may differ slightly due to the time-stepping routine used
   allocate (plot_times(0:np-1))
   plot_times = target_plot_times

   write(*,*) 'Allocating arrays'
   !Allocated array sizes
   allocate(h_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(h_res_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(P_array(0:nx-1, 0:ny-1, 0:np-1))

   ! -- Run Solver ---------------------------------------------------------------------------------------------------

   write(*,*) '-- Running solver --'

   call solver_run(h_array,h_res_array,P_array,plot_times)

   write(*,*) '-- Solver run complete --'


   ! -- Save Outputs -------------------------------------------------------------------------------------------------

   write(*,*) 'Saving outputs'

   !Save arrays that contain a different slice for each plot time
   write(*,*) ' - saving mobile current profiles'
   do j_p = 0, np-1

      write(h_output_file, '( A , "/Current_Thickness/h", I2.2, ".txt")') trim(output_folder), j_p

      open (newunit=io, file=h_output_file, action='write')
      do j_y = 0, ny-1
         write(io, *) h_array(:,j_y,j_p)
      end do
      close(io)

   end do


   write(*,*) ' - saving trapped current profiles'
   do j_p = 0, np-1

      write(h_res_output_file, '( A , "/Current_Thickness/h_res", I2.2, ".txt")') trim(output_folder), j_p

      open (newunit=io, file=h_res_output_file, action='write')
      do j_y = 0, ny-1
         write(io, *) h_res_array(:,j_y,j_p)
      end do
      close(io)

   end do


   write(*,*) ' - saving ambient pressure profiles'
   do j_p = 0, np-1

      write(P_output_file, '( A , "/Current_Pressure/P", I2.2, ".txt")') trim(output_folder), j_p

      open (newunit=io, file=P_output_file, action='write')
      do j_y = 0, ny-1
         write(io, *) P_array(:,j_y,j_p)
      end do
      close(io)

   end do


   write(*,*) ' - saving plot times'
   write(output_times_file, '( A , "/Other/plot_times.txt")') trim(output_folder)
   open (newunit=io, file=output_times_file, action='write')
   do j_p = 0, np-1
      write(io, *) plot_times(j_p)
   end do
   close(io)

   
   write(*,*) ' - saving run parameters'
   write(param_output_file, '( A , "/Other/parameter_values.txt")') trim(output_folder)
   open (newunit=io, file=param_output_file, action='write')
   write(io, '(I5)') nx
   write(io, '(I5)') ny
   write(io, '(F10.3)') dx
   write(io, '(F10.3)') dy
   write(io, '(F10.3)') M
   write(io, '(F10.3)') Gamma_val
   write(io, '(F10.3)') s_c_r
   write(io, '(F10.3)') s_a_i
   write(io, '(F10.3)') C_sat
   write(io, '(F10.3)') q_dissolve
   close (io)


   write(*,*) ' - saving injection locations'
   write(inj_locs_output_file, '( A , "/Other/injection_locations.txt")') trim(output_folder)
   open (newunit=io, file=inj_locs_output_file, action='write')
   do j_p = 0, n_inj_locs-1
      write(io, *) Q_inj_locs(j_p,0), Q_inj_locs(j_p,1)
   end do
   close(io)


   write(*,*) 'Run complete.'

end program CO2GraVISim_single_run
