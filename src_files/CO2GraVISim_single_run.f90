program CO2GraVISim_single_run
   ! This program performs a single run of CO2GraVISim.
   ! Various input parameters and arrays are read into memory, the main solver program
   ! is called, and then the corresponding output files are written.

   use CO2GraVISim_global !wp, iter_max, tmax, dt_init, dt_min, errmax_t_step, errmax_P_iter, P_iter_omega, t0
   use CO2GraVISim_input_parameters  !H0, B0, D0, perm_h, poro_h
   use CO2GraVISim_injection_profiles !n_inj_locs, n_flux_times, Q_inj_locs, Q_flux_vals, t_flux_vals
   use CO2GraVISim_solver

   implicit none
   integer :: j_x, j_y, j_p, io
   real(wp), dimension(:,:,:), allocatable :: h_array, h_res_array, P_array
   real(wp), dimension(:),     allocatable :: target_plot_times, plot_times

   !Input files read in directly here
   character(50)  ::  target_plot_times_file = "./Input/target_plot_times.txt"

   !Output files
   character(50)  ::  h_output_file          = "./Output/Current_Thickness/h_0.txt"
   character(50)  ::  h_res_output_file      = "./Output/Current_Thickness/h_res_0.txt"
   character(50)  ::  P_output_file          = "./Output/Current_Pressure/P_0.txt"
   character(50)  ::  output_times_file      = "./Output/Other/plot_times.txt"
   character(50)  ::  param_output_file      = "./Output/Other/parameter_values.txt"
   character(50)  ::  inj_locs_output_file   = "./Output/Other/injection_locations.txt"



   ! -- Read in Inputs ----------------------------------------------------------------------------------------------

   call read_flow_parameters  !from the -confined_current_input_parameters- module
   call read_reservoir_inputs !from the -confined_current_input_parameters- module
   call read_injection_data   !from the -Injection_profiles- module


   write(*,*) 'Reading in plot times'
   !Read in plot times
   open (newunit=io, file=target_plot_times_file, status='old', action='read')
   read(io,*)     !Skip the first line, as it's formatting instructions
   read(io, *) np !Read in the number of plot times
   allocate (target_plot_times(0:np-1))
   allocate (plot_times(0:np-1))
   do j_x = 0, np-1
      read(io, *) target_plot_times(j_x)
   end do
   close (io)

   !Set some temporary values in the array of actual plot_times
   !These will be close to the target_plot_times values, but 
   !may differ slightly due to the time-stepping routine used
   plot_times = target_plot_times

   write(*,*) 'Allocating arrays'
   !Allocated array sizes
   allocate(h_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(h_res_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(P_array(0:nx-1, 0:ny-1, 0:np-1))

   ! -- Run Solver ---------------------------------------------------------------------------------------------------

   write(*,*) '-- Running solver --'

   call solver_run(h_array,h_res_array,P_array,plot_times,target_plot_times)

   write(*,*) '-- Solver run complete --'


   ! -- Save Outputs -------------------------------------------------------------------------------------------------

   write(*,*) 'Saving outputs'


   !Save arrays that have contain a different slice for each plot time
   write(*,*) ' - saving mobile current profiles'
   do j_p = 0, np-1

      write(h_output_file, '("./Output/Current_Thickness/h", I2.2, ".txt")') j_p

      open (newunit=io, file=h_output_file, action='write')
      do j_y = 0, ny-1
         write(io, *) h_array(:,j_y,j_p)
      end do
      close(io)

   end do

   write(*,*) ' - saving trapped current profiles'
   do j_p = 0, np-1

      write(h_res_output_file, '("./Output/Current_Thickness/h_res", I2.2, ".txt")') j_p

      open (newunit=io, file=h_res_output_file, action='write')
      do j_y = 0, ny-1
         write(io, *) h_res_array(:,j_y,j_p)
      end do
      close(io)

   end do

   write(*,*) ' - saving ambient pressure profiles'
   do j_p = 0, np-1

      write(P_output_file, '("./Output/Current_Pressure/P", I2.2, ".txt")') j_p

      open (newunit=io, file=P_output_file, action='write')
      do j_y = 0, ny-1
         write(io, *) P_array(:,j_y,j_p)
      end do
      close(io)

   end do


   write(*,*) ' - saving plot times'
   open (newunit=io, file=output_times_file, action='write')
   do j_p = 0, np-1
      write(io, *) plot_times(j_p)
   end do
   close(io)


   write(*,*) ' - saving run parameters'
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
   open (newunit=io, file=inj_locs_output_file, action='write')
   do j_p = 0, n_inj_locs-1
      write(io, *) Q_inj_locs(j_p,0), Q_inj_locs(j_p,1)
   end do
   close(io)


   write(*,*) 'Run complete.'

end program CO2GraVISim_single_run
