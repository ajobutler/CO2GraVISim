module CO2GraVISim_timer
   !This module is used to time the various subroutines involved in the code,
   !in order to determine where the solver is spending most of its time

   use CO2GraVISim_global

   implicit none
   !runtimes
   real(wp) :: ADI_X_time = 0._wp
   real(wp) :: ADI_Y_time = 0._wp
   real(wp) :: h_coeff_time = 0._wp
   real(wp) :: h_res_time = 0._wp
   real(wp) :: Ilin_Coeff_time = 0._wp
   real(wp) :: gtri_time = 0._wp
   real(wp) :: err_time = 0._wp
   real(wp) :: save_iter_time = 0._wp
   real(wp) :: Gen_flux_time = 0._wp
   real(wp) :: Res_amb_time = 0._wp
   real(wp) :: Array_grad_time = 0._wp
   real(wp) :: Array_Lap_time = 0._wp
   real(wp) :: vert_interp_time = 0._wp
   real(wp) :: vert_interp_time_a = 0._wp
   real(wp) :: vert_interp_time_b = 0._wp
   real(wp) :: vert_interp_time_c = 0._wp
   real(wp) :: vert_interp_time_d = 0._wp
   real(wp) :: vert_interp_time_e = 0._wp
   real(wp) :: P_solve_time = 0._wp
   real(wp) :: P_SOR_time = 0._wp
   real(wp) :: P_BCs_time = 0._wp
   !call numbers
   integer :: ADI_X_calls = 0
   integer :: ADI_Y_calls = 0
   integer :: h_coeff_calls = 0
   integer :: h_res_calls = 0
   integer :: Ilin_Coeff_calls = 0
   integer :: gtri_calls = 0
   integer :: err_calls = 0
   integer :: save_iter_calls = 0
   integer :: Gen_flux_calls = 0
   integer :: Res_amb_calls = 0
   integer :: Array_grad_calls = 0
   integer :: Array_Lap_calls = 0
   integer :: vert_interp_calls = 0
   integer :: P_solve_calls = 0
   integer :: P_SOR_calls = 0
   integer :: P_BCs_calls = 0

   integer :: fid_runtimes

   Private :: fid_runtimes

contains

   subroutine routine_runtimes(total_runtime,Output_folder)

      implicit none
      real(wp), intent(in) :: total_runtime
      character(len=*), intent(in) :: Output_folder
      character(1000) :: runtimes_file

      runtimes_file = trim(Output_folder) // "/Other/runtimes.txt"

      open (newunit=fid_runtimes, file=runtimes_file, action='write')

      write(fid_runtimes,*) ''
      write(fid_runtimes,*) ''
      write(fid_runtimes,'(A,F16.3,A)') 'Total runtime = ', total_runtime, ' seconds.'
      write(fid_runtimes,*) ''
      write(fid_runtimes,*) ''
      write(fid_runtimes,'(A)') 'Routine runtimes:'
      write(fid_runtimes,*) ''
      write(fid_runtimes,'(A13,A3,A17,A3,A18)') ' No. of calls', ' | ', 'runtime (seconds)', ' | ', '% of total runtime'

      call timing_string('ADI_X_solve', ADI_X_calls, ADI_X_time, total_runtime)
      call timing_string('ADI_Y_solve', ADI_Y_calls, ADI_Y_time, total_runtime)
      call timing_string('h_coefficients', h_coeff_calls, h_coeff_time, total_runtime)
      call timing_string('h_res_update', h_res_calls, h_res_time, total_runtime)
      call timing_string('Ilin_Coefficients', Ilin_Coeff_calls, Ilin_Coeff_time, total_runtime)
      call timing_string('gtri', gtri_calls, gtri_time, total_runtime)
      call timing_string('err', err_calls, err_time, total_runtime)
      call timing_string('save_iter_output', save_iter_calls, save_iter_time, total_runtime)
      call timing_string('generate_flux_array', Gen_flux_calls, Gen_flux_time, total_runtime)
      call timing_string('Residual_ambient_flux', Res_amb_calls, Res_amb_time, total_runtime)
      call timing_string('Array_gradient', Array_grad_calls, Array_grad_time, total_runtime)
      call timing_string('Array_Laplacian', Array_Lap_calls, Array_Lap_time, total_runtime)
      call timing_string('vertical_interpolation', vert_interp_calls, vert_interp_time, total_runtime)
      !call timing_string('vertical_interpolation_a', vert_interp_calls, vert_interp_time_a, total_runtime)
      !call timing_string('vertical_interpolation_b', vert_interp_calls, vert_interp_time_b, total_runtime)
      !call timing_string('vertical_interpolation_c', vert_interp_calls, vert_interp_time_c, total_runtime)
      !call timing_string('vertical_interpolation_d', vert_interp_calls, vert_interp_time_d, total_runtime)
      !call timing_string('vertical_interpolation_e', vert_interp_calls, vert_interp_time_e, total_runtime)
      call timing_string('Pressure_solve', P_solve_calls, P_solve_time, total_runtime)
      call timing_string('Pressure_SOR', P_SOR_calls, P_SOR_time, total_runtime)
      call timing_string('Pressure_BCs', P_BCs_calls, P_BCs_time, total_runtime)

      write(*,*) ''

      close(fid_runtimes)

   end subroutine routine_runtimes

   subroutine timing_string(subroutine_name, call_total, subroutine_time, total_time)

      implicit none
      character(len=*), intent(in) :: subroutine_name
      real(wp), intent(in) :: subroutine_time, total_time
      integer, intent(in) :: call_total

      write(fid_runtimes,'(A)') repeat('-',54)
      write(fid_runtimes,'(A)') trim(subroutine_name)
      write(fid_runtimes,'(I13,A3,F17.3,A3,F6.3)') &
      & call_total, ' | ', subroutine_time , ' | ', 100._wp * (subroutine_time/total_time)
      !write(*,'(A)'), repeat('-',54)


   end subroutine timing_string


end module CO2GraVISim_timer
