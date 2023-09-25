module CO2GraVISim_global
   !This module contains parameter values used by the main solver that are not 
   !expected to be changed between individual runs.

   implicit none
   integer, parameter  :: wp=selected_real_kind(8)

   !Solver parameters
   integer, parameter  :: itermax         = 1e6       ! max number of iterations
   real(wp), parameter :: tmax            = 1e4_wp    ! max run time
   real(wp), parameter :: dt_init         = 1e-5_wp   ! initial timestep
   real(wp), parameter :: dt_min          = 1e-8_wp   ! min timestep
   real(wp), parameter :: errmax_t_step   = 1e-5_wp   ! max error allowed in timestep check
   real(wp), parameter :: errmax_P_iter   = 1e-4_wp   ! max error allowed in Pressure solver
   real(wp), parameter :: P_iter_omega    = 1.5_wp    ! overrelaxation parameter for Pressure solver
   real(wp), parameter :: t0              = 0._wp     ! start time

end module CO2GraVISim_global
