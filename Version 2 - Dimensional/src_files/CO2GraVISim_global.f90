module CO2GraVISim_global
   !This module contains parameter values used by the main solver that are not 
   !expected to be changed between individual runs.

   implicit none
   integer, parameter  :: wp=selected_real_kind(8)

   !Solver parameters
   ! integer, parameter  :: itermax         = 1e9       ! max number of iterations
   ! real(wp), parameter :: tmax            = 1e9_wp    ! max run time
   ! real(wp), parameter :: dt_init         = 1e-5_wp   ! initial timestep
   ! real(wp), parameter :: dt_min          = 1e-8_wp   ! min timestep
   ! real(wp), parameter :: errmax_t_step   = 1e-4_wp   ! max error allowed in timestep check
   ! real(wp), parameter :: errmax_P_iter   = 1e-4_wp   ! max error allowed in Pressure solver
   ! real(wp), parameter :: P_iter_omega    = 1.96_wp   ! overrelaxation parameter for Pressure solver
   ! real(wp), parameter :: t0              = 0._wp     ! start time
   ! integer, parameter  :: P_iter_check    = 10        ! How often to check pressure convergence (# of iterations)

   !Solver parameters
   ! - Fixed values
   real(wp), parameter :: t0  = 0._wp                          ! start time
   ! - Default values that can be altered via the XML input
   integer  :: itermax        = 1e9                            ! max number of iterations
   real(wp) :: tmax           = 1e9_wp                         ! max run time
   real(wp) :: dt_init        = 1e-5_wp                        ! initial timestep
   real(wp) :: dt_min         = 1e-8_wp                        ! min timestep
   real(wp) :: errmax_t_step  = 1e-5_wp                        ! max error allowed in timestep check
   real(wp) :: errmax_P_iter  = 5e-4_wp                        ! max error allowed in Pressure solver
   real(wp) :: P_iter_omega   = 1.96_wp                        ! overrelaxation parameter for Pressure solver
   integer  :: P_iter_check   = 10                             ! How often to check pressure convergence (# of iterations)



   !Dimensional conversion factors
   real(wp), parameter :: seconds_to_days = 1._wp / (24._wp*60._wp*60._wp) !Time conversion: seconds to days
   real(wp), parameter :: days_to_seconds = (24._wp*60._wp*60._wp)         !Time conversion: days to seconds
   real(wp), parameter :: mD_to_m2 = 9.869e-16_wp                          !Permeability conversion: milliDarcy to m^2

end module CO2GraVISim_global
