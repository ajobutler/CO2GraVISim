module CO2GraVISim_solver
   !This is the main solver for CO2GraVISim, calculating the mobile current thickness h,
   !the resdiually trapped current thickness h_res, and the ambient pressure Pa, for
   !given inputs and at the specified plot times.


   use CO2GraVISim_global !wp, iter_max, tmax, dt_init, dt_min, errmax_t_step, errmax_P_iter, P_iter_omega, t0
   use CO2GraVISim_input_parameters !H0, B0, D0, perm_h, poro_h
   !use CO2GraVISim_injection_profiles !n_inj_locs, n_flux_times, Q_inj_locs, Q_flux_vals, t_flux_vals
   use CO2GraVISim_vertical_structure !poro_cur
   use CO2GraVISim_iter_calculations
   use CO2GraVISim_logging
   use CO2GraVISim_h_routines
   use CO2GraVISim_h_res_routines
   use CO2GraVISim_pressure_routines
   use CO2GraVISim_err_calc
   use CO2GraVISim_timer
   use OMP_LIB

   implicit none
   ! integer :: fid_Volume !, fid_P_iter, fid_Iter_info !file IDs for writing to output files
   ! integer :: iter_save_old = 0
   real(wp) :: start, finish, runtime, prev_rel_err_val, dt_prev

contains

   ! ------------------------------------------------------------------------
   subroutine solver_run(h_array,h_res_array,P_array,V_mobile_array,V_trapped_array,&
   & Max_mobile_thickness,Max_pressure,Min_pressure,plot_times,Output_folder,Confined_flag)

      implicit none
      real(wp), dimension(0:nx-1,0:ny-1,0:np-1), intent(out) :: h_array, h_res_array, P_array
      real(wp), dimension(0:nx-1,0:ny-1,0:np-1), intent(out) :: V_mobile_array, V_trapped_array
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: Max_mobile_thickness, Max_pressure, Min_pressure
      real(wp), dimension(0:np-1), intent(out) :: plot_times
      character(1000), intent(in) :: Output_folder
      logical, intent(in) :: Confined_flag
      real(wp), dimension(0:nx-1,0:ny-1) :: h, h_res, h_old, dh_dt, P, P_old, Q
      real(wp) :: t, dt, dt_ideal, dt_new, tplot, V_injected, V_tot_dsln_loss, V_boundary
      integer  :: it, plot_count, flux_idx, P_iter_count
      logical :: flux_change_flag = .False.
      integer, dimension(0:1) :: maxloc_idxs

      !output file for the various calculated volumes after each successful iteration
      ! character(1000) :: Volumes_file, P_iter_file, Iter_info_file

      ! Volumes_file   = trim(Output_folder) // "/Other/Volumes.txt"
      ! P_iter_file    = trim(Output_folder) // "/Other/P_iter_info.txt"
      ! Iter_info_file = trim(Output_folder) // "/Other/Iter_info.txt"

      ! start timer
      !call cpu_time(start)
      start = omp_get_wtime()

      call log_solver_mode(Confined_flag)

      !Open a text file to store the volumes at each iteration
      ! open (newunit=fid_Volume,     file=Volumes_file,   action='write')
      ! open (newunit=fid_P_iter,     file=P_iter_file,    action='write')
      ! open (newunit=fid_Iter_info,  file=Iter_info_file, action='write')

      call iter_file_open(Output_folder)


      !Write first lines listing useful information
      !Volume_data = [ t, V_mobile_sum, V_trapped_sum, V_injected, V_mobile_local_max, V_trapped_local_max, V_dsln ]
      ! write (fid_Volume, '(A5,2X,7(A15,2X))') &
      ! & 'it', 't', 'V_mob_sum', 'V_trp_sum', 'V_inj', 'V_mob_lcl_max', 'V_trp_lcl_max', 'V_dsln'


      ! write (fid_Volume, '(A5,2X,8(A15,2X))') &
      ! write (fid_Volume, '(A5,2X,(A9,2X),7(A15,2X))') &
      ! write (fid_Volume, '(A5,1X,(A9,1X),8(A15,1X))') &
      ! ! write (fid_Volume, '(A5,1X,(A9,1X),7(A12,1X))') &
      ! & 'it', 't [days]', 'V_inj [m^3]', 'V_tot', 'V_mob_sum', 'V_trp_sum', 'V_dsln_sum', 'V_boundary', 'V_mob_max', 'V_trp_max'
      

      ! write(fid_Iter_info,'(A5,2X,7(A15,2X))') &
      ! & 'it', 't', 'dt_old', 'dt_ideal_old', 'err_h', 'err_P', 'max_h', 'max_P'

      !Initialise current profiles and arrays
      h     = h_init        ! Initial current profile
      h_res = h_res_init       ! Initial residually trapped profile
      h_old = h            ! Bookkeeping of previous profile, to check the sign of dh/dt
      dh_dt = 0._wp
      P     = P_init       ! Initial pressure profile
      Q     = 0._wp        ! Initial injection array

      D0    = B0 - H0      ! Reservoir thickness

      h_array     = 0._wp  ! Current Thickness profiles array
      h_res_array = 0._wp  ! Maximum Current thickness arrray
      P_array     = 0._wp  ! Pressure array

      h_array(:,:,0)     = h
      h_res_array(:,:,0) = h_res
      P_array(:,:,0)     = P

      ! Calculate the volumes corresponding to the initial profiles
      call initial_volumes(V_injected, V_tot_dsln_loss, V_mobile_array, V_trapped_array)


      V_boundary = 0._wp
      V_tot_dsln_loss = 0._wp

      plot_times = 0._wp

      P_iter_count = 0

      flux_idx = 0

      Max_mobile_thickness = 0._wp
      Max_pressure = 0._wp
      Min_pressure = 0._wp

      prev_rel_err_val = 0._wp
      dt_prev = 0._wp

      ! initialise variables
      t           = t0                    ! start time
      dt          = dt_init               ! initial time step
      dt_ideal    = dt
      plot_count  = 0                     ! Data output counter
      tplot       = target_plot_times(0)  !First time at which to plot


      call log_run_parameters(V_injected)

      ! Iteration counter
      it = 1


      !Calculate intial pressure (to avoid it having to be done twice by the
      !two timestepping branches for the first step)
      call generate_flux_array(Q,flux_idx,t,flux_change_flag)
      call change_of_flux(t,Q,h,h_res,P,P_old,P_iter_count,Confined_flag)


      ! -- Main calculation loop --------------------------------------------------------------

      do while ((t < tmax) .and. (it < itermax))

         !Check if there are any NaNs
         if ( any( isnan(h) ) ) then
            write(*,'(A, I0, A, F10.5)') 'h has NaN values! it = ', it, ', t = ', t
            write(*,*) 'Halting.'
            ERROR STOP
         end if

         if ( any( isnan(h_res) ) ) then
            write(*,'(A, I0, A, F10.5)') 'h_res has NaN values! it = ', it, ', t = ', t
            write(*,*) 'Halting.'
            ERROR STOP
         end if

         if ( any( isnan(P) ) ) then
            write(*,'(A, I0, A, F10.5)') 'P has NaN values! it = ', it, ', t = ', t
            write(*,*) 'Halting.'
            ERROR STOP
         end if
         
         ! Check if the current has passed through the basement.
         ! This only applies for the Confined case; skip the check if Unconfined
         if ( any( (h+h_res) > D0 ) ) then
            maxloc_idxs = maxloc(h+h_res) - 1 !Subtract one to convert to zero-indexing
            if ( Confined_flag ) then 
               print '("!!!  t = ", f10.2, ": at (", I0, ",", I0, ") &
               & max(h+h_res) = ", f10.5, ",  D0 = ", f10.5, ", h = ", f10.5, ", h_res = ", f10.5, " !!!")', &
               & t*Time_scale*seconds_to_days, &
               & maxloc_idxs, &
               & maxval(h+h_res)*Length_scale, &
               & D0(maxloc_idxs(0), maxloc_idxs(1))*Length_scale, &
               & h(maxloc_idxs(0), maxloc_idxs(1))*Length_scale, &
               & h_res(maxloc_idxs(0), maxloc_idxs(1))
               ! &  t, maxval(h+h_res), maxloc(h+h_res)
               !stop

               !If h has exceeded the reservoir depth, pull it back
               h = merge( h , D0 , h <= D0 )
               !If the residually trapped region is more than the gap between the mobile region and the
               !base of the reservoir, pull it back
               h_res = merge( h_res , D0 - h , (h+h_res)<=D0 )
               !Recalculate pressure
               P_old = P
               call Pressure_calculation(P, P_old, h, h_res, Q, P_iter_count, Confined_flag)
            end if
         end if


         ! Update the injection flux array
         call generate_flux_array(Q,flux_idx,t,flux_change_flag)
         if (flux_change_flag) then
            ! ! ! !Flux has changed, recalculate pressure
            call change_of_flux(t,Q,h,h_res,P,P_old,P_iter_count,Confined_flag)
         endif

         !Perform a calculation for this timestep
         call adaptive_timestepping(h_array, h_res_array, P_array, V_mobile_array, V_trapped_array, &
         & Max_mobile_thickness, Max_pressure, Min_pressure, h, h_res, h_old, dh_dt, P, &
         & t, dt_ideal, dt_new, tplot, plot_times, Q, V_injected, V_tot_dsln_loss, V_boundary, dt, plot_count, it, Confined_flag)

         !Increment iteration counter
         it = it + 1

         !Update dt for the next iteration
         dt = max(dt_min, dt_new)

         !If we've passed the final plot time then stop, otherwise prepare for the next loop
         if (plot_count .ge. np) then
            ! stop timer
            !call cpu_time(finish)
            finish = omp_get_wtime()
            runtime = finish - start

            call log_runtime(runtime)

            call routine_runtimes(runtime,Output_folder)

            call iter_file_close
      
            !Close record of Volumes
            ! close(fid_Volume)
            ! !Close Iteration info file
            ! close(fid_Iter_info)

            return

         end if

      end do


      call iter_file_close

      !Close record of Volumes
      ! close(fid_Volume)
      ! !Close record of Pressure iterations
      ! close(fid_P_iter)
      ! !Close Iteration info file
      ! close(fid_Iter_info)


   end subroutine solver_run


   ! ------------------------------------------------------------------------
   subroutine adaptive_timestepping(h_array, h_res_array, P_array, V_mobile_array, V_trapped_array, &
   &   Max_mobile_thickness, Max_pressure, Min_pressure, h, h_res, h_old, dh_dt, P, &
   &   t, dt_ideal, dt_new, tplot, plot_times, Q, V_injected, V_tot_dsln_loss, V_boundary, &
   &   dt_old, plot_count, iter_num, Confined_flag)
      ! This subroutine performs the calculation twice - from t to t+dt in one step, or in
      ! two half steps. The results are then kept and time advance if these two answers agree
      ! sufficiently, otherwise we will repeat the calculation with a smaller timestep (if possible)

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: Q
      logical, intent(in) :: Confined_flag
      real(wp), dimension(0:nx-1, 0:ny-1, 0:np-1), intent(inout) :: h_array, h_res_array, P_array
      real(wp), dimension(0:nx-1, 0:ny-1, 0:np-1), intent(inout) :: V_mobile_array, V_trapped_array
      real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) :: h, h_res, h_old, dh_dt, P
      real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) :: Max_mobile_thickness, Max_pressure, Min_pressure
      real(wp), dimension(0:np-1), intent(inout) :: plot_times
      real(wp), intent(inout) :: V_injected, V_tot_dsln_loss, V_boundary
      real(wp), intent(in) :: dt_old
      real(wp), intent(inout) :: dt_ideal
      real(wp), intent(out) :: t, tplot
      real(wp), intent(out) :: dt_new
      integer, intent(inout) :: plot_count
      integer, intent(in) :: iter_num
      real(wp), dimension(0:nx-1, 0:ny-1) :: htest1, htest2, h_res_test1, h_res_test2, Ptest1, Ptest2
      real(wp), dimension(0:nx-1, 0:ny-1) :: dh_dt_test1, dh_dt_test2, V_dsln_loss_1, V_dsln_loss_2
      real(wp) :: V_boundary_1, V_boundary_2
      real(wp) :: err_h, err_P, toterr, dt_2, dt_ideal_old
      integer, dimension(0:3) :: P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b

      integer :: iter_check = 1


      !Dummy variables for the one-step and two-step calculations, which will then be compared
      htest1 = h
      htest2 = h
      h_res_test1 = h_res
      h_res_test2 = h_res
      Ptest1 = P
      Ptest2 = P
      dh_dt_test1 = dh_dt
      dh_dt_test2 = dh_dt

      dt_2 = dt_old/2._wp

      !The previous step produced a target step size, dt_ideal. We may have ended up taking a smaller step,
      !dt, if we were sufficiently close to a plot time. We want to remember the remaining part of the
      !target step when deciding how big a step to take after this iteration.
      dt_ideal_old = max( dt_min, dt_ideal - dt_old )

      ! Do one full step with current stepsize (t -> t+dt)
      call timestep(htest1, Ptest1, h_res_test1, dh_dt_test1, Q, dt_old, P_iter_vals_1, V_dsln_loss_1, V_boundary_1, Confined_flag)

      if ( (modulo(iter_num,iter_check) .ne. 0) .and. (iter_num > 10) ) then

         !Keep the output of the latest run
         h_old             = h
         h                 = htest1
         h_res             = h_res_test1
         dh_dt             = dh_dt_test1
         P                 = Ptest1
         V_injected        = V_injected + sum(Q)*dt_old
         V_tot_dsln_loss   = V_tot_dsln_loss + sum(V_dsln_loss_1)*dx*dy
         V_boundary        = V_boundary + V_boundary_1

         !Update the time
         t = t + dt_old

         !Store the output from the current iteration and check if it's time to save a profile
         call save_iter_output(h_array, h_res_array, P_array, V_mobile_array, V_trapped_array, &
         & Max_mobile_thickness, Max_pressure, Min_pressure, h, h_res, P, &
         & t, tplot, plot_count, plot_times, V_injected, V_tot_dsln_loss, V_boundary, iter_num)

         !Set the new time step
         dt_ideal = dt_old

         !Set nominal error values
         err_h = 0._wp
         err_P = 0._wp
         P_iter_vals_2a(:) = 0
         P_iter_vals_2b(:) = 0

         ! ! Pressure iteration info
         ! write(fid_P_iter,'(2E12.4E3,1X,3(4I5,2X),1X,E12.4E3)') &
         ! & t, dt_old, P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b, err_P

         ! !General iteration info: iter, t, dt, err_h, err_P, max_h, max_P
         ! write(fid_Iter_info,'(I5,1X,7(E16.8E3,1X))') &
         ! & iter_num, t-dt_old, dt_old, dt_ideal_old, err_h, err_P, maxval(h), maxval(P)
         ! !Saving t-dt_old here since t is updated before this here, as opposed to below

         call iter_calculations(h, h_res, P, iter_num, t-dt_old, dt_old, dt_ideal, err_h, err_P, &
         & P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b)

      else

         ! Do two half steps with current stepsize (t -> t+dt/2 -> t+dt)
         call timestep(htest2, Ptest2, h_res_test2, dh_dt_test2, Q, &
         &                    dt_2, P_iter_vals_2a, V_dsln_loss_2, V_boundary_2, Confined_flag)
         call timestep(htest2, Ptest2, h_res_test2, dh_dt_test2, Q, &
         &                    dt_2, P_iter_vals_2b, V_dsln_loss_2, V_boundary_2, Confined_flag)


         ! Calculate the error between the two results
         err_h = err(htest1, htest2, nx, ny)
         err_P = err(Ptest1, Ptest2, nx, ny)
         toterr = err_h

         ! ! Pressure iteration info
         ! write(fid_P_iter,'(2E12.4E3,1X,3(4I5,2X),1X,E12.4E3)') &
         ! & t, dt_old, P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b, err_P

      !   !General iteration info: t, dt, err_h, err_P, max_h, max_P
      !    write(fid_Iter_info,'(I5,1X,7(E16.8E3,1X))') &
      !    & iter_num, t, dt_old, dt_ideal, err_h, err_P, maxval(htest2), maxval(Ptest2)

         call iter_calculations(htest2, h_res_test2, Ptest2, iter_num, t, dt_old, dt_ideal, err_h, err_P, &
         & P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b)


         ! Decide whether and how to update based on the error value
         if (toterr <= errmax_t_step) then
            !Keep the output of the latest run
            h_old             = h
            h                 = htest2
            h_res             = h_res_test2
            dh_dt             = dh_dt_test2
            P                 = Ptest2
            V_injected        = V_injected + sum(Q)*dt_old
            V_tot_dsln_loss   = V_tot_dsln_loss + sum(V_dsln_loss_2)*dx*dy
            V_boundary        = V_boundary + V_boundary_2

            !Update the time
            t = t + dt_old

            !Store the output from the current iteration and check if it's time to save a profile
            call save_iter_output(h_array, h_res_array, P_array, V_mobile_array, V_trapped_array, &
            & Max_mobile_thickness, Max_pressure, Min_pressure, h, h_res, P, &
            & t, tplot, plot_count, plot_times, V_injected, V_tot_dsln_loss, V_boundary, iter_num)

            !Set the new time step
            call ideal_timestep(dt_ideal, dt_old, toterr, iter_num)

         else
            !The error is too large - either try again with a reduced
            !step-size, or if that's not possible then keep what we have but
            !reduce the step size to the minimum

            call ideal_timestep(dt_ideal, dt_old, toterr, iter_num)

            if (dt_ideal > dt_min) then
               !Don't keep the output
               !Don't update the time

               ! call ideal_timestep(dt_ideal, dt_old, toterr, iter_num)
               ! dt_ideal = min( dt_ideal , dt_2 )

            else
               !Keep output
               h_old             = h
               h                 = htest2
               h_res             = h_res_test2
               dh_dt             = dh_dt_test2
               P                 = Ptest2
               V_injected        = V_injected + sum(Q)*dt_old
               V_tot_dsln_loss   = V_tot_dsln_loss + sum(V_dsln_loss_2)*dx*dy
               V_boundary        = V_boundary + V_boundary_2

               !Update time
               t = t + dt_old

               !Store the output from the current iteration and check if it's time to save a profile
               call save_iter_output(h_array, h_res_array, P_array, V_mobile_array, V_trapped_array, &
               & Max_mobile_thickness, Max_pressure, Min_pressure, h, h_res, P, &
               & t, tplot, plot_count, plot_times, V_injected, V_tot_dsln_loss, V_boundary, iter_num)

               !Set the new step size
               dt_ideal = dt_min
            end if

         end if

      end if


      !If the next plot time is within the intended step, go to the next plot time instead
      dt_new = min(dt_ideal, tplot - t)

   end subroutine adaptive_timestepping


   ! ------------------------------------------------------------------------
   subroutine ideal_timestep(dt_ideal, dt_old, rel_err_val, it_val)

      implicit none
      real(wp), intent(out) :: dt_ideal
      real(wp), intent(in) :: dt_old, rel_err_val
      integer, intent(in) :: it_val
      real(wp) :: beta_val_1, beta_val_2, alpha_val, k_val
      real(wp) :: eps_val, r_cur_val, r_prev_val

      k_val       = 2._wp

      ! !Elementary Controller
      ! alpha_val   = 0._wp
      ! beta_val_1  = 1._wp / k_val
      ! beta_val_2  = 0._wp

      ! !PI.4.2 Controller
      ! alpha_val   = 0._wp
      ! beta_val_1  = ( 3._wp/5._wp) / k_val
      ! beta_val_2  = (-1._wp/5._wp) / k_val

      !H211b digital filter
      alpha_val   = (1._wp/4._wp)
      beta_val_1  = (1._wp/4._wp) / k_val
      beta_val_2  = (1._wp/4._wp) / k_val

      !Target relative error value
      eps_val     = 0.1_wp * errmax_t_step

      !Recent measured rel. error values
      r_cur_val   = rel_err_val
      r_prev_val  = prev_rel_err_val



      ! dt_ideal = ( ( eps_val / (2._wp*rel_err_val) )**(1._wp/k_val) ) * dt_old

      if (it_val == 1) then
         dt_ideal = ( ( eps_val / r_cur_val )**(1._wp/k_val) ) * dt_old
      else
         dt_ideal = ( ( eps_val / r_cur_val  )**(beta_val_1) ) &
         &        * ( ( eps_val / r_prev_val )**(beta_val_2) ) &
         &        * ( ( dt_prev / dt_old )**(alpha_val) ) &
         &        * dt_old
      end if

      prev_rel_err_val = rel_err_val
      dt_prev = dt_old

   end subroutine ideal_timestep


   ! ------------------------------------------------------------------------
   subroutine timestep(ht1, Pt1, h_res, dh_dt, Q, dt_1, P_iter_vals, V_dissolution_loss, V_boundary, Confined_flag)
      ! This subroutine performs a single timestep calculation, updating h and P
      ! based on their current values. This is done in 4 steps of size dt_1/4,
      ! as a combination of a Predictor-Corrector method (since the h equation is
      ! nonlinear), and the ADI scheme (to solve the h equation in both x and y).
      ! The Pressure P, and where Convective Dissolution directly affects h, are
      ! updated at the start of each quarter step. The thickness of the
      ! residually trapped region, h_res, is updated after every half step.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: Q
      logical, intent(in) :: Confined_flag
      real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) :: ht1, Pt1, h_res, dh_dt
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) :: V_dissolution_loss
      real(wp), intent(out) :: V_boundary
      integer, dimension(0:3), intent(out) :: P_iter_vals
      real(wp), intent(in) :: dt_1
      real(wp), dimension(0:nx-1, 0:ny-1) :: ht2, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Pt2
      real(wp), dimension(0:nx-1, 0:ny-1) :: Forcing_term_h, h_initial, h_res_1, h_res_2
      real(wp), dimension(0:nx-1, 0:ny-1) :: F_a_H, F_a_V, sigma_a_V
      real(wp) :: dt_2, dt_4

      ! Intialise values
      V_dissolution_loss = 0._wp
      V_boundary = 0._wp

      P_iter_vals = 0

      dt_2 = dt_1/2._wp
      dt_4 = dt_1/4._wp

      h_initial = ht1

      ht2 = ht1
      Pt2 = Pt1

      h_res_1 = h_res
      h_res_2 = h_res


      ! -- x: t -> t + dt/4 -----------------------------------------------------------
      ! Old solution:            ht1, h_res_1, Pt1  at t
      ! Intermediate solution:   ht1, h_res_1, Pt1  at t
      ! Predicted solution:      ht2, h_res_2, Pt2  at t+dt/4

      ! Fluxes of ambient fluid within the residually trapped region
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht1, h_res_1, Pt1, dh_dt)
      ! Total flux term for h - combination of injection term and dissolution terms
      call calculate_h_forcing_term(Forcing_term_h, ht1, h_res_1, dh_dt, Q, sigma_a_V, F_a_H)

      ! Calculate coefficients for h equation and solve in x
      call h_coefficients(ht1, Pt1, h_res_1, dh_dt, sigma_a_V, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_x_solve(ht2, ht1, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_4)

      !Calculate Pressure
      ! call Pressure_solve(Pt2, Pt1, ht2, h_res_1, Q, P_iter_check, P_iter_vals(0))
      call Pressure_calculation(Pt2, Pt1, ht2, h_res_1, Q, P_iter_vals(0), Confined_flag)

      ! Calculate change in h_res
      dh_dt = (ht2-ht1)/(dt_4)
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht2, h_res_1, Pt2, dh_dt)
      call h_res_update(h_res_2,h_res_1,ht2,dh_dt,F_a_H,F_a_V,dt_4,Confined_flag)



      ! Calculate volume lost through East and West boundaries
      call calculate_boundary_flux(ht2, Pt2, dt_2, V_boundary, 'x')


      ! -- x: t -> t + dt/2 -----------------------------------------------------------
      ! Old solution:            ht1, h_res_1, Pt1  at t
      ! Intermediate solution:   ht2, h_res_2, Pt2  at t+dt/4
      ! New solution:            ht2, h_res_2, Pt2  at t+dt/2

      ! Fluxes of ambient fluid within the residually trapped region
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht2, h_res_2, Pt2, dh_dt)
      ! Total flux term for h - combination of injection term and dissolution terms
      call calculate_h_forcing_term(Forcing_term_h, ht2, h_res_2, dh_dt, Q, sigma_a_V, F_a_H)

      ! Calculate coefficients for h equation and solve in x
      call h_coefficients(ht2, Pt2, h_res_2, dh_dt, sigma_a_V, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_x_solve(ht2, ht1, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_2)

      !Calculate Pressure
      ! call Pressure_solve(Pt2, Pt1, ht2, h_res_2, Q, P_iter_check, P_iter_vals(1))
      call Pressure_calculation(Pt2, Pt1, ht2, h_res_2, Q, P_iter_vals(1), Confined_flag)


      ! Calculate change in h_res
      dh_dt = (ht2-ht1)/(dt_2)
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht2, h_res_2, Pt2, dh_dt)
      call h_res_update(h_res_2,h_res_1,ht2,dh_dt,F_a_H,F_a_V,dt_2,Confined_flag)

      ! Calculate volume lost to dissolution
      call calculate_dissolution_loss(V_dissolution_loss,ht2,ht1,h_res_2,h_res_1,F_a_V,dt_2)

      ! ! Calculate volume lost through East and West boundaries
      ! call calculate_boundary_flux(ht2, Pt2, dt_2, V_boundary, 'x')


      ! -- y: t+dt/2 -> t+dt/2 + dt/4 ------------------------------------------------
      ! Old solution:            ht2, h_res_2, Pt2  at t+dt/2
      ! Intermediate solution:   ht2, h_res_2, Pt2  at t+dt/2
      ! Predicted solution:      ht1, h_res_1, Pt1  at t+dt/2+dt/4

      ! Fluxes of ambient fluid within the residually trapped region
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht2, h_res_2, Pt2, dh_dt)
      ! Total flux term for h - combination of injection term and dissolution terms
      call calculate_h_forcing_term(Forcing_term_h, ht2, h_res_2, dh_dt, Q, sigma_a_V, F_a_H)

      ! Calculate coefficients for h equation and solve in y
      call h_coefficients(ht2, Pt2, h_res_2, dh_dt, sigma_a_V, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_y_solve(ht1, ht2, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_4)


      !Calculate Pressure
      ! call Pressure_solve(Pt1, Pt2, ht1, h_res_2, Q, P_iter_check, P_iter_vals(2))
      call Pressure_calculation(Pt1, Pt2, ht1, h_res_2, Q, P_iter_vals(2), Confined_flag)

      ! Calculate change in h_res
      dh_dt = (ht1-ht2)/(dt_4)
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht1, h_res_2, Pt1, dh_dt)
      call h_res_update(h_res_1,h_res_2,ht1,dh_dt,F_a_H,F_a_V,dt_4,Confined_flag)


      ! Calculate volume lost through North and South boundaries
      call calculate_boundary_flux(ht1, Pt1, dt_2, V_boundary, 'y')


      ! -- y: t+dt/2 -> t+dt/2 + dt/2 ------------------------------------------------
      ! Old solution:            ht2, h_res_2, Pt2  at t+dt/2
      ! Intermediate solution:   ht1, h_res_1, Pt1  at t+dt/2+dt/4
      ! New solution:            ht1, h_res_1, Pt1  at t+dt

      ! Fluxes of ambient fluid within the residually trapped region
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht1, h_res_1, Pt1, dh_dt)
      ! Total flux term for h - combination of injection term and dissolution terms
      call calculate_h_forcing_term(Forcing_term_h, ht1, h_res_1, dh_dt, Q, sigma_a_V, F_a_H)

      ! Calculate coefficients for h equation and solve in y
      call h_coefficients(ht1, Pt1, h_res_1, dh_dt, sigma_a_V, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_y_solve(ht1, ht2, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_2)

      !Calculate Pressure
      ! call Pressure_solve(Pt1, Pt2, ht1, h_res_1, Q, P_iter_check, P_iter_vals(3))
      call Pressure_calculation(Pt1, Pt2, ht1, h_res_1, Q, P_iter_vals(3), Confined_flag)


      ! Calculate change in h_res
      dh_dt = (ht1 - ht2)/(dt_2)
      call Residual_ambient_flux(F_a_H, F_a_V, sigma_a_V, ht1, h_res_1, Pt1, dh_dt)
      call h_res_update(h_res_1,h_res_2,ht1,dh_dt,F_a_H,F_a_V,dt_2,Confined_flag)

      ! Calculate volume lost to dissolution
      call calculate_dissolution_loss(V_dissolution_loss,ht1,ht2,h_res_1,h_res_2,F_a_V,dt_2)

      ! ! Calculate volume lost through North and South boundaries
      ! call calculate_boundary_flux(ht1, Pt1, dt_2, V_boundary, 'y')


      V_boundary = V_boundary + V_boundary

      ! Save value for output from this step
      h_res = h_res_1

      ! Update time derivative for the start of the next step
      dh_dt = (ht1 - h_initial)/(dt_1)

   end subroutine timestep


   ! ------------------------------------------------------------------------
   subroutine initial_volumes(V_injected, V_tot_dsln_loss, V_mobile_array, V_trapped_array)
      !This subroutine calculates the various volume values and arrays needed at the start
      !of the calculation based on the initial profiles

      implicit none
      real(wp), intent(out) :: V_injected, V_tot_dsln_loss
      real(wp), dimension(0:nx-1,0:ny-1,0:np-1), intent(out) :: V_mobile_array, V_trapped_array
      real(wp), dimension(0:nx-1,0:ny-1) :: V_mobile, V_trapped, V_mobile_dissolved , V_trapped_dissolved


      call calculate_volumes(h_init, h_res_init, V_mobile, V_trapped, V_mobile_dissolved, V_trapped_dissolved)

      ! Running total of CO2 injected
      V_injected      = sum(V_mobile + V_trapped + V_mobile_dissolved + V_trapped_dissolved)
      ! Running total of dissolved CO2 lost to the ambient
      V_tot_dsln_loss = sum(V_mobile_dissolved + V_trapped_dissolved)

      ! Arrays to store the local mobile and trapped volumes at each plot time
      V_mobile_array  = 0._wp
      V_trapped_array = 0._wp

      V_mobile_array(:,:,0)   = V_mobile
      V_trapped_array(:,:,0)  = V_trapped


   end subroutine initial_volumes

   ! ------------------------------------------------------------------------
   subroutine generate_flux_array(Q_array,flux_idx_cur,t_val,flux_change_flag)
      !This subroutine calculates the current flux profile, based on the
      !piecewise flux values in Q_flux_vals and the injection locations in
      !Q_inj_locs. Each row of Q_flux_vals specifies a new flux value for
      !each of the injection locations at a given time in t_flux_vals.
      !
      !n_inj_locs, n_flux_times, Q_inj_locs, Q_flux_vals, and t_flux_vals
      !are set in the -Injection_profiles- module.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) :: Q_array
      integer, intent(inout) :: flux_idx_cur
      real(wp), intent(in) :: t_val
      logical, intent(out) :: flux_change_flag
      real(wp), dimension(0:n_inj_locs-1) :: Q_vals
      integer :: k, flux_idx_new
      !real(wp) :: start_time, end_time

      ! start_time = omp_get_wtime()

      Q_array = 0._wp

      if (t_val < t_flux_vals(1)) then
         !Injection hasn't started yet
         flux_idx_new = 0
         Q_vals(:) = 0._wp

      else if (t_val >= t_flux_vals(n_flux_times-1)) then
         !After the last change in flux - use the index of the last specified value
         flux_idx_new = n_flux_times-1

      else
         !Calculate the index of the most recent specified flux time
         !in t_flux_vals before the current time t.
         !Start at the previous index value to save excessive counting.
         k = flux_idx_cur
         do while ( t_flux_vals(k) <= t_val )
            k = k+1
         end do
         flux_idx_new = k-1

      end if

      !Flag whether the flux index has changed
      if (flux_idx_new .ne. flux_idx_cur) then
         flux_change_flag = .True.
      else
         flux_change_flag = .False.
      endif
      !flux_change_flag = (flux_idx_new .ne. flux_idx_cur)

      !Update the flux index
      flux_idx_cur = flux_idx_new

      !Use the previous flux value - piecewise-constant interpolation
      !This reads a flux value for each injection location
      Q_vals = Q_flux_vals(flux_idx_cur,:)

      !Apply these flux values at their respective injection locations
      do k = 0, n_inj_locs-1
         Q_array( Q_inj_locs(k,0) , Q_inj_locs(k,1) ) = Q_vals(k)
      end do

      ! end_time = omp_get_wtime()
      ! ! Add runtime to running total
      ! Gen_flux_time = Gen_flux_time + (end_time - start_time)
      ! ! Increase count of calls
      ! Gen_flux_calls = Gen_flux_calls + 1

   end subroutine generate_flux_array

   ! ------------------------------------------------------------------------
   subroutine Residual_ambient_flux(F_H, F_V, sigma, h, h_res, P, dh_dt)
      ! This function calculates the term due to ambient fluid moving between
      ! vertical columns within the residually trapped region, i.e.
      ! H+h <= z <= H+h+h_res:
      !  F_H = int_{H+h}^{H+h+h_res} ( div(u_a) ) dz
      !      = int_{H+h}^{H+h+h_res} ( -M * k * grad(P) ) dz
      ! This can be rewritten as
      ! F_H  = M*[ -div( int_{H+h}^{H+h+h_res} k dz * grad P )
      !         + k|_{H+h+h_res} grad(P) \cdot grad(h_res)
      !         + (k|_{H+h+h_res} - k|_{H+h}) grad(P) \cdot grad(H+h) ]
      !
      ! This integral is necessarily zero wherever there is no residually
      ! trapped fluid, i.e. where h_res = 0.
      !
      ! Also calculated is sigma, which is the sign of the vertical flux of
      ! ambient fluid through the lower boundary of this region within each
      ! vertical column:
      ! sigma = sign( phi|_{H+h}*(1-s_a_i-s_c_r)*dh/dt - F_H )

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) :: F_H, F_V, sigma
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: h, h_res, P, dh_dt
      real(wp), dimension(0:nx-1, 0:ny-1) :: k_int_temp, k_int_res, k_res, k_h
      real(wp), dimension(0:nx-1, 0:ny-1) :: Dx_P, Dx_H_h, Dx_h_res, Dx_k_int_res
      real(wp), dimension(0:nx-1, 0:ny-1) :: Dy_P, Dy_H_h, Dy_h_res, Dy_k_int_res, D2_P
      real(wp), dimension(0:nx-1, 0:ny-1) :: phi_h
      logical :: extrapolate_flag
      !real(wp) :: start_time, end_time

      ! start_time = omp_get_wtime()
    
      !Initialise arrays
      sigma           = 1._wp
      F_H             = 0._wp
      F_V             = 0._wp
      k_int_temp      = 0._wp
      k_int_res       = 0._wp
      k_res           = 0._wp
      k_h             = 0._wp
      phi_h           = 0._wp
      Dx_P            = 0._wp
      Dx_H_h          = 0._wp
      Dx_h_res        = 0._wp
      Dx_k_int_res    = 0._wp
      Dy_P            = 0._wp
      Dy_H_h          = 0._wp
      Dy_h_res        = 0._wp
      Dy_k_int_res    = 0._wp
      D2_P            = 0._wp

      !Calculate the integral of the permeability over the trapped region
      extrapolate_flag = .True.
      call vertical_interpolation(perm_CInt, h      , D0, k_int_temp, extrapolate_flag)
      call vertical_interpolation(perm_CInt, h+h_res, D0, k_int_res , extrapolate_flag)
      k_int_res = k_int_res - k_int_temp
      !Scale by relative permeability value
      k_int_res = krw_residual * k_int_res

      !Calculate the permeability at the upper and lower boundaries
      call get_permeability_value(h      , k_h  )
      call get_permeability_value(h+h_res, k_res)
      !Scale by relative permeability values
      k_h = krn_mobile * k_h
      k_res = krw_residual * k_res
      
      !Calculate the porosity at the mobile interface (z=H0+h)
      call get_porosity_value(h, phi_h)

      !Calculate the gradients involved in F
      call Array_gradient(P        , Dx_P        , Dy_P        )
      call Array_gradient(H0+h     , Dx_H_h      , Dy_H_h      )
      call Array_gradient(h_res    , Dx_h_res    , Dy_h_res    )
      call Array_gradient(k_int_res, Dx_k_int_res, Dy_k_int_res)

      !Calculate the Laplacian of P
      call Array_Laplacian(P,D2_P)

      !Combine these all into F as described above
      F_H = M * ( &
      & - k_int_res * D2_P &
      & - ( Dx_k_int_res * Dx_P + Dy_k_int_res * Dy_P ) &
      & + k_res           * ( Dx_P * Dx_h_res + Dy_P * Dy_h_res ) &
      & + ( k_res - k_h ) * ( Dx_P * Dx_H_h   + Dy_P * Dy_H_h   ) &
         )

      !Set F_H = 0. outside of h_res > 0.
      F_H = merge( F_H , 0._wp , h_res > 0._wp )

      !Vertical flux of saturated ambient fluid out of the CO2 region  
      F_V = phi_h*(1._wp-s_a_i-s_c_r)*dh_dt - F_H

      ! Calculate sigma
      sigma(:,:) = 1._wp
      sigma = sign(sigma , F_V)

    
      ! end_time = omp_get_wtime()
      ! ! Add runtime to running total
      ! Res_amb_time = Res_amb_time + (end_time - start_time)
      ! ! Increase count of calls
      ! Res_amb_calls = Res_amb_calls + 1

   end subroutine Residual_ambient_flux

   ! ------------------------------------------------------------------------
   subroutine Array_gradient(F,DF_x,DF_y)

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: F
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) :: DF_x, DF_y
      integer :: i, j
      !real(wp) :: start_time, end_time

      ! start_time = omp_get_wtime()

      DF_x = 0._wp
      DF_y = 0._wp

      do i = 1,nx-2
         !I can ignore the edges for what I need
         DF_x(i,:) = ( F(i+1,:) - F(i-1,:) )/(2*dx)
      end do

      do j = 1,ny-2
         DF_y(:,j) = ( F(:,j+1) - F(:,j-1) )/(2*dy)
      end do
    
   !  end_time = omp_get_wtime()
   !  ! Add runtime to running total
   !  Array_grad_time = Array_grad_time + (end_time - start_time)
   !  ! Increase count of calls
   !  Array_grad_calls = Array_grad_calls + 1

   end subroutine Array_gradient

   ! ------------------------------------------------------------------------
   subroutine Array_Laplacian(F,D2F)

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: F
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) :: D2F
      integer :: i, j
      !real(wp) :: start_time, end_time

      ! start_time = omp_get_wtime()

      D2F = 0._wp

      !d2/dx^2 first
      do i = 1,nx-2
         !I can ignore the edges for what I need
         D2F(i,:) = ( F(i+1,:) - 2*F(i,:) + F(i-1,:) )/(dx**2)
      end do

      !Add d2/dy^2 to this
      do j = 1,ny-2
         D2F(:,j) = D2F(:,j) + &
         &      ( F(:,j+1) - 2*F(:,j) + F(:,j-1) )/(dy**2)
      end do
    
      ! end_time = omp_get_wtime()
      ! ! Add runtime to running total
      ! Array_Lap_time = Array_Lap_time + (end_time - start_time)
      ! ! Increase count of calls
      ! Array_Lap_calls = Array_Lap_calls + 1

   end subroutine Array_Laplacian

   ! ------------------------------------------------------------------------
   subroutine change_of_flux(t,Q,h,h_res,P,P_old,P_iter_count,Confined_flag)

      implicit none
      real(wp), intent(in) :: t
      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: Q, h, h_res
      real(wp), dimension(0:nx-1,0:ny-1), intent(inout) :: P
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: P_old
      integer, intent(out) :: P_iter_count
      logical, intent(in) :: Confined_flag
      character(100) :: flux_str


      !Flux has changed, recalculate pressure

      write(flux_str,'(A,F9.2,A)') '[ Flux change at t = ', t * Time_scale * seconds_to_days, ' days ]'

      call log_details_line(' ','')
      call log_details_line('*',flux_str)
      ! call log_details_line(' ','')
      write(*,'(A)') ' - Injection flux updated'
      write(*,'(6X,A,F10.5,1X,A)') 'Q_inj_total = ', sum(Q) * Flux_scale , '[m^3 s^-1]'

      write(*,'(A)') ' - Recalculating pressure'
      P_old = P
      call Pressure_calculation(P, P_old, h, h_res, Q, P_iter_count, Confined_flag)
      call P_flux_change_record(t, P_iter_count)

      write(*,'(A)') ' - Pressure calculation complete'

      call log_details_line('*','')
      call log_details_line(' ','')
      call log_title_row




   end subroutine change_of_flux





end module CO2GraVISim_solver
