module CO2GraVISim_solver
   !This is the main solver for CO2GraVISim, calculating the mobile current thickness h, 
   !the resdiually trapped current thickness h_res, and the ambient pressure Pa, for 
   !given inputs and at the specified plot times.
   

   use CO2GraVISim_global !wp, iter_max, tmax, dt_init, dt_min, errmax_t_step, errmax_P_iter, P_iter_omega, t0
   use CO2GraVISim_input_parameters !H0, B0, D0, perm_h, poro_h
   use CO2GraVISim_injection_profiles !n_inj_locs, n_flux_times, Q_inj_locs, Q_flux_vals, t_flux_vals

   implicit none

contains

   ! ------------------------------------------------------------------------
   subroutine solver_run(h_array,h_res_array,P_array,plot_times,target_plot_times)

      implicit none
      real(wp), dimension(0:nx-1,0:ny-1,0:np-1), intent(out) :: h_array, h_res_array, P_array
      real(wp), dimension(0:nx-1,0:ny-1) :: h, h_res, h_old, P, Q
      real(wp), dimension(0:np-1), intent(out) :: plot_times
      real(wp), dimension(0:np-1), intent(in) :: target_plot_times
      real(wp) :: t, dt, dt_ideal, dt_new, tplot, start, finish, V_injected
      integer  :: it, plot_count, flux_idx

      !output file for the various calculated volumes after each successful iteration
      character(50)  ::  Volumes_file = "./Output/Other/Volumes.txt"

      ! start timer
      call cpu_time(start)

      !Open a text file to store the volumes at each iteration
      open (99, file=Volumes_file, action='write')

      !Initialise current profiles and arrays
      h     = 0._wp        ! Initial current profile
      h_res = 0._wp        ! Initial residually trapped profile
      h_old = h            ! Bookkeeping of previous profile, to check the sign of dh/dt
      P     = 0._wp        ! Initial pressure profile
      Q     = 0._wp        ! Initial injection array

      D0    = B0 - H0      ! Reservoir thickness
      V_injected = 0._wp   ! Running total of CO2 injected

      h_array     = 0._wp  ! Current Thickness profiles array
      h_res_array = 0._wp  ! Maximum Current thickness arrray
      P_array     = 0._wp  ! Pressure array


      ! initialise variables
      t           = t0                    ! start time
      dt          = dt_init               ! initial time step
      plot_count  = 0                     ! Data output counter
      tplot       = target_plot_times(0)  !First time at which to plot

      write (*, *) 'variables intialised'
      print '("Initial time step is dt = ",e8.3 )', dt


      !Display parameter values used
      print '("-------------------------------------")'
      print '(" nx         = ", i4 )', nx
      print '(" ny         = ", i4 )', ny
      print '(" dx         = ", f8.3 )', dx
      print '(" dy         = ", f8.3 )', dy
      print '(" M          = ", f8.3 )', M
      print '(" Gamma_val  = ", f8.3 )', Gamma_val
      print '(" s_c_r      = ", f8.3 )', s_c_r
      print '(" s_a_i      = ", f8.3 )', s_a_i
      print '(" C_sat      = ", f8.3 )', C_sat
      print '(" q_dissolve = ", e8.3 )', q_dissolve
      print '("-------------------------------------")'

      ! Iteration counter
      it = 1

      print '("-------------------------------------")'


      ! -- Main calculation loop --------------------------------------------------------------

      do while ((t < tmax) .and. (it < itermax))

         ! Check if the current has passed through the basement
         if ( any( h > 1.01_wp * D0 ) ) then
            print '("!!!  t = ", f10.5, ", h_max = ", f10.5, ", P_max = ", e10.5, " !!!")', &
            &  t, maxval(h), maxval(abs(P))
            stop
         end if

         ! Update the injection flux array
         call generate_flux_array(Q,flux_idx,t)

         !Perform a calculation for this timestep
         call adaptive_timestepping(h_array, h_res_array, P_array, h, h_res, h_old, P, &
         & t, dt_ideal, dt_new, tplot, target_plot_times, plot_times,&
         & Q, V_injected, dt, dt_min, plot_count)

         !Increment iteration counter
         it = it + 1

         !Update dt for the next iteration
         dt = max(dt_min, dt_new)

         !If we've passed the final plot time then stop, otherwise prepare for the next loop
         if (plot_count .ge. np) then
            ! stop timer
            call cpu_time(finish)

            print '("-------------------------------------")'
            print '("Time taken to run = ",f8.3," seconds.")', finish - start
            print '("-------------------------------------")'

            return

         end if

      end do

      !Close record of Volumes
      close(99)

   end subroutine solver_run


   ! ------------------------------------------------------------------------
   subroutine adaptive_timestepping(h_array, h_res_array, P_array,&
   & h, h_res, h_old, P, t, dt_ideal, dt_new, tplot, target_plot_times, plot_times, &
   & Q, V_injected, dt_old, dt_min_val, plot_count)
      ! This subroutine performs the calculation twice - from t to t+dt in one step, or in
      ! two half steps. The results are then kept and time advance if these two answers agree
      ! sufficiently, otherwise we will repeat the calculation with a smaller timestep (if possible)

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: Q
      real(wp), dimension(0:nx-1, 0:ny-1, 0:np-1), intent(inout) :: h_array, h_res_array, P_array
      real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) :: h, h_res, h_old, P
      real(wp), dimension(0:np-1), intent(inout) :: plot_times
      real(wp), dimension(0:np-1), intent(in) :: target_plot_times
      real(wp), intent(inout) :: V_injected
      real(wp), intent(in) :: dt_old, dt_min_val
      real(wp), intent(inout) :: dt_ideal
      real(wp), intent(out) :: t, tplot
      real(wp), intent(out) :: dt_new
      real(wp), dimension(0:nx-1, 0:ny-1) :: htest1, htest2, h_res_test1, h_res_test2, Ptest1, Ptest2
      real(wp) :: toterr, dt_2
      integer, intent(inout) :: plot_count


      !Dummy variables for the one-step and two-step calculations, which will then be compared
      htest1 = h
      htest2 = h
      h_res_test1 = h_res
      h_res_test2 = h_res
      Ptest1 = P
      Ptest2 = P

      dt_2 = dt_old/2._wp

      ! Do one full step with current stepsize (t -> t+dt)
      call timestep(htest1, Ptest1, h_res_test1, h_old,  Q, dt_old)
      ! Do two half steps with current stepsize (t -> t+dt/2 -> t+dt)
      call timestep(htest2, Ptest2, h_res_test2, h_old,  Q, dt_2)
      call timestep(htest2, Ptest2, h_res_test2, htest2, Q, dt_2)

      ! Calculate the error between the two results
      toterr = err(htest1, htest2)

      ! Decide whether and how to update based on the error value
      if (toterr <= errmax_t_step) then
         !Keep the output of the latest run
         h_old       = h
         h           = htest2
         h_res       = h_res_test2
         P           = Ptest2
         V_injected  = V_injected + sum(Q)*dt_old

         !Update the time
         t = t + dt_old

         !Store the output from the current iteration and check if it's time to save a profile
         call save_iter_output(h_array, h_res_array, P_array,&
         & h, h_res, P, t, tplot, plot_count, target_plot_times, plot_times, V_injected)

         !Set the new time step
         if (toterr < errmax_t_step/2._wp) then
            !Try increasing the step size
            dt_ideal = 2.1_wp*dt_old
         else
            dt_ideal = dt_old
         end if

      else
         !The error is too large - either try again with a reduced
         !step-size, or if that's not possible then keep what we have but
         !reduce the step size to the minimum

         if (dt_old > 2._wp*dt_min_val) then
            !Don't keep the output
            !Don't update the time
            !Halve dt
            dt_ideal = dt_2
         else
            !Keep output
            h_old = h
            h     = htest2
            h_res = h_res_test2
            P     = Ptest2
            V_injected = V_injected + sum(Q)*dt_old

            !Update time
            t = t + dt_old

            !Store the output from the current iteration and check if it's time to save a profile
            call save_iter_output(h_array, h_res_array, P_array,&
            & h, h_res, P, t, tplot, plot_count, target_plot_times, plot_times, V_injected)

            !Set the new step size
            dt_ideal = dt_min_val
         end if

      end if

      !If the next plot time is within the intended step, go to the next plot time instead
      dt_new = min(dt_ideal, tplot - t)

   end subroutine adaptive_timestepping


   ! ------------------------------------------------------------------------
   subroutine timestep(ht1, Pt1, h_res, h_old, Q, dt_1)
      ! This subroutine performs a single timestep calculation, updating h and P
      ! based on their current values. This is done in 4 steps of size dt_1/4,
      ! as a combination of a Predictor-Corrector method (since the h equation is
      ! nonlinear), and the ADI scheme (to solve the h equation in both x and y).
      ! The Pressure P, and where Convective Dissolution directly affects h, are
      ! updated at the start of each quarter step. The thickness of the
      ! residually trapped region, h_res, is updated after every half step.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: h_old, Q
      real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) :: ht1, Pt1, h_res
      real(wp), intent(in) :: dt_1
      real(wp), dimension(0:nx-1, 0:ny-1) :: ht2, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Pt2
      real(wp), dimension(0:nx-1, 0:ny-1) :: Forcing_term_h, h_cur, Dissolution_array
      real(wp) :: dt_2, dt_4
      integer, dimension(0:3) :: P_iter_vals

      dt_2 = dt_1/2._wp
      dt_4 = dt_1/4._wp

      h_cur = ht1

      ! -- x: t -> t + dt/4 -----------------------------------------------------------
      ! Points where Convective Dissolution occurs for h
      Dissolution_array = merge(1._wp, 0._wp, (ht1 > 0._wp) .and. (h_res <= 0._wp) )
      !Total flux term - combination of injection term and convective dissolution term
      Forcing_term_h = Q/(dx*dy) - q_dissolve*Dissolution_array

      call Pressure_solve(Pt2, Pt1, ht1, Q, P_iter_vals(0))
      call h_coefficients(ht1, Pt2, h_res, h_old, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_x_solve(ht2, ht1, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_4)


      ! -- x: t -> t + dt/2 -----------------------------------------------------------
      Dissolution_array = merge(1._wp, 0._wp, (ht1 > 0._wp) .and. (h_res <= 0._wp) )
      Forcing_term_h = Q/(dx*dy) - q_dissolve*Dissolution_array

      call Pressure_solve(Pt2, Pt1, ht2, Q, P_iter_vals(1))
      call h_coefficients(ht2, Pt2, h_res, h_old, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_x_solve(ht2, ht1, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_2)

      call h_res_update(h_res,ht2,ht1,dt_2)



      ! -- y: t+dt/2 -> t+dt/2 + dt/4 ------------------------------------------------
      Dissolution_array = merge(1._wp, 0._wp, (ht1 > 0._wp) .and. (h_res <= 0._wp) )
      Forcing_term_h = Q/(dx*dy) - q_dissolve*Dissolution_array

      call Pressure_solve(Pt1, Pt2, ht2, Q, P_iter_vals(2))
      call h_coefficients(ht2, Pt1, h_res, h_cur, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_y_solve(ht1, ht2, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_4)



      ! -- y: t+dt/2 -> t+dt/2 + dt/2 ------------------------------------------------
      Dissolution_array = merge(1._wp, 0._wp, (ht1 > 0._wp) .and. (h_res <= 0._wp) )
      Forcing_term_h = Q/(dx*dy) - q_dissolve*Dissolution_array

      call Pressure_solve(Pt1, Pt2, ht1, Q, P_iter_vals(3))
      call h_coefficients(ht1, Pt1, h_res, h_cur, cxp_h, cxm_h, cyp_h, cym_h, ct_h)
      call ADI_y_solve(ht1, ht2, cxp_h, cxm_h, cyp_h, cym_h, ct_h, Forcing_term_h, dt_2)

      call h_res_update(h_res,ht1,ht2,dt_2)


   end subroutine timestep


   ! ------------------------------------------------------------------------
   subroutine Pressure_solve(P_new, P_old, h_cur, Q, iter_count)
      ! This subroutine calculates the coefficients that go into the Poisson equation
      ! for the ambient pressure, which is then solved by Pressure_SOR.
      ! This equation is of the form
      ! Div( A * grad(P) ) = - Div( B * grad(G) ) - Q
      ! This is ultimately discretised as
      !     cxm(i,j)*P(i-1,j  ) + cx0(i,j)*P(i,j) + cxp(i,j)*P(i+1,j  )
      !   + cym(i,j)*P(i  ,j-1) + cy0(i,j)*P(i,j) + cyp(i,j)*P(i  ,j+1)
      !   = Forcing_term(i,j)

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  P_old, h_cur, Q
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  P_new
      real(wp), dimension(0:nx-1, 0:ny-1) ::  Axavg, Cxp, Cxm, Cx0
      real(wp), dimension(0:nx-1, 0:ny-1) ::  Ayavg, Cyp, Cym, Cy0
      real(wp), dimension(0:nx-1, 0:ny-1) :: Afn, Bfn, Gfn, Forcing_term
      integer, intent(out) :: iter_count

      Afn = h_cur*perm_h + M*(D0 - h_cur)*perm_h
      Bfn = Gamma_val*h_cur*perm_h
      Gfn = H0 + h_cur

      Axavg(:, :) = 0._wp
      Ayavg(:, :) = 0._wp

      ! Calculate the values of A at the interfaces between the grid points
      Axavg(0:nx-2, 0:ny-1) = (Afn(1:nx-1, 0:ny-1) + Afn(0:nx-2, 0:ny-1))/(2._wp)
      Ayavg(0:nx-1, 0:ny-2) = (Afn(0:nx-1, 1:ny-1) + Afn(0:nx-1, 0:ny-2))/(2._wp)

      ! x coefficients
      Cxp(:, :) = 0._wp
      Cxm(:, :) = 0._wp
      Cx0(:, :) = 0._wp

      Cxp(1:nx-2, 1:ny-2) = Axavg(1:nx-2, 1:ny-2)/(dx**2)
      Cxm(1:nx-2, 1:ny-2) = Axavg(0:nx-3, 1:ny-2)/(dx**2)
      Cx0(1:nx-2, 1:ny-2) =  -Cxp(1:nx-2, 1:ny-2) - Cxm(1:nx-2, 1:ny-2)

      ! y coefficients
      Cyp(:, :) = 0._wp
      Cym(:, :) = 0._wp
      Cy0(:, :) = 0._wp

      Cyp(1:nx-2, 1:ny-2) = Ayavg(1:nx-2, 1:ny-2)/(dy**2)
      Cym(1:nx-2, 1:ny-2) = Ayavg(1:nx-2, 0:ny-3)/(dy**2)
      Cy0(1:nx-2, 1:ny-2) =  -Cyp(1:nx-2, 1:ny-2) - Cym(1:nx-2, 1:ny-2)


      ! Build the forcing term up from B, G, and Q
      Forcing_term(:, :) = 0._wp

      Forcing_term(1:nx - 2, 1:ny - 2) = Q(1:nx - 2, 1:ny - 2)/(dx*dy) &
      & + (Bfn(2:nx-1, 1:ny-2) + Bfn(1:nx-2, 1:ny-2))*(Gfn(2:nx-1, 1:ny-2) - Gfn(1:nx-2, 1:ny-2))/(2._wp*(dx**2)) &
      & - (Bfn(1:nx-2, 1:ny-2) + Bfn(0:nx-3, 1:ny-2))*(Gfn(1:nx-2, 1:ny-2) - Gfn(0:nx-3, 1:ny-2))/(2._wp*(dx**2)) &
      & + (Bfn(1:nx-2, 2:ny-1) + Bfn(1:nx-2, 1:ny-2))*(Gfn(1:nx-2, 2:ny-1) - Gfn(1:nx-2, 1:ny-2))/(2._wp*(dy**2)) &
      & - (Bfn(1:nx-2, 1:ny-2) + Bfn(1:nx-2, 0:ny-3))*(Gfn(1:nx-2, 1:ny-2) - Gfn(1:nx-2, 0:ny-3))/(2._wp*(dy**2))


      !Solve for the new pressure
      call Pressure_SOR(P_new, P_old, cxp, cxm, cx0, cyp, cym, cy0, Forcing_term, iter_count)

      return
   end subroutine Pressure_solve


   ! ------------------------------------------------------------------------
   subroutine Pressure_SOR(P_new, P_old, cxp, cxm, cx0, cyp, cym, cy0, Forcing_term, iter_count)
      !This subroutine solves the discretised Poisson equation for the pressure:
      !     cxm(i,j)*P(i-1,j  ) + cx0(i,j)*P(i,j) + cxp(i,j)*P(i+1,j  )
      !   + cym(i,j)*P(i  ,j-1) + cy0(i,j)*P(i,j) + cyp(i,j)*P(i  ,j+1)
      !   = Forcing_term(i,j)
      !These coefficients are adjusted if Neumann boundary conditions are imposed on
      !the respective edges of the domain.
      !This equation is solved iteratively via Successive Over-Relaxation (SOR).
      !The relaxation parameter P_iter_omega and the convergence threshold errmax_P_iter
      !are set in the -global- module.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  P_old, cxp, cxm, cx0, cyp, cym, cy0, Forcing_term
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  P_new
      real(wp), dimension(0:nx-1, 0:ny-1) ::  P_prev, dxy0, dxm, dxp, dym, dyp
      real(wp) :: error_val
      integer ::  i, j
      integer, intent(out) :: iter_count

      P_prev = P_old

      !Coefficients in the pressure equation below. Introduce new arrays so that entries can be
      !modified if Neumann BCs are used.
      dxy0  = cx0 + cy0
      dxm   = cxm
      dxp   = cxp
      dym   = cym
      dyp   = cyp

      !If we have Neumann BCs, then we need to modify these to
      !include the contributions from P(0,:)=P(1,:), etc.

      !North side
      if (P_BC_params(0) == 2) then
         !Neumann BC
         !P(:,ny-1) = P(:,ny-2), so add (i,j+1) coefficient into the (i,j) coefficient
         dxy0(:,ny-2) = dxy0(:,ny-2) + dyp(:,ny-2)
         dyp(:,ny-2) = 0._wp
      end if

      !East side
      if (P_BC_params(1) == 2) then
         !Neumann BC
         !P(nx-1,:) = P(nx-2,:), so add (i+1,j) coefficient into the (i,j) coefficient
         dxy0(nx-2,:) = dxy0(nx-2,:) + dxp(nx-2,:)
         dxp(nx-2,:) = 0._wp
      end if

      !South side
      if (P_BC_params(2) == 2) then
         !Neumann BC
         !P(:,0) = P(:,1), so add (i,j-1) coefficient into the (i,j) coefficient
         dxy0(:,1) = dxy0(:,1) + dym(:,1)
         dym(:,1) = 0._wp
      end if

      !West side
      if (P_BC_params(3) == 2) then
         !Neumann BC
         !P(0,:) = P(1,:), so add (i-1,j) coefficient into the (i,j) coefficient
         dxy0(1,:) = dxy0(1,:) + dxm(1,:)
         dxm(1,:) = 0._wp
      end if


      !Apply the BCs to P_prev before we start
      call Pressure_BCs(P_prev)


      !Solve the equation for the pressure iteratively via SOR, iterating until the solution
      !has converged sufficiently.
      P_new = P_prev

      !Set the initial error value to something above the threshold
      error_val = 2*errmax_P_iter

      iter_count = 0

      do while (error_val > errmax_P_iter)
         do i = 1, nx-2
            do j = 1, ny-2

               P_new(i, j) = (1 - P_iter_omega)*P_prev(i, j) + &
               & (P_iter_omega/dxy0(i, j))* &
               &            (- dxm(i, j)*P_new(i - 1, j    ) - dxp(i, j)*P_new(i + 1, j    ) &
               &             - dym(i, j)*P_new(i    , j - 1) - dyp(i, j)*P_new(i    , j + 1) &
               &         - Forcing_term(i, j))

            end do
         end do

         !Apply the BCs to the edges of P_new
         call Pressure_BCs(P_new)

         error_val = err(P_prev, P_new)
         P_prev = P_new
         iter_count = iter_count + 1
      end do

      return

   end subroutine Pressure_SOR


   ! ------------------------------------------------------------------------
   subroutine Pressure_BCs(P)
      !This subroutine imposed the boundary conditions -- either Dirichlet or
      !Neumann -- on the edges of the Pressure solution array

      implicit none
      real(wp), dimension(0:nx-1,0:ny-1), intent(inout):: P


      !North side
      if (P_BC_params(0) == 1) then
         !Dirichlet BC
         P(:, ny-1) = 0._wp
      else if (P_BC_params(0) == 2) then
         !Neumann BC
         P(:, ny-1) = P(:, ny-2)
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if

      !East side
      if (P_BC_params(1) == 1) then
         !Dirichlet BC
         P(nx-1, :) = 0._wp
      else if (P_BC_params(1) == 2) then
         !Neumann BC
         P(nx-1, :) = P(nx-2, :)
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if

      !South side
      if (P_BC_params(2) == 1) then
         !Dirichlet BC
         P(:, 0) = 0._wp
      else if (P_BC_params(2) == 2) then
         !Neumann BC
         P(:, 0) = P(:, 1)
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if

      !West side
      if (P_BC_params(3) == 1) then
         !Dirichlet BC
         P(0, :) = 0._wp
      else if (P_BC_params(3) == 2) then
         !Neumann BC
         P(0, :) = P(1, :)
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if


   end subroutine


   ! ------------------------------------------------------------------------
   subroutine h_coefficients(h_cur, P_cur, h_res, h_prev, cxp, cxm, cyp, cym, ct)
      !This subroutine calculates the coefficients to be used as part of the ADI method
      !to solve for the current thickness h.
      !The equation for h is of the form
      !            ct*dh/dt + Div( v * h ) = Div( Diff * grad(h) ) + Q
      !Here the discretised forms of v and Diff are calculated and passed to
      !Ilin_Coefficients, which calculates the corresponding coefficients for the ADI method.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  h_cur, P_cur, h_res, h_prev
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  cxp, cxm, cyp, cym, ct
      real(wp), dimension(0:nx-2, 0:ny-2) :: Diff_x, vi_x, Diff_y, vi_y
      real(wp), dimension(0:nx-1, 0:ny-1) :: Afn, Bfn, R_trap_array, h_inc_array, h_dec_array, u_array
      real(wp), dimension(0:nx-1, 0:ny-1) :: ct_inc, ct_dec


      ! -- Calculate ct --------------------------------------------------------------------------------

      !Locations where h is increasing, h is decreasing,
      !where residually trapped CO2 is present, and where it isn't present
      h_inc_array    = merge(1._wp, 0._wp, h_cur >= h_prev)
      h_dec_array    = merge(1._wp, 0._wp, h_cur <  h_prev)
      R_trap_array   = merge(1._wp, 0._wp, h_res >  0._wp )
      u_array        = merge(1._wp, 0._wp, h_res <= 0._wp )


      ! dh/dt prefactor - depends on sign of dh/dt, and where h_res>0
      ct_inc = poro_h * h_inc_array *( (1._wp - C_sat)*(1._wp - s_a_i - s_c_r)*R_trap_array &
      &                                      + (1._wp - s_a_i + C_sat*s_a_i)*u_array )

      ct_dec = poro_h * h_dec_array *( (1._wp - s_a_i - s_c_r) )

      ct = ct_inc + ct_dec


      ! -- Calculate v and Diff ---------------------------------------------------------------------------------------
      Afn = Gamma_val*h_cur*perm_h
      Bfn = -(P_cur + Gamma_val*H0)

      vi_x = (perm_h(1:nx-1, 0:ny-2) + perm_h(0:nx-2, 0:ny-2))*(Bfn(1:nx-1, 0:ny-2) - Bfn(0:nx-2, 0:ny-2))/(2._wp*dx)

      vi_y = (perm_h(0:nx-2, 1:ny-1) + perm_h(0:nx-2, 0:ny-2))*(Bfn(0:nx-2, 1:ny-1) - Bfn(0:nx-2, 0:ny-2))/(2._wp*dy)

      Diff_x = (Afn(1:nx-1, 0:ny-2) + Afn(0:nx-2, 0:ny-2))/2._wp

      Diff_y = (Afn(0:nx-2, 1:ny-1) + Afn(0:nx-2, 0:ny-2))/2._wp


      !Convert v and Diff into the coefficients for the ADI method via the Il'in scheme
      call Ilin_Coefficients(vi_x, Diff_x, vi_y, Diff_y, cxp, cxm, cyp, cym)

      return
   end subroutine h_coefficients


   ! ------------------------------------------------------------------------
   subroutine h_res_update(h_res,h_cur,h_prev,dt_val)
      !This subroutine updates the thickness of the residually trapped region, h_res.
      !The rate of change depends on whether h is increasing or decreasing, and whether
      !h_res is zero or positive.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) ::  h_res
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) ::  h_cur, h_prev
      real(wp), dimension(0:nx-1, 0:ny-1) ::  C_dissolve, R_trap_array, h_inc_array, h_dec_array
      real(wp), intent(in) :: dt_val
      real(wp) :: sigma

      sigma = C_sat*(1._wp-s_a_i-s_c_r)/( s_c_r + C_sat*(1._wp-s_c_r))
      C_dissolve = 1._wp / (poro_h*( s_c_r + C_sat*(1._wp-s_c_r) ))

      !Locations where h is increasing, h is decreasing,
      !and where residually trapped CO2 is present
      h_inc_array  = merge(1._wp, 0._wp, h_cur >= h_prev)
      h_dec_array  = merge(1._wp, 0._wp, h_cur <  h_prev)
      R_trap_array = merge(1._wp, 0._wp, h_res >  0._wp )


      h_res = h_res + &
      &        (h_cur-h_prev)*( -h_inc_array*R_trap_array + (sigma-1._wp)*h_dec_array ) + &
      &        ( - C_dissolve * q_dissolve * dt_val * R_trap_array )

      ! Impose that h_res = 0 at the boundaries.
      ! ** This is current a fix to avoid occasional spurious behaviour at the edges! **
      h_res(0,:)     = 0._wp
      h_res(nx-1,:)  = 0._wp
      h_res(:,0)     = 0._wp
      h_res(:,ny-1)  = 0._wp

      !Make sure that h_res >= 0
      h_res = merge( h_res , 0._wp , h_res > 0._wp )

   end subroutine h_res_update


   ! ------------------------------------------------------------------------
   subroutine Ilin_Coefficients(vi_x, D_x, vi_y, D_y, cxp, cxm, cyp, cym)
      !This subroutine takes the advective velocity (vi_x,vi_y) and the diffusive term Diff
      !for and Advection-Diffusion equation, and applies the Il'in scheme to determine
      !the appropriate coefficients for the corresponding matrix problem so that
      !advection and diffusion are appropriately captured.

      implicit none
      real(wp), dimension(0:nx-2, 0:ny-2), intent(in) :: vi_x, D_x, vi_y, D_y
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  cxp, cxm, cyp, cym
      real(wp), dimension(0:nx-2, 0:ny-2) :: qi_x, ai_x, qi_y, ai_y, Diff_x, Diff_y

      !Dummy variables, as we may change the values (to make sure they're >=0)
      Diff_x = D_x
      Diff_y = D_y

      ! -- x coefficients ----------------------------------------------------------------------

      where (Diff_x <= 0._wp)
         !Purely advective here
         ai_x = sign(1._wp, vi_x) !This applies the sign of vi_x(i) to 1._wp
         Diff_x = 0._wp
      elsewhere
         ! Il'in upwinding parameters q, alpha
         qi_x = (vi_x*dx)/(2._wp*Diff_x)

         where (dabs(qi_x) .le. 1e-1_wp)
            ai_x = (1._wp/3._wp)*qi_x + (-1._wp/45._wp)*qi_x**3 + (2._wp/945._wp)*qi_x**5
         elsewhere
            ai_x = 1._wp/dtanh(qi_x) - 1._wp/qi_x
         end where

      end where

      cxp(0:nx-2, 0:ny-2) = Diff_x(0:nx-2, 0:ny-2)/dx &
      &                    + (vi_x(0:nx-2, 0:ny-2)/2._wp)*(1._wp + ai_x(0:nx-2, 0:ny-2))
      cxm(0:nx-2, 0:ny-2) = Diff_x(0:nx-2, 0:ny-2)/dx &
      &                    - (vi_x(0:nx-2, 0:ny-2)/2._wp)*(1._wp - ai_x(0:nx-2, 0:ny-2))

      cxp(nx-1, :) = 0._wp
      cxm(nx-1, :) = 0._wp

      ! -- y coefficients ----------------------------------------------------------------------

      where (Diff_y <= 0._wp)
         !Purely advective here
         ai_y = sign(1._wp, vi_y) !This applies the sign of vi_y(i) to 1._wp
         Diff_y = 0._wp
      elsewhere
         ! Il'in upwinding parameters q, alpha
         qi_y = (vi_y*dy)/(2._wp*Diff_y)

         where (dabs(qi_y) .le. 1e-1_wp)
            ai_y = (1._wp/3._wp)*qi_y + (-1._wp/45._wp)*qi_y**3 + (2._wp/945._wp)*qi_y**5
         elsewhere
            ai_y = 1._wp/dtanh(qi_y) - 1._wp/qi_y
         end where

      end where

      cyp(0:nx-2, 0:ny-2) = Diff_y(0:nx-2, 0:ny-2)/dy &
      &                    + (vi_y(0:nx-2, 0:ny-2)/2._wp)*(1._wp + ai_y(0:nx-2, 0:ny-2))
      cym(0:nx-2, 0:ny-2) = Diff_y(0:nx-2, 0:ny-2)/dy &
      &                    - (vi_y(0:nx-2, 0:ny-2)/2._wp)*(1._wp - ai_y(0:nx-2, 0:ny-2))

      cyp(:, ny-1) = 0._wp
      cym(:, ny-1) = 0._wp

      return

   end subroutine Ilin_Coefficients


   ! ------------------------------------------------------------------------
   subroutine ADI_x_solve(F_new, F_old, cxp, cxm, cyp, cym, ct, Forcing_term, dtb)
      !This subroutine formulates the ADI update stage in x for given coefficents.
      !The resulting matrix problem is tridiagonal; gtri is then called to solve this.
      !This is done iteratively, solving over x for each y value in the domain.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  F_old, cxp, cxm, cyp, cym, ct, Forcing_term
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  F_new
      real(wp), dimension(0:nx-1) :: ab, al, ad, au, ax
      real(wp), intent(in) :: dtb
      integer ::  j

      !Apply the appropriate BCs to h based on the values specified in h_BC_params

      !West side
      if (h_BC_params(3) == 1) then
         ! Dirichlet Boundary conditions
         ab(0) = 0._wp
         al(0) = 0._wp
         ad(0) = 1._wp
         au(0) = 0._wp
      else if (h_BC_params(3) == 2) then
         ! Neumann Boundary conditions - zero flux
         ! ab(0) = - Gamma_val * (H0(1,0) - H0(0,0))
         ab(0) = 0._wp
         al(0) = 0._wp
         ad(0) = -1._wp
         au(0) = 1._wp
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if

      !East side
      if (h_BC_params(1) == 1) then
         ! Dirichlet Boundary conditions
         ab(nx-1) = 0._wp
         al(nx-1) = 0._wp
         ad(nx-1) = 1._wp
         au(nx-1) = 0._wp
      else if (h_BC_params(1) == 2) then
         ! Neumann Boundary conditions
         ab(nx-1) = 0._wp
         al(nx-1) = -1._wp
         ad(nx-1) = 1._wp
         au(nx-1) = 0._wp
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if


      !Formulate the lower diagonal (al), main diagonal (ad), upper diagonal (au), and the
      !right-hand side (ab) and solve the corresponding tridiagonal matrix problem,
      !iterating over the index for y
      do j = 1, ny-2

         ab(1:nx-2) = ct(1:nx-2, j)*F_old(1:nx-2, j) &
         &           + (dtb/dy)*cym(1:nx-2,j)*F_old(1:nx-2,j+1) &
         &           - (dtb/dy)*(cyp(1:nx-2,j)+cym(1:nx-2,j-1))*F_old(1:nx-2,j) &
         &           + (dtb/dy)*cyp(1:nx-2,j-1)*F_old(1:nx-2,j-1) &
         &           + dtb*Forcing_term(1:nx-2, j)

         al(1:nx-2) = -(dtb/dx)*cxp(0:nx-3, j)
         ad(1:nx-2) = ct(1:nx-2, j) + (dtb/dx)*(cxp(1:nx-2, j) + cxm(0:nx-3, j))
         au(1:nx-2) = -(dtb/dx)*cxm(1:nx-2, j)

         ! tridiagonal solver
         call gtri(al, ad, au, ab, ax, 0, nx-1)
         F_new(:, j) = ax

      end do

      return
   end subroutine ADI_x_solve


   ! ------------------------------------------------------------------------
   subroutine ADI_y_solve(F_new, F_old, cxp, cxm, cyp, cym, ct, Forcing_term, dtb)
      !This subroutine formulates the ADI update stage in y for given coefficents.
      !The resulting matrix problem is tridiagonal; gtri is then called to solve this.
      !This is done iteratively, solving over y for each x value in the domain.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  F_old, cxp, cxm, cyp, cym, ct, Forcing_term
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  F_new
      real(wp), dimension(0:ny-1) :: ab, al, ad, au, ax
      real(wp), intent(in) :: dtb
      integer ::  i

      !Apply the appropriate BCs to h based on the values specified in h_BC_params

      !North side
      if (h_BC_params(0) == 1) then
         ! Dirichlet Boundary conditions
         ab(0) = 0._wp
         al(0) = 0._wp
         ad(0) = 1._wp
         au(0) = 0._wp
      else if (h_BC_params(0) == 2) then
         ! Neumann Boundary conditions
         ab(0) = 0._wp
         al(0) = 0._wp
         ad(0) = -1._wp
         au(0) = 1._wp
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if

      !South side
      if (h_BC_params(2) == 1) then
         ! Dirichlet Boundary conditions
         ab(ny-1) = 0._wp
         al(ny-1) = 0._wp
         ad(ny-1) = 1._wp
         au(ny-1) = 0._wp
      else if (h_BC_params(2) == 2) then
         ! Neumann Boundary conditions
         ab(ny-1) = 0._wp
         al(ny-1) = -1._wp
         ad(ny-1) = 1._wp
         au(ny-1) = 0._wp
      else
         print '("!!! -- Incorrect BC Flag! -- !!!")'
         stop
      end if


      !Formulate the lower diagonal (al), main diagonal (ad), upper diagonal (au), and the
      !right-hand side (ab) and solve the corresponding tridiagonal matrix problem,
      !iterating over the index for x
      do i = 1, nx-2

         ab(1:ny-2) = ct(i, 1:ny-2)*F_old(i, 1:ny-2) &
         &           + (dtb/dx)*cxm(i,1:ny-2)*F_old(i+1,1:ny-2) &
         &           - (dtb/dx)*(cxp(i,1:ny-2)+cxm(i-1,1:ny-2))*F_old(i,1:ny-2) &
         &           + (dtb/dx)*cxp(i-1,1:ny-2)*F_old(i-1,1:ny-2) &
         &           + dtb*Forcing_term(i, 1:ny-2)

         al(1:ny-2) = -(dtb/dy)*cyp(i, 0:ny-3)
         ad(1:ny-2) = ct(i, 1:ny-2) + (dtb/dy)*(cyp(i, 1:ny-2) + cym(i, 0:ny-3))
         au(1:ny-2) = -(dtb/dy)*cym(i, 1:ny-2)

         ! tridiagonal solver
         call gtri(al, ad, au, ab, ax, 0, ny-1)
         F_new(i, :) = ax

      end do

      return
   end subroutine ADI_y_solve


   ! ------------------------------------------------------------------------
   subroutine gtri(l, d, u, b, x, n1, n2)
      ! This subroutine solves the tridiagonal system LnXn-1 + DnXn + UnXn+1=Bn
      !via the Thomas algorithm

      implicit none
      integer, intent(in) :: n1, n2
      real(wp), dimension(n1:n2), intent(in) :: l, d, u, b
      real(wp), dimension(n1:n2), intent(out) :: x
      real(wp), dimension(n1:n2) :: wb
      real(wp) :: s
      integer :: i

      x(n1)  = d(n1)
      wb(n1) = b(n1)

      do i = n1 + 1, n2
         s = l(i)/x(i - 1)
         x(i) = d(i) - u(i - 1)*s
         wb(i) = b(i) - wb(i - 1)*s
      end do

      x(n2) = wb(n2)/x(n2)

      do i = n2 - 1, n1, -1
         x(i) = (wb(i) - x(i + 1)*u(i))/x(i)
      end do

   end subroutine gtri


   ! ------------------------------------------------------------------------
   function err(Fa, Fb)
      !This function calculates the L2 error between two arrays. 
      !This is normalised (if possible) by the maximum value of 
      !one of the two arrays to give a measure of the relative error. 

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: Fa, Fb
      real(wp) :: err, scale

      !Determine the scale to normalise the error by. If one of the two
      !arrays is sufficiently close to 0, then set the scale to 1 - then 
      !the error value will only be small if the other array is also 
      !close to zero
      scale = maxval(abs(Fb))
      if (abs(scale) .le. 1e-8_wp) then
         scale = 1._wp
      end if

      err = sum((Fa - Fb)**2.)/scale
      err = dsqrt(err/(nx*ny))

   end function err


   ! ------------------------------------------------------------------------
   subroutine save_iter_output(h_array, h_res_array, P_array, h_cur, h_res_cur, P_cur, &
   & t, tplot, plot_count, target_plot_times, plot_times, V_injected)
      !This subroutine calculates the the active and trapped volumes of CO2 after each
      !successful timestep, and saves them to an output file.
      !If the timestepping routine has also passed a planned output time, the calculated 
      !profiles h, h_res, and P are saved into their respective output arrays.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1, 0:np-1), intent(inout) :: h_array, h_res_array, P_array
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: h_cur, h_res_cur, P_cur
      real(wp), dimension(0:np-1), intent(inout) :: plot_times
      real(wp), dimension(0:np-1), intent(in) :: target_plot_times
      real(wp), intent(inout) :: tplot
      real(wp), intent(in) :: t, V_injected
      integer, intent(inout) :: plot_count
      real(wp) :: V_active, V_trapped

      ! calculate and record the volumes of Active CO2, Trapped CO2, and total injected CO2
      V_active    = sum( poro_h(:, :)*h_cur(:, :) )*(1._wp - s_a_i)*dx*dy
      V_trapped   = sum( poro_h(:, :)*(h_res_cur(:,:)) )*s_c_r*dx*dy
      write (99, *) t, V_active, V_trapped, V_injected


      !Check if we've reached a plot time
      if (plot_count .le. np-1 .and. t >= tplot) then

         !Save data for current plot time
         print '("[Save #", i4 ,"] t = ", f15.5 ,"")', plot_count, t

         h_array(:,:,plot_count)       = h_cur(:,:)
         h_res_array(:,:,plot_count)   = h_res_cur(:,:)
         P_array(:,:,plot_count)       = P_cur(:,:)
         plot_times(plot_count)        = t


         if (plot_count .lt. np-1) then
            !Update next target plot time
            tplot = target_plot_times(plot_count+1)
         end if

         !Update plot number
         plot_count = plot_count + 1

      end if

   end subroutine save_iter_output

   ! ------------------------------------------------------------------------
   subroutine generate_flux_array(Q_array,flux_idx_cur,t_val)
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
      real(wp), dimension(0:n_inj_locs-1) :: Q_vals
      integer :: k

      Q_array = 0._wp

      if (t_val < t_flux_vals(1)) then
         !Injection hasn't started yet
         flux_idx_cur = 0
         Q_vals(:) = 0._wp

      else if (t_val >= t_flux_vals(n_flux_times-1)) then
         !After the last change in flux - use the index of the last specified value
         flux_idx_cur = n_flux_times-1

      else
         !Calculate the index of the most recent specified flux time 
         !in t_flux_vals before the current time t.
         !Start at the previous index value to save excessive counting.
         k = flux_idx_cur
         do while ( t_flux_vals(k) <= t_val )
            k = k+1
         end do
         flux_idx_cur = k-1

      end if


      !Use the previous flux value - piecewise-constant interpolation
      !This reads a flux value for each injection location
      Q_vals = Q_flux_vals(flux_idx_cur,:)

      !Apply these flux values at their respective injection locations
      do k = 0, n_inj_locs-1
         Q_array( Q_inj_locs(k,0) , Q_inj_locs(k,1) ) = Q_vals(k)
      end do



   end subroutine generate_flux_array

end module CO2GraVISim_solver
