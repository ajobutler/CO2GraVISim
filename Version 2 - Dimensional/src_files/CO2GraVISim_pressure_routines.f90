module CO2GraVISim_pressure_routines
    ! This contains the Pressure subroutines for CO2GraVISim, used to solve for 
    ! the non-hydrostatic component of the ambient pressure P

    use CO2GraVISim_global !wp, iter_max, tmax, dt_init, dt_min, errmax_t_step, errmax_P_iter, P_iter_omega, t0
    use CO2GraVISim_input_parameters !H0, B0, D0, perm_h, poro_h
    ! use CO2GraVISim_injection_profiles !n_inj_locs, n_flux_times, Q_inj_locs, Q_flux_vals, t_flux_vals
    use CO2GraVISim_vertical_structure !poro_cur
   !  use CO2GraVISim_iter_calculations
    use CO2GraVISim_err_calc
    use CO2GraVISim_timer
   !  use OMP_LIB

    implicit none

    contains


   subroutine Pressure_calculation(P_new, P_old, h, h_res, Q, P_iter_count, Confined_flag)
      ! This is a wrapper subroutine that chooses the behaviour based on whether the main code
      ! is being run in Confined or Unconfined mode. If Unconfined, then the pressure calculations 
      ! can be skipped entirely, setting P=0. Otherwise, the pressure calculations need to be performed.

      implicit none
      real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  P_old, h, h_res, Q
      logical, intent(in) :: Confined_flag
      integer, intent(out) :: P_iter_count
      real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  P_new

      if (Confined_flag) then
         !Confined case
         call Pressure_solve(P_new, P_old, h, h_res, Q, P_iter_count)
      else
         !Unconfined case
         P_new = 0._wp
         P_iter_count = 0
      end if


   end subroutine Pressure_calculation

    ! ------------------------------------------------------------------------
   subroutine Pressure_solve(P_new, P_old, h_cur, h_res_cur, Q, iter_count)
    ! This subroutine calculates the coefficients that go into the Poisson equation
    ! for the ambient pressure, which is then solved by Pressure_SOR.
    ! This equation is of the form
    ! Div( A * grad(P) ) = - Div( B * grad(G) ) - Q
    ! This is ultimately discretised as
    !     cxm(i,j)*P(i-1,j  ) + cx0(i,j)*P(i,j) + cxp(i,j)*P(i+1,j  )
    !   + cym(i,j)*P(i  ,j-1) + cy0(i,j)*P(i,j) + cyp(i,j)*P(i  ,j+1)
    !   = Forcing_term(i,j)

    implicit none
    real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  P_old, h_cur, h_res_cur, Q
    real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  P_new
    real(wp), dimension(0:nx-1, 0:ny-1) ::  Axavg, Cxp, Cxm, Cx0
    real(wp), dimension(0:nx-1, 0:ny-1) ::  Ayavg, Cyp, Cym, Cy0
    real(wp), dimension(0:nx-1, 0:ny-1) :: Afn, Bfn, Gfn, Forcing_term
    real(wp), dimension(0:nx-1, 0:ny-1) :: kappa_c, kappa_a
    integer, intent(out) :: iter_count
    !real(wp) :: start_time, end_time

    ! start_time = omp_get_wtime()

    ! Initialise arrays
    Axavg = 0._wp
    Ayavg = 0._wp
    Cxp = 0._wp
    Cxm = 0._wp
    Cx0 = 0._wp
    Cyp = 0._wp
    Cym = 0._wp
    Cy0 = 0._wp
    Forcing_term = 0._wp


    ! Calculate permeability averages across the mobile CO2 (kappa_c) and the rest of the column (kappa_a)  
    call get_permeability_average_mobile(h_cur ,kappa_c)
    call get_permeability_average_ambient(h_cur,h_res_cur,kappa_a)


    ! Build relevant terms from the governing pressure equation
    Afn = h_cur*kappa_c + M*(D0 - h_cur)*kappa_a
    Bfn = Gamma_val*h_cur*kappa_c
    Gfn = H0 + h_cur


    ! Calculate the values of A at the interfaces between the grid points
    Axavg(0:nx-2, 0:ny-1) = (Afn(1:nx-1, 0:ny-1) + Afn(0:nx-2, 0:ny-1))/(2._wp)
    Ayavg(0:nx-1, 0:ny-2) = (Afn(0:nx-1, 1:ny-1) + Afn(0:nx-1, 0:ny-2))/(2._wp)


    ! x coefficients
    Cxp(1:nx-2, 1:ny-2) = Axavg(1:nx-2, 1:ny-2)/(dx**2)
    Cxm(1:nx-2, 1:ny-2) = Axavg(0:nx-3, 1:ny-2)/(dx**2)
    Cx0(1:nx-2, 1:ny-2) =  -Cxp(1:nx-2, 1:ny-2) - Cxm(1:nx-2, 1:ny-2)

    ! y coefficients
    Cyp(1:nx-2, 1:ny-2) = Ayavg(1:nx-2, 1:ny-2)/(dy**2)
    Cym(1:nx-2, 1:ny-2) = Ayavg(1:nx-2, 0:ny-3)/(dy**2)
    Cy0(1:nx-2, 1:ny-2) =  -Cyp(1:nx-2, 1:ny-2) - Cym(1:nx-2, 1:ny-2)


    ! Build the forcing term up from B, G, and Q
    Forcing_term(1:nx - 2, 1:ny - 2) = Q(1:nx-2, 1:ny-2)/(dx*dy) &
    & + (Bfn(2:nx-1, 1:ny-2) + Bfn(1:nx-2, 1:ny-2))*(Gfn(2:nx-1, 1:ny-2) - Gfn(1:nx-2, 1:ny-2))/(2._wp*(dx**2)) &
    & - (Bfn(1:nx-2, 1:ny-2) + Bfn(0:nx-3, 1:ny-2))*(Gfn(1:nx-2, 1:ny-2) - Gfn(0:nx-3, 1:ny-2))/(2._wp*(dx**2)) &
    & + (Bfn(1:nx-2, 2:ny-1) + Bfn(1:nx-2, 1:ny-2))*(Gfn(1:nx-2, 2:ny-1) - Gfn(1:nx-2, 1:ny-2))/(2._wp*(dy**2)) &
    & - (Bfn(1:nx-2, 1:ny-2) + Bfn(1:nx-2, 0:ny-3))*(Gfn(1:nx-2, 1:ny-2) - Gfn(1:nx-2, 0:ny-3))/(2._wp*(dy**2))



    ! Solve for the new pressure
    call Pressure_SOR(P_new, P_old, cxp, cxm, cx0, cyp, cym, cy0, Forcing_term, iter_count)



    ! end_time = omp_get_wtime()
    ! ! Add runtime to running total
    ! P_solve_time = P_solve_time + (end_time - start_time)
    ! ! Increase count of calls
    ! P_solve_calls = P_solve_calls + 1

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
    integer, intent(out) :: iter_count
    real(wp), dimension(0:nx-1, 0:ny-1) ::  P_prev, dxy0, dxm, dxp, dym, dyp!, P_temp
    real(wp), dimension(0:nx-1, 0:ny-1) :: Inverse_Factor, Scaled_Forcing
    real(wp), dimension(0:nx-1, 0:ny-1) :: Axm, Axp, Aym, Ayp
    real(wp) :: R, error_val, sigma
    integer ::  i, j
    !real(wp) :: start_time, end_time

    ! start_time = omp_get_wtime()

    ! !Arrays to vectorise inner loops in parallel section
    ! !The sizes here use integer division to take the appropriate floor value
    ! real(wp), dimension( 0:((nx-1)/2 - 1) ) :: V_1 !Odd indices between 1 and (nx-2) (then 0-indexed)
    ! real(wp), dimension( 0: ( (nx-2) - (nx-1)/2 - 1 ) ) :: V_2 !Even indices between 2 and (nx-2) (then 0-indexed)

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

    sigma = (1._wp - P_iter_omega)

    Inverse_Factor(1:nx-2, 1:ny-2) = (P_iter_omega/dxy0(1:nx-2, 1:ny-2))
    Axm(1:nx-2, 1:ny-2) = Inverse_Factor(1:nx-2, 1:ny-2) * dxm(1:nx-2, 1:ny-2)
    Axp(1:nx-2, 1:ny-2) = Inverse_Factor(1:nx-2, 1:ny-2) * dxp(1:nx-2, 1:ny-2)
    Aym(1:nx-2, 1:ny-2) = Inverse_Factor(1:nx-2, 1:ny-2) * dym(1:nx-2, 1:ny-2)
    Ayp(1:nx-2, 1:ny-2) = Inverse_Factor(1:nx-2, 1:ny-2) * dyp(1:nx-2, 1:ny-2)
    Scaled_Forcing(1:nx-2, 1:ny-2) = Inverse_Factor(1:nx-2, 1:ny-2) * Forcing_term(1:nx-2, 1:ny-2)

    ! ! !$OMP PARALLEL SHARED(P_new,P_prev,Axm,Axp,Aym,Ayp,P_temp,Scaled_Forcing,error_val) PRIVATE(j,V_1,V_2) NUM_THREADS(4)

    ! ! do while (error_val > errmax_P_iter)

    ! !    !$OMP SINGLE
    ! !    P_temp = P_new
    ! !    !$OMP END SINGLE


    ! !    !Black squares, even rows
    ! !    !$OMP DO SCHEDULE(STATIC)
    ! !    do j = 1, ny-2, 2
    ! !       !i = 1, nx-2, 2

    ! !       V_1(:) =               sigma*P_prev( 1:nx-2:2 , j     ) &
    ! !          &   - Axm( 1:nx-2:2 , j )*P_temp( 0:nx-3:2 , j     ) &
    ! !          &   - Axp( 1:nx-2:2 , j )*P_temp( 2:nx-1:2 , j     ) &
    ! !          &   - Aym( 1:nx-2:2 , j )*P_temp( 1:nx-2:2 , j - 1 ) &
    ! !          &   - Ayp( 1:nx-2:2 , j )*P_temp( 1:nx-2:2 , j + 1 ) &
    ! !          &   - Scaled_Forcing( 1:nx-2:2 , j )

    ! !       P_new( 1:nx-2:2 , j ) = V_1(:)

    ! !    end do
    ! !    !$OMP END DO NOWAIT


    ! !    !Black squares, odd rows
    ! !    !$OMP DO SCHEDULE(STATIC)
    ! !    do j = 2, ny-2, 2
    ! !       !i = 2, nx-2, 2

    ! !       V_2(:) =                sigma*P_prev( 2:nx-2:2 , j     ) &
    ! !          &   - Axm( 2:nx-2:2 , j )*P_temp( 1:nx-3:2 , j     ) &
    ! !          &   - Axp( 2:nx-2:2 , j )*P_temp( 3:nx-1:2 , j     ) &
    ! !          &   - Aym( 2:nx-2:2 , j )*P_temp( 2:nx-2:2 , j - 1 ) &
    ! !          &   - Ayp( 2:nx-2:2 , j )*P_temp( 2:nx-2:2 , j + 1 ) &
    ! !          &   - Scaled_Forcing( 2:nx-2:2 , j )

    ! !       P_new( 2:nx-2:2 , j ) = V_2(:)

    ! !    end do
    ! !    !$OMP END DO

    ! !    !$OMP BARRIER

    ! !    !$OMP SINGLE
    ! !    P_temp = P_new
    ! !    !$OMP END SINGLE


    ! !    !Red squares, even rows
    ! !    !$OMP DO SCHEDULE(STATIC)
    ! !    do j = 2, ny-2, 2
    ! !       !i = 1, nx-2, 2

    ! !       V_1(:) =               sigma*P_prev( 1:nx-2:2 , j     ) &
    ! !          &   - Axm( 1:nx-2:2 , j)*P_temp( 0:nx-3:2 , j     ) &
    ! !          &   - Axp( 1:nx-2:2 , j)*P_temp( 2:nx-1:2 , j     ) &
    ! !          &   - Aym( 1:nx-2:2 , j)*P_temp( 1:nx-2:2 , j - 1 ) &
    ! !          &   - Ayp( 1:nx-2:2 , j)*P_temp( 1:nx-2:2 , j + 1 ) &
    ! !          &   - Scaled_Forcing( 1:nx-2:2 , j )

    ! !       P_new( 1:nx-2:2 , j ) = V_1(:)

    ! !    end do
    ! !    !$OMP END DO NOWAIT

    ! !    !Red squares, odd rows
    ! !    !$OMP DO SCHEDULE(STATIC)
    ! !    do j = 1, ny-2, 2
    ! !       !i = 2, nx-2, 2

    ! !       V_2(:) =                sigma*P_prev( 2:nx-2:2 , j     ) &
    ! !          &   - Axm( 2:nx-2:2 , j )*P_temp( 1:nx-3:2 , j     ) &
    ! !          &   - Axp( 2:nx-2:2 , j )*P_temp( 3:nx-1:2 , j     ) &
    ! !          &   - Aym( 2:nx-2:2 , j )*P_temp( 2:nx-2:2 , j - 1 ) &
    ! !          &   - Ayp( 2:nx-2:2 , j )*P_temp( 2:nx-2:2 , j + 1 ) &
    ! !          &   - Scaled_Forcing( 2:nx-2:2 , j )

    ! !       P_new( 2:nx-2:2 , j ) = V_2(:)

    ! !    end do
    ! !    !$OMP END DO

    ! !    !$OMP BARRIER

    ! !    !Apply the BCs to the edges of P_new
    ! !    !$OMP SINGLE
    ! !    call Pressure_BCs(P_new)

    ! !    if ( mod(iter_count, P_iter_check) == 0 ) then
    ! !       error_val = err(P_prev, P_new, nz, ny)
    ! !    end if

    ! !    P_prev = P_new
    ! !    iter_count = iter_count + 1
    ! !    !$OMP END SINGLE


    ! ! end do

    ! ! !$OMP END PARALLEL

    do while (error_val > errmax_P_iter)
       do j = 1, ny-2
          do i = 1, nx-2

             R =   - Axm(i, j)*P_new(i - 1, j    )
             R = R - Axp(i, j)*P_new(i + 1, j    )
             R = R - Aym(i, j)*P_new(i    , j - 1)
             R = R - Ayp(i, j)*P_new(i    , j + 1)
             R = R - Scaled_Forcing(i, j)
             R = R + (1._wp - P_iter_omega)*P_prev(i, j)

             P_new(i, j) = R

             ! P_new(i, j) = (1._wp - P_iter_omega)*P_prev(i, j) + &
             ! & (P_iter_omega/dxy0(i, j))* &
             ! &            (- dxm(i, j)*P_new(i - 1, j    ) - dxp(i, j)*P_new(i + 1, j    ) &
             ! &             - dym(i, j)*P_new(i    , j - 1) - dyp(i, j)*P_new(i    , j + 1) &
             ! &         - Forcing_term(i, j))

          end do
       end do

       !!!
       !! TODO: I need to add a way of dealing with NaNs appearing in the pressure
       !!!

       !Apply the BCs to the edges of P_new
       call Pressure_BCs(P_new)

       if ( mod(iter_count, P_iter_check) == 0 ) then
          error_val = err(P_prev, P_new, nx, ny)
       end if

       P_prev = P_new
       iter_count = iter_count + 1

    end do
  
    if ( count(isnan(P_new)) > 0 ) then
       write(*,'(I0,A,I0)') count(isnan(P_new)), ' NaNs at P_iter = ', iter_count     
   end if
   
   
 !   end_time = omp_get_wtime()
 !   ! Add runtime to running total
 !   P_SOR_time = P_SOR_time + (end_time - start_time)
 !   ! Increase count of calls
 !   P_SOR_calls = P_SOR_calls + 1

    return

 end subroutine Pressure_SOR


 ! ------------------------------------------------------------------------
 subroutine Pressure_BCs(P)
    !This subroutine imposed the boundary conditions -- either Dirichlet or
    !Neumann -- on the edges of the Pressure solution array

    implicit none
    real(wp), dimension(0:nx-1,0:ny-1), intent(inout):: P
    !real(wp) :: start_time, end_time

    ! start_time = omp_get_wtime()

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
  
    ! end_time = omp_get_wtime()
    ! ! Add runtime to running total
    ! P_BCs_time = P_BCs_time + (end_time - start_time)
    ! ! Increase count of calls
    ! P_BCs_calls = P_BCs_calls + 1


 end subroutine Pressure_BCs


end module CO2GraVISim_pressure_routines