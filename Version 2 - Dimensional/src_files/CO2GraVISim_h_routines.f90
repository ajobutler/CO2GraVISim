module CO2GraVISim_h_routines

    use CO2GraVISim_global
    use CO2GraVISim_input_parameters
    use CO2GraVISim_vertical_structure
    use CO2GraVISim_timer

    implicit none
    
    contains

    ! ------------------------------------------------------------------------
    subroutine h_coefficients(h_cur, P_cur, h_res, dh_dt, sigma_a_V, cxp, cxm, cyp, cym, ct)
        !This subroutine calculates the coefficients to be used as part of the ADI method
        !to solve for the current thickness h.
        !The equation for h is of the form
        !            ct*dh/dt + Div( v * h ) = Div( Diff * grad(h) ) + Q
        !Here the discretised forms of v and Diff are calculated and passed to
        !Ilin_Coefficients, which calculates the corresponding coefficients for the ADI method.

        implicit none
        real(wp), dimension(0:nx-1, 0:ny-1), intent(in)  ::  h_cur, P_cur, h_res, dh_dt, sigma_a_V
        real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  cxp, cxm, cyp, cym, ct
        real(wp), dimension(0:nx-2, 0:ny-2) :: Diff_x, vi_x, Diff_y, vi_y
        real(wp), dimension(0:nx-1, 0:ny-1) :: Afn, Bfn !, Trapping_array, h_inc_array, h_dec_array, No_trapping_array
        real(wp), dimension(0:nx-1, 0:ny-1) :: poro_cur, kappa_c
        !real(wp) :: start_time, end_time
    
        ! start_time = omp_get_wtime()

        ! -- Current Porosity and Permeability --------------------------------------------------------------------------

        call get_porosity_value(h_cur,poro_cur)
        call get_permeability_average_mobile(h_cur,kappa_c)

        ! -- Calculate ct -----------------------------------------------------------------------------------------------
        !! Time derivative prefactor - depends on sign of dh/dt, and where h_res>0

        ct = 0._wp


        where ( (h_res <= 0._wp) .and. (dh_dt >= 0._wp) )

            ct = poro_cur * ( (1._wp-s_a_i) + C_sat*s_a_i )

        elsewhere

            ct = poro_cur * (1._wp - sigma_a_V*C_sat) * (1._wp - s_a_i - s_c_r)

        end where




        ! where ( (h_res <= 0._wp) .and. (dh_dt >= 0._wp ))
        
        !    ct = poro_cur * ( (1._wp-s_a_i) + C_sat*s_a_i )

        ! elsewhere ( (h_res <= 0._wp) .and. (dh_dt < 0._wp ))

        !    ct = poro_cur * (1._wp - s_a_i - s_c_r) * &
        !    & ( (1._wp - 2._wp*C_sat) + C_sat*( s_c_r + C_sat*(1._wp-s_c_r))/s_c_r )             

        ! elsewhere ( (h_res > 0._wp) .and. sigma_a_V < 0._wp )

        !    ct = poro_cur * (1._wp - s_a_i - s_c_r) * &
        !    & ( (1._wp - 2._wp*C_sat) + C_sat*( s_c_r + C_sat*(1._wp-s_c_r))/s_c_r )

        ! elsewhere ( (h_res > 0._wp) .and. sigma_a_V >= 0._wp )

        !    ct = poro_cur * ( 1._wp - s_a_i - s_c_r )
        
        ! end where




        ! ! where (dh_dt < 0._wp )

        ! !     ct = poro_cur * ( 1._wp - s_a_i - s_c_r )

        ! !     elsewhere

        ! !     where (h_res > 0._wp)

        ! !         ct = poro_cur * (1._wp - C_sat) * ( 1._wp - s_a_i - s_c_r )

        ! !     elsewhere  
        ! !         ! dh_dt>=0 , h_res = 0
        ! !         ct = ( 1._wp - s_a_i +C_sat*s_a_i )

        ! !     end where

        ! ! end where



        ! -- Calculate v and Diff ---------------------------------------------------------------------------------------
        Afn = Gamma_val*h_cur*kappa_c
        Bfn = -(P_cur + Gamma_val*H0)

        vi_x = (kappa_c(1:nx-1, 0:ny-2) + kappa_c(0:nx-2, 0:ny-2))*(Bfn(1:nx-1, 0:ny-2) - Bfn(0:nx-2, 0:ny-2))/(2._wp*dx)

        vi_y = (kappa_c(0:nx-2, 1:ny-1) + kappa_c(0:nx-2, 0:ny-2))*(Bfn(0:nx-2, 1:ny-1) - Bfn(0:nx-2, 0:ny-2))/(2._wp*dy)

        Diff_x = (Afn(1:nx-1, 0:ny-2) + Afn(0:nx-2, 0:ny-2))/2._wp

        Diff_y = (Afn(0:nx-2, 1:ny-1) + Afn(0:nx-2, 0:ny-2))/2._wp


        !Convert v and Diff into the coefficients for the ADI method via the Il'in scheme
        call Ilin_Coefficients(vi_x, Diff_x, vi_y, Diff_y, cxp, cxm, cyp, cym)
    
    
        ! end_time = omp_get_wtime()
        ! ! Add runtime to running total
        ! h_coeff_time = h_coeff_time + (end_time - start_time)
        ! ! Increase count of calls
        ! h_coeff_calls = h_coeff_calls + 1

        return
    end subroutine h_coefficients


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
        !real(wp) :: start_time, end_time
    
        ! start_time = omp_get_wtime()
    
        !Initialise arrays
        cxp     = 0._wp
        cxm     = 0._wp
        cyp     = 0._wp
        cym     = 0._wp
        qi_x    = 0._wp
        ai_x    = 0._wp
        qi_y    = 0._wp
        ai_y    = 0._wp

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

            where (abs(qi_x) .le. 1e-1_wp)
                ai_x = (1._wp/3._wp)*qi_x + (-1._wp/45._wp)*qi_x**3 + (2._wp/945._wp)*qi_x**5
            elsewhere
                ai_x = 1._wp/tanh(qi_x) - 1._wp/qi_x
            end where

        end where

        cxp(0:nx-2, 0:ny-2) = Diff_x(0:nx-2, 0:ny-2)/dx &
        &                    + (vi_x(0:nx-2, 0:ny-2)/2._wp)*(1._wp + ai_x(0:nx-2, 0:ny-2))
        cxm(0:nx-2, 0:ny-2) = Diff_x(0:nx-2, 0:ny-2)/dx &
        &                    - (vi_x(0:nx-2, 0:ny-2)/2._wp)*(1._wp - ai_x(0:nx-2, 0:ny-2))

        !!Already set when initialised
        !cxp(nx-1, :) = 0._wp
        !cxm(nx-1, :) = 0._wp

        ! -- y coefficients ----------------------------------------------------------------------

        where (Diff_y <= 0._wp)
            !Purely advective here
            ai_y = sign(1._wp, vi_y) !This applies the sign of vi_y(i) to 1._wp
            Diff_y = 0._wp
        elsewhere
            ! Il'in upwinding parameters q, alpha
            qi_y = (vi_y*dy)/(2._wp*Diff_y)

            where (abs(qi_y) .le. 1e-1_wp)
                ai_y = (1._wp/3._wp)*qi_y + (-1._wp/45._wp)*qi_y**3 + (2._wp/945._wp)*qi_y**5
            elsewhere
                ai_y = 1._wp/tanh(qi_y) - 1._wp/qi_y
            end where

        end where

        cyp(0:nx-2, 0:ny-2) = Diff_y(0:nx-2, 0:ny-2)/dy &
        &                    + (vi_y(0:nx-2, 0:ny-2)/2._wp)*(1._wp + ai_y(0:nx-2, 0:ny-2))
        cym(0:nx-2, 0:ny-2) = Diff_y(0:nx-2, 0:ny-2)/dy &
        &                    - (vi_y(0:nx-2, 0:ny-2)/2._wp)*(1._wp - ai_y(0:nx-2, 0:ny-2))

        !!Already set when initialised
        !cyp(:, ny-1) = 0._wp
        !cym(:, ny-1) = 0._wp
    
    
        ! end_time = omp_get_wtime()
        ! ! Add runtime to running total
        ! Ilin_Coeff_time = Ilin_Coeff_time + (end_time - start_time)
        ! ! Increase count of calls
        ! Ilin_Coeff_calls = Ilin_Coeff_calls + 1

        return

    end subroutine Ilin_Coefficients


    ! ------------------------------------------------------------------------
    subroutine calculate_h_forcing_term(Forcing_term_h, h, h_res, dh_dt, Q, sigma_a_V, F_a_H)

        implicit none
        real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: Forcing_term_h
        real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h, h_res, Q, dh_dt, sigma_a_V, F_a_H
        real(wp), dimension(0:nx-1,0:ny-1) :: Dissolution_array

        !Locations where convective dissolution is applied to the mobile CO2
        Dissolution_array = merge(1._wp, 0._wp, (h > 0._wp) .and. (h_res <= 0._wp) )

        Forcing_term_h = 0._wp

        ! Total flux term for h - combination of injection term and dissolution terms
        ! 1) Injection
        Forcing_term_h = Q/(dx*dy) 
        ! 2) Horizontal flux of dissolved CO2 in trapped region
        !     - Only contributes for sigma_a_V = sign(F_a_V) < 0 
        ! ! where (sigma_a_V < 0._wp)
        ! ! Forcing_term_h = Forcing_term_h +&
        ! ! &   - C_sat * ( -2._wp + (1._wp + C_sat*(1._wp-s_c_r)/s_c_r) ) * F_a_H
        ! ! endwhere
        where ( (h_res <= 0._wp) .and. (dh_dt >= 0._wp) )
            !Forcing_term_h is as above
        elsewhere
            Forcing_term_h = Forcing_term_h - C_sat * (1._wp - sigma_a_V) * F_a_H
        end where
        ! 3) Convective dissolution
        Forcing_term_h = Forcing_term_h - q_dissolve*Dissolution_array

        ! ! Forcing_term_h = Q/(dx*dy) - q_dissolve*Dissolution_array + (1._wp - sigma_a_V)*C_sat*F_a_H

    end subroutine calculate_h_forcing_term


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
        !real(wp) :: start_time, end_time
        integer ::  j
        
        ! start_time = omp_get_wtime()

        !Initialise arrays
        F_new   = 0._wp
        ab      = 0._wp
        al      = 0._wp
        ad      = 0._wp
        au      = 0._wp
        ax      = 0._wp

        !Apply the appropriate BCs to h based on the values specified in h_BC_params

        !East side
        select case (h_BC_params(1))
        case (1)
            ! Dirichlet Boundary conditions
            ab(nx-1) = 0._wp
            al(nx-1) = 0._wp
            ad(nx-1) = 1._wp
            au(nx-1) = 0._wp     
        case (2)
            ! Neumann Boundary conditions
            ab(nx-1) = 0._wp
            al(nx-1) = -1._wp
            ad(nx-1) = 1._wp
            au(nx-1) = 0._wp        
        case default 
            print '("!!! -- Incorrect BC Flag! -- !!!")'
            stop
        end select

        !West side
        select case (h_BC_params(3))
        case (1)
            ! Dirichlet Boundary conditions
            ab(0) = 0._wp
            al(0) = 0._wp
            ad(0) = 1._wp
            au(0) = 0._wp     
        case (2)
            ! Neumann Boundary conditions
            ab(0) = 0._wp !-(H0(0,ny-1)-H0(0,ny-2))
            al(0) = 0._wp
            ad(0) = -1._wp
            au(0) = 1._wp        
        case default 
            print '("!!! -- Incorrect BC Flag! -- !!!")'
            stop
        end select


        !Formulate the lower diagonal (al), main diagonal (ad), upper diagonal (au), and the
        !right-hand side (ab) and solve the corresponding tridiagonal matrix problem,
        !iterating over the index for y
        do j = 1, ny-2

            ab(1:nx-2) = ct(1:nx-2, j)*F_old(1:nx-2, j) &
            &       + (dtb/dy)*                    cym(1:nx-2,j  )  *F_old(1:nx-2,j+1) &
            &       - (dtb/dy)*( cyp(1:nx-2,j  ) + cym(1:nx-2,j-1) )*F_old(1:nx-2,j  ) &
            &       + (dtb/dy)*  cyp(1:nx-2,j-1)                    *F_old(1:nx-2,j-1) &
            &       + dtb*Forcing_term(1:nx-2, j)

            al(1:nx-2) =               - (dtb/dx)* cxp(0:nx-3, j)
            ad(1:nx-2) = ct(1:nx-2, j) + (dtb/dx)*(cxp(1:nx-2, j) + cxm(0:nx-3, j))
            au(1:nx-2) =               - (dtb/dx)*                  cxm(1:nx-2, j)

            ! tridiagonal solver
            call gtri(al, ad, au, ab, ax, 0, nx-1)
            F_new(:, j) = ax

        end do    
    
        ! end_time = omp_get_wtime()
        ! ! Add runtime to running total
        ! ADI_X_time = ADI_X_time + (end_time - start_time)
        ! ! Increase count of calls
        ! ADI_X_calls = ADI_X_calls + 1

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
        !real(wp) :: start_time, end_time
        integer ::  i
        
        ! start_time = omp_get_wtime()

        !Initialise arrays
        F_new   = 0._wp
        ab      = 0._wp
        al      = 0._wp
        ad      = 0._wp
        au      = 0._wp
        ax      = 0._wp

        !Apply the appropriate BCs to h based on the values specified in h_BC_params

        !North side
        select case (h_BC_params(0))
        case (1)
            ! Dirichlet Boundary conditions
            ab(ny-1) = 0._wp
            al(ny-1) = 0._wp
            ad(ny-1) = 1._wp
            au(ny-1) = 0._wp     
        case (2)
            ! Neumann Boundary conditions
            ab(ny-1) = 0._wp !-(H0(0,ny-1)-H0(0,ny-2))
            al(ny-1) = -1._wp
            ad(ny-1) = 1._wp
            au(ny-1) = 0._wp        
        case default 
            print '("!!! -- Incorrect BC Flag! -- !!!")'
            stop
        end select

        !South side
        select case (h_BC_params(2))
        case (1)
            ! Dirichlet Boundary conditions
            ab(0) = 0._wp
            al(0) = 0._wp
            ad(0) = 1._wp
            au(0) = 0._wp     
        case (2)
            ! Neumann Boundary conditions
            ab(0) = 0._wp !-(H0(0,2)-H0(0,1))
            al(0) = 0._wp
            ad(0) = -1._wp
            au(0) = 1._wp        
        case default 
            print '("!!! -- Incorrect BC Flag! -- !!!")'
            stop
        end select


        !Formulate the lower diagonal (al), main diagonal (ad), upper diagonal (au), and the
        !right-hand side (ab) and solve the corresponding tridiagonal matrix problem,
        !iterating over the index for x
        do i = 1, nx-2

            ab(1:ny-2) = ct(i, 1:ny-2)*F_old(i, 1:ny-2) &
            &       + (dtb/dx)*                    cxm(i  ,1:ny-2)  *F_old(i+1,1:ny-2) &
            &       - (dtb/dx)*( cxp(i  ,1:ny-2) + cxm(i-1,1:ny-2) )*F_old(i  ,1:ny-2) &
            &       + (dtb/dx)*  cxp(i-1,1:ny-2)                    *F_old(i-1,1:ny-2) &
            &       + dtb*Forcing_term(i, 1:ny-2)
        
            al(1:ny-2) =               - (dtb/dy)* cyp(i, 0:ny-3)
            ad(1:ny-2) = ct(i, 1:ny-2) + (dtb/dy)*(cyp(i, 1:ny-2) + cym(i, 0:ny-3))
            au(1:ny-2) =               - (dtb/dy)*                  cym(i, 1:ny-2)

            ! tridiagonal solver
            call gtri(al, ad, au, ab, ax, 0, ny-1)
            F_new(i, :) = ax

        end do
    
        ! end_time = omp_get_wtime()
        ! ! Add runtime to running total
        ! ADI_Y_time = ADI_Y_time + (end_time - start_time)
        ! ! Increase count of calls
        ! ADI_Y_calls = ADI_Y_calls + 1

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
        !real(wp) :: start_time, end_time

        ! start_time = omp_get_wtime()

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
    
    
        ! end_time = omp_get_wtime()
        ! ! Add runtime to running total
        ! gtri_time = gtri_time + (end_time - start_time)
        ! ! Increase count of calls
        ! gtri_calls = gtri_calls + 1
        
    end subroutine gtri



end module CO2GraVISim_h_routines