module CO2GraVISim_h_res_routines

    use CO2GraVISim_global
    use CO2GraVISim_input_parameters
    use CO2GraVISim_vertical_structure
    use CO2GraVISim_timer
    
    implicit none
    
    contains
    

    ! ------------------------------------------------------------------------
    subroutine h_res_update(h_res_new,h_res_old,h_cur,dh_dt,F_a_H,F_aV,dt_val,Confined_flag)
        !This subroutine updates the thickness of the residually trapped region, h_res.
        !The rate of change depends on whether h is increasing or decreasing, and whether
        !h_res is zero or positive.

        implicit none
        real(wp), dimension(0:nx-1, 0:ny-1), intent(out) ::  h_res_new
        real(wp), dimension(0:nx-1, 0:ny-1), intent(in) ::  h_res_old, h_cur, dh_dt, F_a_H, F_aV
        logical, intent(in) :: Confined_flag
        ! real(wp), dimension(0:nx-1, 0:ny-1) :: Trapping_array, No_trapping_array, h_inc_array, h_dec_array
        real(wp), dimension(0:nx-1, 0:ny-1) :: phi_res, sigma_aV, dh_res_dt, phi_h, A !, poro_cur
        real(wp), intent(in) :: dt_val
        !real(wp) :: start_time, end_time
    
        ! start_time = omp_get_wtime()

        dh_res_dt = 0._wp

        sigma_aV(:,:) = 1._wp
        sigma_aV = sign(sigma_aV , F_aV)

        ! -- Current Porosity -------------------------------------------------------------------------------------------

        call get_porosity_value(h_cur           , phi_h  )
        call get_porosity_value(h_cur+h_res_old , phi_res)

        A = 1._wp / (phi_res*( s_c_r + C_sat*(1._wp-s_c_r) ))

        where ( (h_res_old <= 0._wp) .and. (dh_dt >= 0._wp) )

            dh_res_dt = 0._wp

        elsewhere

            dh_res_dt = (-1._wp + A*C_sat*(1._wp-sigma_aV)*phi_h*(1._wp-s_a_i-s_c_r)) * dh_dt + &
            & - A*C_sat*(1._wp-sigma_aV)*F_a_H

            ! dh_res_dt = (-1._wp) * dh_dt + A*C_sat*(1._wp-sigma_aV)*F_aV

        end where




        ! ! ! where (dh_dt >= 0._wp)

        ! ! ! where (h_res_old > 0._wp)

        ! ! !     dh_res_dt = (-1._wp) * dh_dt

        ! ! ! elsewhere

        ! ! !     dh_res_dt = 0._wp

        ! ! ! endwhere

        ! ! ! elsewhere
        ! ! ! !dh_dt < 0

        ! ! !     dh_res_dt = ( -1._wp + A*C_sat*phi_h*(1._wp-s_a_i-s_c_r) ) * dh_dt

        ! ! ! end where
        

        ! ! ! ! ! ! ! !Locations where h is increasing, h is decreasing,
        ! ! ! ! ! ! ! !and where residually trapped CO2 is present
        ! ! ! ! ! ! ! h_inc_array       = merge(1._wp, 0._wp, dh_dt     >  0._wp)
        ! ! ! ! ! ! ! h_dec_array       = merge(1._wp, 0._wp, dh_dt     <= 0._wp)
        ! ! ! ! ! ! ! Trapping_array    = merge(1._wp, 0._wp, h_res_old >  0._wp)
        ! ! ! ! ! ! ! No_trapping_array = merge(1._wp, 0._wp, h_res_old <= 0._wp)

        ! ! ! ! ! F_aV = phi_h*(1._wp-s_a_i-s_c_r)*dh_dt - F_aH




      ! !   where ( (h_res_old <= 0._wp) .and. (dh_dt >= 0._wp) )

      ! !      dh_res_dt = 0._wp

      ! !   elsewhere ( (h_res_old <= 0._wp) .and. (dh_dt < 0._wp) )

      ! !      dh_res_dt = ( -1._wp + C_sat*(1._wp-s_a_i-s_c_r)/s_c_r ) * dh_dt
        
      ! !   elsewhere ( (h_res_old > 0._wp) .and. (F_aV < 0._wp) ) 

      ! !      dh_res_dt = (-1._wp) * dh_dt + (C_sat/(phi_res*s_c_r)) * F_aV 

      ! !   elsewhere ( (h_res_old > 0._wp) .and. (F_aV >= 0._wp) ) 

      ! !      dh_res_dt = (-1._wp) * dh_dt

      ! !   end where

        !Add in convective dissolution
        where (h_res_old > 0._wp)

        dh_res_dt = dh_res_dt - q_dissolve

        endwhere


        h_res_new = h_res_old + dh_res_dt*dt_val

        ! Impose that h_res = 0 at the boundaries.
        ! ** This is currently a fix to avoid occasional spurious behaviour at the edges! **
        h_res_new(0,:)     = 0._wp
        h_res_new(nx-1,:)  = 0._wp
        h_res_new(:,0)     = 0._wp
        h_res_new(:,ny-1)  = 0._wp


        !Make sure that h+h_res <= D0 - only applicable in the Confined case
        if ( Confined_flag ) then
            h_res_new = merge( h_res_new, D0 - h_cur , (h_cur + h_res_new) < D0 )
        end if

        !Make sure that h_res >= 0
        h_res_new = merge( h_res_new , 0._wp , h_res_new > 0._wp )
        
        
    
        ! end_time = omp_get_wtime()
        ! ! Add runtime to running total
        ! h_res_time = h_res_time + (end_time - start_time)
        ! ! Increase count of calls
        ! h_res_calls = h_res_calls + 1

    end subroutine h_res_update


end module CO2GraVISim_h_res_routines