module CO2GraVISim_iter_calculations

    use CO2GraVISim_global
    use CO2GraVISim_input_parameters
    use CO2GraVISim_logging
    use CO2GraVISim_vertical_structure

    implicit none
    integer :: fid_Iter_info, fid_P_iter, fid_Volume, fid_Plume_size

    integer :: iter_save_old = 0

    character(1000) :: iter_title_format = '(A5,2X,7(A15,2X))'
    character(1000) :: iter_data_format = '(I5,1X,2(G16.8,1X),5(E16.8E3,1X))'

    character(1000) :: P_title_format = '(2(A16,1X),3(A20,2X),1X,A16)'
    character(1000) :: P_data_format = '(2(G16.8,1X),3(4I5,2X),1X,E16.8E3)'

    character(1000) :: Volume_title_format = '(A5,1X,(A9,1X),8(A16,1X))'
    character(1000) :: Volume_data_format = '(I5,1X,(F9.2,1X),8(E16.8E3,1X))'

    character(1000) :: Plume_title_format = '(A5,2X,4(A15,2X))'
    character(1000) :: Plume_data_format = '(I5,1X,(G16.8,1X),3(E16.8E3,1X))'

contains

    ! ------------------------------------------------------------------------
    subroutine iter_file_open(Output_folder)

        implicit none
        character(1000), intent(in) :: Output_folder
        character(1000) :: Iter_info_file, P_iter_file, V_file, Plume_size_file


        ! Record of general simulation information
        Iter_info_file = trim(Output_folder) // "/Other/Iter_info.txt"

        open (newunit=fid_Iter_info,  file=Iter_info_file, action='write')

        write(fid_Iter_info,iter_title_format) &
        & 'it', 't [days]', 'dt_old [days]', 'err_h', 'err_P', 'max_h [m]', 'max_h_res [m]', 'max_P [Pa]'


        ! Record of Pressure iterations
        P_iter_file = trim(Output_folder) // "/Other/P_iter_info.txt"

        open (newunit=fid_P_iter,  file=P_iter_file, action='write')

        write(fid_P_iter,P_title_format) &
        & 't [days]', 'dt [days]', 'P_iter_vals_1', 'P_iter_vals_2a', 'P_iter_vals_2b', 'err_P [Pa]'


        ! Record of volume information
        V_file   = trim(Output_folder) // "/Other/Volumes.txt"

        open (newunit=fid_Volume,     file=V_file,   action='write')

        write (fid_Volume, Volume_title_format) &
        & 'it', 't [days]', 'V_inj [m^3]', 'V_tot [m^3]', 'V_mob_sum [m^3]', 'V_trp_sum [m^3]', 'V_dsln_sum [m^3]', &
        &   'V_boundary [m^3]', 'V_mob_max [m^3]', 'V_trp_max [m^3]'


        ! Record of the plume size
        Plume_size_file = trim(Output_folder) // "/Other/Plume_size.txt"

        open (newunit=fid_Plume_size,  file=Plume_size_file, action='write')

        write(fid_Plume_size,Plume_title_format) &
        & 'it', 't [days]', 'max_h [m]', 'width_x [m]', 'width_y [m]'

  

    end subroutine iter_file_open
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine iter_calculations(h,h_res,P,iter_num,t,dt,dt_ideal,err_h,err_P, &
        & P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b)

        implicit none
        real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: h, h_res, P
        integer, intent(in) :: iter_num
        real(wp), intent(in) :: t, dt, dt_ideal, err_h, err_P
        integer, dimension(0:3), intent(in) :: P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b
        real(wp) :: t_scaled, dt_scaled
        real(wp) :: width_x, width_y, max_h, max_h_res, max_P

        max_h       = maxval(h)     * Length_scale
        max_h_res   = maxval(h_res) * Length_scale
        max_P       = maxval(P)     * Pressure_scale

        t_scaled    = t  * Time_scale * seconds_to_days
        dt_scaled   = dt * Time_scale * seconds_to_days

        ! Record general iteration information
        write(fid_Iter_info, iter_data_format) &
        & iter_num, t_scaled, dt_scaled, err_h, err_P, max_h, max_h_res, max_P
        
        ! Record Pressure iteration information
        write(fid_P_iter, P_data_format) &
        & t_scaled, dt_scaled, P_iter_vals_1, P_iter_vals_2a, P_iter_vals_2b, err_P

        ! Record plume width in x and y (useful for when it's one cohesive mass)
        call plume_size(h, width_x, width_y)
        write(fid_Plume_size, Plume_data_format) &
        & iter_num, t_scaled, maxval(h), width_x, width_y


    
    end subroutine iter_calculations
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine P_flux_change_record(t, P_iter_count)

        implicit none
        real(wp), intent(in) :: t
        integer, intent(in) :: P_iter_count
        real(wp) :: t_scaled

        t_scaled = t * Time_scale * seconds_to_days
        
        write(fid_P_iter,'(G16.8, 1X, A16, 1X, I5)') t_scaled, '', P_iter_count

    end subroutine P_flux_change_record
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine iter_file_close

        implicit none

        !Close the various files that have been opened here
        close(fid_Iter_info)
        close(fid_P_iter)
        close(fid_Volume)
        close(fid_Plume_size)

    end subroutine iter_file_close
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine plume_size(h, width_x, width_y)

        implicit none
        real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: h
        real(wp), intent(out) :: width_x, width_y
        real(wp) :: threshold
        integer :: kx_min, kx_max, kx_min_temp, kx_max_temp
        integer :: ky_min, ky_max, ky_min_temp, ky_max_temp
        integer ::count

        threshold = 1e-10

        ! -- Scan over y to find x width -------------------------------------------------
        kx_min = 0
        kx_max = 0
        kx_min_temp = 0
        kx_max_temp = 0
        width_x = 0._wp

        do count = 0, ny-1

            kx_min_temp = findloc(h(:,count) .ge. threshold, .TRUE. , dim=1, back=.false.)
            kx_max_temp = findloc(h(:,count) .ge. threshold, .TRUE. , dim=1, back=.true. )

            kx_min = max(kx_min_temp, kx_min)
            kx_max = max(kx_max_temp, kx_max)

        end do

        width_x = dx*(kx_max - kx_min)

        ! -- Scan over x to find y width -------------------------------------------------
        ky_min = 0
        ky_max = 0
        ky_min_temp = 0
        ky_max_temp = 0
        width_y = 0._wp

        do count = 0, nx-1

            ky_min_temp = findloc(h(count,:) .ge. threshold, .TRUE. , dim=1, back=.false.)
            ky_max_temp = findloc(h(count,:) .ge. threshold, .TRUE. , dim=1, back=.true. )

            ky_min = max(ky_min_temp, ky_min)
            ky_max = max(ky_max_temp, ky_max)

        end do

        width_y = dy*(ky_max - ky_min)

        ! Redimensionalise
        width_x = width_x * Length_scale
        width_y = width_y * Length_scale


    end subroutine plume_size
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine save_iter_output(h_array, h_res_array, P_array, V_mobile_array, V_trapped_array, &
    &   Max_mobile_thickness, Max_pressure, Min_pressure, h_cur, h_res_cur, P_cur, &
    &   t, tplot, plot_count, plot_times, V_injected, V_total_dsln_loss, V_boundary, iter)
        !This subroutine calculates the mobile and trapped volumes of CO2 after each
        !successful timestep, and saves them to an output file.
        !If the timestepping routine has also passed a planned output time, the calculated
        !profiles h, h_res, and P are saved into their respective output arrays.
 
        implicit none
        real(wp), dimension(0:nx-1, 0:ny-1, 0:np-1), intent(inout) :: h_array, h_res_array, P_array
        real(wp), dimension(0:nx-1, 0:ny-1, 0:np-1), intent(inout) :: V_mobile_array, V_trapped_array
        real(wp), dimension(0:nx-1, 0:ny-1), intent(inout) :: Max_mobile_thickness, Max_pressure, Min_pressure
        real(wp), dimension(0:nx-1, 0:ny-1), intent(in) :: h_cur, h_res_cur, P_cur
        real(wp), dimension(0:np-1), intent(inout) :: plot_times
        real(wp), intent(inout) :: tplot
        real(wp), intent(in) :: t, V_injected, V_total_dsln_loss, V_boundary
        integer, intent(in) :: iter
        integer, intent(inout) :: plot_count
        real(wp) :: V_mobile_sum, V_trapped_sum, V_mobile_local_max, V_trapped_local_max, V_total_dsln, V_total
        real(wp) :: V_mob_dslv_sum, V_trp_dslv_sum
        ! real(wp), dimension(0:nx-1,0:ny-1) :: poro_integral_cur, poro_integral_res
        real(wp), dimension(0:nx-1,0:ny-1) :: V_mobile, V_trapped, V_mobile_dissolved , V_trapped_dissolved
        real(wp), dimension(0:7) :: Volume_data
        real(wp) :: t_dim, avg_timestep
        integer :: i
        !real(wp) :: start_time, end_time



        ! start_time = omp_get_wtime()

        t_dim = t * Time_scale * seconds_to_days !Express in days rather than seconds


        call calculate_volumes(h_cur, h_res_cur, V_mobile, V_trapped, V_mobile_dissolved, V_trapped_dissolved)

        V_mobile_sum  = sum( V_mobile  )
        V_trapped_sum = sum( V_trapped )

        V_mob_dslv_sum = sum( V_mobile_dissolved  )
        V_trp_dslv_sum = sum( V_trapped_dissolved )

        V_mobile_local_max   = maxval( V_mobile  )
        V_trapped_local_max  = maxval( V_trapped )

        V_total_dsln = V_mob_dslv_sum + V_trp_dslv_sum + V_total_dsln_loss

        V_total = V_mobile_sum + V_trapped_sum + V_total_dsln + V_boundary

        ! Volume_data = [ t, V_mobile_sum, V_trapped_sum, V_injected, V_mobile_local_max, V_trapped_local_max, V_total_dsln ]
        Volume_data = [ V_injected, V_total, V_mobile_sum, V_trapped_sum, V_total_dsln, &
        &                                      V_boundary, V_mobile_local_max, V_trapped_local_max]

        ! Dimensionalise
        Volume_data = Volume_data * Volume_scale

        ! Write volume data to output text file
        write (fid_Volume, Volume_data_format) iter, t_dim, (Volume_data(i), i=0, size(Volume_data)-1)

        ! Update the presence array, tracking which cells have encountered CO2
        Max_mobile_thickness = merge( h_cur, Max_mobile_thickness, h_cur > Max_mobile_thickness )

        Max_pressure = merge( P_cur, Max_pressure, P_cur > Max_pressure )
        Min_pressure = merge( P_cur, Min_pressure, P_cur < Min_pressure )

 
        ! Check if we've reached a plot time
        if (plot_count .le. np-1 .and. t >= tplot) then

            !Save data for current plot time
            !These get redimensionalised outside of the solver

            h_array(:,:,plot_count)       = h_cur(:,:)
            h_res_array(:,:,plot_count)   = h_res_cur(:,:)
            P_array(:,:,plot_count)       = P_cur(:,:)
            plot_times(plot_count)        = t

            V_mobile_array(:,:,plot_count)  = V_mobile
            V_trapped_array(:,:,plot_count) = V_trapped

            !Calculate the average (dimensional) length of timestep between now and the previous plot time
            if (plot_count == 0) then
                avg_timestep = 0._wp
            else
                avg_timestep = ((t - target_plot_times(plot_count-1)) * Time_scale * seconds_to_days) / (iter - iter_save_old)
            end if

            call log_save_info(plot_count, t, iter, avg_timestep, V_injected, V_total)





            if (plot_count .lt. np-1) then
                !Update next target plot time
                tplot = target_plot_times(plot_count+1)
            end if

            !Update plot number
            plot_count = plot_count + 1

            !Update iteration counter for 'previous' plot time
            iter_save_old = iter

        end if
     
       ! end_time = omp_get_wtime()
       ! ! Add runtime to running total
       ! save_iter_time = save_iter_time + (end_time - start_time)
       ! ! Increase count of calls
       ! save_iter_calls = save_iter_calls + 1
 
    end subroutine save_iter_output
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine calculate_volumes(h_cur, h_res_cur, V_mobile, V_trapped, V_mobile_dissolved, V_trapped_dissolved)
        !This subroutine calculates the local volumes of free-phase and dissolved CO2
        !in the mobile (H<=z<=H+h) and residually trapped (H+h<=z<=H+h+h_res) regions

        implicit none
        real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h_cur, h_res_cur
        real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: V_mobile, V_trapped, V_mobile_dissolved , V_trapped_dissolved
        real(wp), dimension(0:nx-1,0:ny-1) :: poro_integral_cur, poro_integral_res, poro_integral_trp

        poro_integral_cur = 0._wp
        poro_integral_res = 0._wp
        poro_integral_trp = 0._wp

        call get_porosity_integral(h_cur            , poro_integral_cur) !Integral from H to H+h
        call get_porosity_integral(h_cur + h_res_cur, poro_integral_res) !Integral from H to H+h+h_res


        !Integral just over the residually trapped region (H+h <= z <= H+h+h_res)
        poro_integral_trp = poro_integral_res - poro_integral_cur

        V_mobile    = ( poro_integral_cur )*(1._wp - s_a_i)*dx*dy
        V_trapped   = ( poro_integral_trp )*s_c_r          *dx*dy
            
        V_mobile_dissolved    = C_sat*( poro_integral_cur )*s_a_i          *dx*dy
        V_trapped_dissolved   = C_sat*( poro_integral_trp )*(1._wp - s_c_r)*dx*dy

    end subroutine calculate_volumes
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine calculate_dissolution_loss(V_dissolution_loss,h_new,h_old,h_res_new,h_res_old,F_a_V,dt)
        !This subroutine calculates the amount of of CO2 lost to dissolution - either via
        !convective dissolution scouring at the base of the CO2 column, or the vertical flux 
        !of saturated ambient fluid at the base of the CO2 column into the ambient bulk, 
        !that is then lost from the record 
        
        implicit none
        real(wp), dimension(0:nx-1,0:ny-1), intent(inout) :: V_dissolution_loss
        real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h_new, h_old, h_res_new, h_res_old, F_a_V
        real(wp), intent(in) :: dt
        real(wp), dimension(0:nx-1,0:ny-1) :: poro_h, poro_h_res
        real(wp), dimension(0:nx-1,0:ny-1) :: V_dh_free, V_dh_sat, V_dhres_free, V_dhres_sat
        real(wp), dimension(0:nx-1,0:ny-1) :: V_mob, V_mob_sat, V_mob_conv, V_res, V_res_sat, V_res_conv

        V_dh_free      = 0._wp
        V_dh_sat       = 0._wp
        V_dhres_free   = 0._wp
        V_dhres_sat    = 0._wp
        
        V_mob          = 0._wp
        V_mob_sat      = 0._wp
        V_mob_conv     = 0._wp

        V_res          = 0._wp
        V_res_sat      = 0._wp
        V_res_conv     = 0._wp

        call get_porosity_value(h_old,     poro_h    )
        call get_porosity_value(h_res_old, poro_h_res)

        !! Residual Trapping Region !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        where (h_res_old > 0._wp)
        V_dhres_free = -poro_h_res*( s_c_r             )*(h_res_new - h_res_old) !Change in free-phase trapped volume 
        V_dhres_sat  = -poro_h_res*(C_sat*(1._wp-s_c_r))*(h_res_new - h_res_old) !Change in dissolved trapped volume 

        V_res_sat  = C_sat * abs(F_a_V) * dt
        V_res_conv = q_dissolve * dt
        
        endwhere
        V_dhres_free = merge( V_dhres_free , 0._wp , V_dhres_free > 0._wp )
        V_dhres_sat  = merge( V_dhres_sat  , 0._wp , V_dhres_sat  > 0._wp )
        V_res_sat  = merge( V_res_sat  , V_dhres_sat  , V_res_sat  < V_dhres_sat  ) !Actual volume change is an upper bound
        V_res_conv = merge( V_res_conv , V_dhres_free , V_res_conv < V_dhres_free ) !Actual volume change is an upper bound

        V_res = V_res_sat + V_res_conv


        !! Region of exposed mobile CO2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        where ( (h_old > 0._wp) .and. (h_res_old < 0._wp) )
        V_dh_free = -poro_h*( (1._wp - s_a_i) )*(h_new - h_old) !Change in free-phase trapped volume 
        ! V_dh_sat  = -poro_h*( C_sat*s_a_i     )*(h_new - h_old) !Change in dissolved trapped volume 

        ! V_mob_sat  = C_sat * abs(F_a_V) * dt
        V_mob_conv = q_dissolve * dt
        
        endwhere
        V_dh_free = merge( V_dh_free , 0._wp , V_dh_free > 0._wp )
        ! V_dh_sat  = merge( V_dh_sat  , 0._wp , V_dh_sat  > 0._wp )
        ! V_mob_sat  = merge( V_mob_sat  , V_dh_sat  , V_mob_sat  < V_dh_sat  ) !Actual volume change is an upper bound
        V_mob_conv = merge( V_mob_conv , V_dh_free , V_mob_conv < V_dh_free ) !Actual volume change is an upper bound

        V_mob = V_mob_sat + V_mob_conv

        V_dissolution_loss = V_mob + V_res 



        ! V_saturated = 0._wp 
        ! !Vertical flux of saturated ambient
        ! where (h_res_old > 0._wp)
        !    V_saturated = C_sat * abs(F_a_V) * dt
        ! endwhere
        ! V_saturated = merge( V_saturated , V_delta_h_res , V_saturated > V_delta_h_res )
        
        ! !Convective dissolution
        ! V_convective = 0._wp
        ! where ((h_old+h_res_old) > 0._wp)
        !    V_convective = q_dissolve * dt
        ! endwhere
        ! V_saturated = merge( V_saturated , V_delta_h_res , V_saturated > V_delta_h_res )

        
        ! V_dissolution_loss = V_dissolution_loss + V_saturated

    
    
    end subroutine calculate_dissolution_loss
    ! ------------------------------------------------------------------------


    ! ------------------------------------------------------------------------
    subroutine calculate_boundary_flux(h, P, dt_val, V_boundary, direction)

       implicit none   

       real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h, P
       real(wp), intent(in) :: dt_val
       real(wp), intent(inout) :: V_boundary
       character(1), intent(in) :: direction
       real(wp), dimension(0:nx-1,0:ny-1) :: F, int_k_c
       real(wp) :: Q_N, Q_E, Q_S, Q_W

       ! The Darcy flux for the CO2, vertically integrated across the mobile region, is
       ! int_{H}^{H+h} u_c(x,y,z) dz = - int_{H}^{H+h} k(x,y,z) * grad_{x,y}( P + Gamma*(H0+h) ) dz
       !                             = - int_k_c * grad_{x,y}( P + Gamma*(H0+h) )


       ! call get_permeability_average_mobile(h,kappa_c)
       call get_permeability_integral(h,int_k_c)

       F = P + Gamma_val*(H0 + h)

       select case (direction)

       case ('x')
          !Solving in the x direction here, so we apply the East and West boundary conditions
          ! East boundary : (i,j) = (nx-1, j)
          Q_E = sum( - 0.5_wp*(int_k_c(nx-2,:) + int_k_c(nx-1,:)) * (F(nx-1,:) - F(nx-2,:))/dx ) * dy
          ! West boundary : (i,j) = (0, j)
          Q_W = sum( - 0.5_wp*(int_k_c(0   ,:) + int_k_c(1   ,:)) * (F(0   ,:) - F(1   ,:))/dx ) * dy

          Q_N = 0._wp
          Q_S = 0._wp

       case ('y')
          !Solving in the y direction here, so we apply the North and South boundary conditions
          ! North boundary : (i,j) = (i, ny-1)
          Q_N = sum( - 0.5_wp*(int_k_c(:,ny-2) + int_k_c(:,ny-1)) * (F(:,ny-1) - F(:,ny-2))/dy ) * dx
          ! South boundary : (i,j) = (i, 0)
          Q_S = sum( - 0.5_wp*(int_k_c(:,0   ) + int_k_c(:,1   )) * (F(:,0   ) - F(:,1   ))/dy ) * dx

          Q_E = 0._wp
          Q_W = 0._wp

       case default
          write(*,'(A)') 'Wrong direction used in calculate_boundary_flux; should be x or y !'
          write(*,'(A,A)') 'direction used: ', trim(direction)
          error stop

       end select

       V_boundary = V_boundary + (Q_N + Q_E + Q_S + Q_W)*dt_val

       ! write(*,'(A,E9.2E3)') 'V_boundary = ', V_boundary



       !! Do I need to incorporate dissolution in here as well??

       !! Do I need to run this for the initial profile?


       ! ! ! Total volume lost through the boundaries
       ! ! ! write(*,'(4(A6,E9.2E3,1X))') 'Q_N = ', Q_N, 'Q_E = ', Q_E, 'Q_S = ', Q_S, 'Q_W = ', Q_W
       ! ! ! write(*,'(2(A,E9.2E3,1X))') 'Q_W = ', sum( int_k_c(0,:) ), 'dF/dn_W = ', sum( (F(0,:)-F(1,:))/dx ) 
       ! ! V_boundary = V_boundary + (Q_N + Q_E + Q_S + Q_W)*dt_val

    
    end subroutine calculate_boundary_flux


    

end module CO2GraVISim_iter_calculations