module CO2GraVISim_logging
    ! This module contains the subroutines used to write logging information (sent to STDOUT by default) 
    ! during a run of CO2GraVISim

    use CO2GraVISim_global
    use CO2GraVISim_input_parameters


    contains

    ! ---------------------------------------------------------------------------------------------
    subroutine log_title_row

        implicit none
        character(100) :: title_str
        integer :: len_title

        write(title_str,'("| ", A7, I0, " | ", A9, " | ", A10, " | ", A14, " | " , A18, " | " , A11 , " |" )')&
        & 'Save #/', np-1, 't [days]', 'timestep #', 'Avg. dt [days]', 'Mass injected [kg]', 'M_tot/M_inj'
    
        len_title = len_trim(title_str)

        write(*,'(A)') REPEAT('-', len_title )
        write(*,'(A)') trim(title_str)
        write(*,'(A)') REPEAT('-', len_title )

    end subroutine log_title_row
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine log_details_line(line_character, text_string)
        
        implicit none
        character(1), intent(in) :: line_character
        character(*), intent(in) :: text_string
        integer :: np_width, n_columns
        integer :: total_length
        character(1000) :: details_string
        
        np_width = floor(log10(real(np-1))) + 1

        n_columns = 6

        total_length = 2*2 + (n_columns - 1)*3 + 7 + np_width + 9 + 10 + 14 + 18 + 11

        call create_details_line(details_string, line_character, text_string, total_length)
    
        write(*,'(A)') details_string(1:total_length)

    end subroutine log_details_line
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine log_runtime(runtime)

        implicit none
        real(wp), intent(in) :: runtime
        character(100) :: runtime_string

        write(runtime_string, '(A,F16.3,A)') 'Run complete. Runtime = ', runtime, ' seconds.'

        call log_details_line('-','')
        call log_details_line(' ','')
        call log_details_line('=','')
        call log_details_line('=', runtime_string)
        call log_details_line('=','')


    end subroutine log_runtime
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine log_solver_mode(Confined_flag)

        implicit none
        logical, intent(in) :: Confined_flag
        character(100) :: mode_string

        if (Confined_flag) then
            mode_string = "Running solver: CONFINED MODE"
         else
            mode_string = "Running solver: UNCONFINED MODE"
         end if


         call log_details_line(' ','')
         call log_details_line(' ','')
         call log_details_line('=','')
         call log_details_line('=',mode_string)
         call log_details_line('=','')
         call log_details_line(' ','')
         call log_details_line(' ','')

    end subroutine log_solver_mode
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine log_run_parameters(V_initial)

        implicit none
        real(wp), intent(in) :: V_initial
        character(100) :: fmt_str_int, fmt_str_real, fmt_str_real_2

        fmt_str_int    = '(A17, " = ", I10)'
        fmt_str_real   = '(A17, " = ",F10.3, 1X, "[", A3, "]")'
        fmt_str_real_2 = '(A17, " = ",E10.3E3, 1X, "[", A3, "]")'


        call log_details_line('-','Run Parameters')

        write(*,fmt_str_int)   'nx           ', nx
        write(*,fmt_str_int)   'ny           ', ny
        write(*,fmt_str_int)   'nz           ', nz
        
        write(*,fmt_str_real) 'dx           ', dx * Length_scale ,  ' m '
        write(*,fmt_str_real) 'dy           ', dy * Length_scale ,  ' m '
        write(*,fmt_str_real) 'M            ', M ,                  ' - '
        write(*,fmt_str_real) 'Gamma_val    ', Gamma_val ,          ' - '
        write(*,fmt_str_real) 's_c_r        ', s_c_r ,              ' - '
        write(*,fmt_str_real) 's_a_i        ', s_a_i ,              ' - '
        write(*,fmt_str_real) 'C_sat        ', C_sat ,              ' - '
        write(*,fmt_str_real) 'q_dissolve   ', q_dissolve ,         ' - '
        write(*,fmt_str_real) 'krn_mobile   ', krn_mobile ,         ' - '
        write(*,fmt_str_real) 'krw_residual ', krw_residual ,       ' - '

        write(*,'(A)') ''
        write(*,fmt_str_real)   't0           ', t0      * Time_scale * seconds_to_days,    'day'
        write(*,fmt_str_real_2) 'dt_init      ', dt_init * Time_scale * seconds_to_days,    'day'
        write(*,fmt_str_real_2) 'dt_min       ', dt_min  * Time_scale * seconds_to_days,    'day'
        write(*,fmt_str_real_2) 'err_max_t    ', errmax_t_step,                             ' - '
        write(*,fmt_str_real_2) 'err_max_P    ', errmax_P_iter,                             ' - '
        write(*,fmt_str_real)   'P_iter_omega ', P_iter_omega,                              ' - '
        write(*,fmt_str_int)    'P_iter_check ', P_iter_check

        write(*,'(A)') ''
        write(*,fmt_str_real_2) 'Initial mass ', V_initial * rho_c * Volume_scale , 'kg'
        call log_details_line('-','')


    end subroutine log_run_parameters
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine log_save_info(plot_count, t, iter, avg_timestep, V_injected, V_total)

        implicit none
        integer, intent(in) :: plot_count, iter
        real(wp), intent(in) :: t, avg_timestep, V_injected, V_total
        real(wp) :: M_injected, Volume_ratio
        character(1000) :: save_format, save_format_0
        integer :: np_width

        !Mass injected
        M_injected = rho_c * V_injected * Volume_scale

        !Ratio of total calculated volume to injected volume
        !(equivalent to mass ratio, since density is constant in this model)
        Volume_ratio = V_total / V_injected
            
        !Repeat the table header every 20 plot times
        if ((plot_count>0) .and. mod(plot_count,20)==0) then
            call log_title_row
        end if


        ! format string for a table row written to standard out
        np_width = floor(log10(real(np-1))) + 1
        write(save_format,'(A,I0,A,I0,A,I0,A)')&
        & "('| ', I", 7+np_width , ", ' | ', F9.2, ' | ', I10, ' | ', F14.5,' | ', E18.5E3, ' | ', F11.5, ' |')"
        !Modified format for when V_injected = 0
        write(save_format_0,'(A,I0,A,I0,A,A,A)')&
            & "('| ', I", 7+np_width , ", ' | ', F9.2, ' | ', I10, ' | ', F14.5,' | ', E18.5E3, ' | ', A11, ' |')"


        if (V_injected > 0._wp) then
            write(*, save_format) &
                &plot_count, &
                &t * Time_scale * seconds_to_days, &
                &iter, &
                &avg_timestep, &
                & M_injected, &
                & Volume_ratio
        else
            !Nothing has been injected yet, so the volume ratio wouldn't make sense
                write(*, save_format_0) &
                &plot_count, &
                &t * Time_scale * seconds_to_days, &
                &iter, &
                &avg_timestep, &
                & M_injected, &
                & '   -N/A-   '
                !& '-'
        end if


    end subroutine log_save_info
    ! ---------------------------------------------------------------------------------------------
    

        ! ---------------------------------------------------------------------------------------------
    subroutine create_details_line(details_string, line_character, text_string, total_length)

        implicit none
        character(1), intent(in) :: line_character
        character(*), intent(in) :: text_string
        integer, intent(in) :: total_length
        integer :: text_length, left_length, right_length
        character(total_length) :: details_string
    
    
        text_length = len_trim(text_string)
    
        if (text_length .eq. 0) then
            !Blank text_string, used to print a line of line_character with no gaps
            write(details_string,'(A)') repeat(line_character, total_length)
        else
            ! Put half (after rounding) of the hyphens on the left, once we've accounted for the title string
            ! and the padding spaces
            left_length = ( total_length - (text_length+2) )/2 !Integer division should round as I want
            !Put the remainder of the hyphens on the right
            right_length = total_length - text_length - left_length - 2
    
            !Print this full string to details_string
            write(details_string,'(A,1X,A,1X,A)') &
            &                   repeat(line_character, left_length) , trim(text_string) , repeat(line_character, right_length)
    
        end if
    
    
        endsubroutine create_details_line
        ! ---------------------------------------------------------------------------------------------
    

end module CO2GraVISim_logging