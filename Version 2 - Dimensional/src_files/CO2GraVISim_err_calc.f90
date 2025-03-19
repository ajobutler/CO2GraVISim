module CO2GraVISim_err_calc

    use CO2GraVISim_global
    use CO2GraVISim_timer

    contains

    ! ------------------------------------------------------------------------
    function err(Fa, Fb, mx, my)
        !This function calculates the L2 error between two arrays.
        !This is normalised (if possible) by the maximum value of
        !one of the two arrays to give a measure of the relative error.

        implicit none
        integer, intent(in) :: mx, my
        real(wp), dimension(0:mx-1, 0:my-1), intent(in) :: Fa, Fb
        real(wp) :: err, scale
        !real(wp) :: start_time, end_time
    
        ! start_time = omp_get_wtime()

        !Determine the scale to normalise the error by. If one of the two
        !arrays is sufficiently close to 0, then set the scale to 1 - then
        !the error value will only be small if the other array is also
        !close to zero
        scale = maxval(abs(Fb))
        if (abs(scale) .le. 1e-8_wp) then
        scale = 1._wp
        end if

        err = sum((Fa - Fb)**2.)/scale
        err = dsqrt(err/(mx*my))
    
    
        ! end_time = omp_get_wtime()
        ! ! Add runtime to running total
        ! err_time = err_time + (end_time - start_time)
        ! ! Increase count of calls
        ! err_calls = err_calls + 1

    end function err

end module CO2GraVISim_err_calc