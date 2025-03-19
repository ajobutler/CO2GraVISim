module CO2GraVISim_vertical_structure

   use CO2GraVISim_global
   use CO2GraVISim_input_parameters
   use CO2GraVISim_timer
   use OMP_LIB

   implicit none

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_porosity_value(h,poro)

      implicit none

      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: poro
      logical :: extrapolate_flag

      extrapolate_flag = .False.
      ! call vertical_interpolation(Porosity, h, D0, poro, extrapolate_flag)

      call piecewise_constant_interpolation(Porosity, h, D0, poro)

   end subroutine get_porosity_value
   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine get_porosity_integral(h,poro_integral)

      implicit none

      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: poro_integral
      logical :: extrapolate_flag

      extrapolate_flag = .True.
      ! call vertical_interpolation(poro_CInt, h, D0, poro_integral, extrapolate_flag)

      call piecewise_constant_integration(Porosity, h, D0, poro_integral)

   end subroutine get_porosity_integral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine get_permeability_value(h,perm)

      implicit none

      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: perm
      logical :: extrapolate_flag

      extrapolate_flag = .False.
      call vertical_interpolation(Permeability, h, D0, perm, extrapolate_flag)

      ! call piecewise_constant_interpolation(Permeability, h, D0, perm)


   end subroutine get_permeability_value

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_permeability_integral(h,perm_integral)

      implicit none

      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: perm_integral
      logical :: extrapolate_flag

      extrapolate_flag = .True.
      call vertical_interpolation(perm_CInt, h, D0, perm_integral,extrapolate_flag)

   end subroutine get_permeability_integral

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine get_permeability_average_mobile(h,kappa_c)

   !! kappa = (1/(z_upper-z_lower)) int_{z_lower}^{z_upper} k dz
   
   implicit none
   
   real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h
   real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: kappa_c
   real(wp), dimension(0:nx-1,0:ny-1) :: int_h, int_approx
   real(wp) :: threshold, dz
   integer :: i,j
   
   ! Threshold value for where we approximate the average value by the function value
   ! (rather than divide by a small value of h)
   threshold = 1e-3_wp
   
   kappa_c = 0._wp
   
   ! Integrate the permeability from 0 to h (relative to the caprock)
   call get_permeability_integral(h,int_h)
   ! Get the approximate value to use if h is too small for the average integral
   call get_permeability_value( (0._wp+h)/2._wp , int_approx )
   
   ! Build the average integral array
   do j = 0,ny-1
       do i = 0,nx-1
           dz = (h(i,j) - 0._wp)
   
           if ( dz < threshold ) then
               ! h is too small to do (1/h)*integral(0->h), so just use the value at h/2 instead
               kappa_c(i,j) = int_approx(i,j)
           else
               ! Calculate the average value via (1/h)*integral(0->h)
               kappa_c(i,j) = ( int_h(i,j) - 0._wp ) / dz
           endif
   
       enddo
   enddo

   ! Scale the permeability in the mobile region by the corresponding
   ! relative permeability value, krn_mobile
   kappa_c = krn_mobile * kappa_c
   
   
   end subroutine get_permeability_average_mobile    
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   subroutine get_permeability_average_ambient(h,h_res,kappa_a)
   
   !! kappa = (1/(z_upper-z_lower)) int_{z_lower}^{z_upper} k dz
   
   implicit none
   
   real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h, h_res
   real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: kappa_a
   real(wp), dimension(0:nx-1,0:ny-1) :: int_res, int_amb, int_approx, int_temp
   real(wp) :: threshold, dz
   integer :: i,j
   
   ! Threshold value for where we approximate the average value by the function value
   ! (rather than divide by a small value of (D0-h))
   threshold = 1e-3_wp
   
   kappa_a = 0._wp
   
   ! Integrate the permeability over the residual trapping region, i.e.
   ! from h to h_res (as [0->h+h_res] - [0->h])
   call get_permeability_integral(h,     int_res )
   call get_permeability_integral(h_res, int_temp) !Store this one for use below
   int_res = int_temp - int_res
   
   ! Integrate the permeability over the ambient region, i.e.
   ! from h+h_res to D0 (as [0->D0] - [0->h+h_res])
   call get_permeability_integral(D0, int_amb)
   int_amb = int_amb - int_temp
   
   ! Scale the integral over the residual trapping region by the corresponding
   ! relative permeability value, krw_residual
   int_res = krw_residual * int_res
   
   ! Combine the two integrals: [h->h+h_res] + [h+h_res->D0]
   int_amb = int_res + int_amb
   
   ! Get the approximate value to use if D0-h is too small for the average integral
   call get_permeability_value( (h+D0)/2._wp , int_approx )
   ! I should actually do something a bit more clever here and approximate it by
   ! (1 - (1-krw_residual*h_res)/(D0-h) ) * k(h+(D0-h)/2) ,
   ! but I'm not yet sure how to deal with the "small"/"small" fraction that could occur here
   ! and which I'm trying to avoid via this process.
   
   ! Build the average integral array
   do j = 0,ny-1
       do i = 0,nx-1
           dz = D0(i,j) - h(i,j)
   
           if ( dz < threshold ) then
               ! (D0-(h+h_res)) is too small to do (1/(D0-(h+h_res)))*integral(h+h_res->D0), 
               ! so just use the value at (D0+(h+h_res))/2 instead
               kappa_a(i,j) = int_approx(i,j)
           else
               ! Calculate the average value via (1/(D0-h))*integral(h+h_res->D0)
               kappa_a(i,j) = int_amb(i,j) / dz
           endif
   
       enddo
   enddo
   
   end subroutine get_permeability_average_ambient


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine vertical_interpolation(Array, h, D, interp_values, extrapolate_flag)

      implicit none
      real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: Array
      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h, D
      logical, intent(in) :: extrapolate_flag
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: interp_values
      real(wp) :: r, d_n, d_p, d_i, z_prev, z_next, Array_n, Array_p
      real(wp) :: start_time, end_time
      integer :: i, j, k_next, k_prev
  
      start_time = omp_get_wtime()
  
  
      !Preallocate values
      z_next          = 0._wp
      z_prev          = 0._wp
      Array_n         = 0._wp
      Array_p         = 0._wp
      d_n             = 0._wp
      d_p             = 0._wp
      d_i             = 0._wp
      interp_values   = 0._wp
 
  
      !Determine relevant reference points and values
      do j=0, ny-1
          do i=0, nx-1

            !! Need to allow for h/D>1 in the Unconfined case.
            r = h(i,j)/D(i,j)
            !   r = max( r , 0._wp )
            !   r = min( r , 1._wp )

            if ( r <= 0 ) then
               !!! Use the value corresponding to the ceiling
               interp_values(i,j) = Array(i, j, 0)

            elseif ( r > 1 ) then
               !!! h>D, which is allowable in the Unconfined case.

               if ( .Not. extrapolate_flag ) then
                  !We don't want to extrapolate (e.g. for the porosity)
                  
                  !Use the final value in the column
                  interp_values(i,j) = Array(i,j,nz-1)

               else
                  !We do want to extrapolate (e.g. for the cumulative integrals)
                  
                  !Linearly extrapolate based on the final interval
                  k_next = nz-1
                  k_prev = nz-2

                  z_next = Z_layers(k_next)
                  z_prev = Z_layers(k_prev)
                  Array_n = Array(i, j, k_next)
                  Array_p = Array(i, j, k_prev)
                  
                  !Calculate distances involved (between interval endpoints and interpolation point)
                  d_n = z_next - r
                  d_p = r      - z_prev
                  d_i = z_next - z_prev
                  
                  interp_values(i,j) = ( d_n*Array_p + d_p*Array_n ) / d_i

               end if
               

            else
               !!! Linearly interpolate based on the interval containing r
  
               !Find the first index where Z_layers >= the current r value
               !!! Need to subtract one because findloc returns a 1-indexed value
               k_next = findloc( Z_layers .ge. r , .TRUE. , Dim=1 ) - 1
               k_prev = max( k_next - 1 , 0 )
   
               z_next = Z_layers(k_next)
               z_prev = Z_layers(k_prev)
               Array_n = Array(i, j, k_next)
               Array_p = Array(i, j, k_prev)
               
               !Calculate distances involved (between interval endpoints and interpolation point)
               d_n = z_next - r
               d_p = r      - z_prev
               d_i = z_next - z_prev
   
               if ( d_i > 0._wp ) then
                     interp_values(i,j) = ( d_n*Array_p + d_p*Array_n ) / d_i
               else
                     interp_values(i,j) = Array_n
               end if

            end if  
  
          end do
      end do
  
      end_time = omp_get_wtime()
      ! Add runtime to running total
      vert_interp_time = vert_interp_time + (end_time - start_time)
      ! Increase count of calls
      vert_interp_calls = vert_interp_calls + 1
  
  
   end subroutine vertical_interpolation

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   subroutine piecewise_constant_interpolation(Array, h, D, interp_values)
    
      !Interpolate by simply returning the value corresponding to the z_layer value just below
  
      implicit none
      real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: Array
      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h, D
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: interp_values
      real(wp) :: r
      integer :: i, j, k_prev
      
      
      do j=0, ny-1
          do i=0, nx-1
              
              r = h(i,j)/D(i,j)
              
              if ( r <= 0 ) then
                  !!! Use the value corresponding to the ceiling
                  k_prev = 0
  
              elseif ( r > 1 ) then
                  !!! h>D, which is allowable in the Unconfined case.
                  !!! Use the value at the basement
                  k_prev = nz-1
                  
              else
                  !!! Find the index of Z_layers for the largest value that lies level with or below r
                  !!! Need to subtract one because findloc returns a 1-indexed value
                  k_prev = findloc( Z_layers .le. r , .TRUE. , Back=.True., Dim=1 ) - 1  
                  
              end if
              
              interp_values(i,j) = Array(i, j, k_prev)
      
          end do
      end do
      
      end subroutine piecewise_constant_interpolation
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      subroutine piecewise_constant_integration(Array, h, D, integral_values)
      
      !Integrate assuming a piecewise-constant form for the integrand
  
      implicit none
      real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: Array
      real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: h, D
      real(wp), dimension(0:nx-1,0:ny-1), intent(out) :: integral_values
      real(wp), dimension(0:nz-1) :: z_vals
      real(wp) :: r, z_prev, delta, Upper_int
      ! real(wp) :: start_time, end_time
      integer :: i, j, p, k_prev
      
      
      do j=0, ny-1
          do i=0, nx-1
              
              !Rescale Z_layers (fraction of the way through the column) by the vertical extent of the column
              z_vals = Z_layers * D(i,j)
              
              r = h(i,j)/D(i,j)
              
              if ( r <= 0 ) then
                  
                  integral_values(i,j) = 0._wp
  
              else
                  
                  !!! Find the index of Z_layers for the largest value that lies level with or below r
                  !!! Need to subtract one because findloc returns a 1-indexed value
                  k_prev = findloc( Z_layers .le. r , .TRUE. , Back=.True., Dim=1 ) - 1  
                  z_prev = z_vals(k_prev) 
                  delta = h(i,j) - z_prev
                  
                  !Integral of values for cells above the one where h resides.
                  !This is just the sum of function values * cell heights
                  Upper_int = 0._wp
                  if (k_prev > 0) then
                      !Do I need to write this as an if?
                      do p=0,k_prev-1
                          Upper_int = Upper_int + Array(i,j,p)*( z_vals(p+1) - z_vals(p) )   
                      end do
                  end if
                  
                  integral_values(i,j) = Upper_int + Array(i,j,k_prev)*delta
                    
              end if
              
          end do
      end do
      
      end subroutine piecewise_constant_integration


end module CO2GraVISim_vertical_structure
