module CO2GraVISim_input_parameters
   ! This module is used by CO2GraVISim to load in various input parameters that are 
   ! specific to an individual run, such as reservoir topography and 
   ! permeabililty arrays, so that other modules and subroutines can access them.

   use CO2GraVISim_global
   use CO2GraVISim_XML_input

   implicit none

   ! integer :: nx, ny, nz, np
   real(wp), dimension(:,:),   allocatable :: H0, B0, D0
   real(wp), dimension(:,:),   allocatable :: h_init, h_res_init, P_init
   real(wp), dimension(:,:,:), allocatable :: permeability, porosity
   real(wp) :: M, Gamma_val, q_dissolve
   real(wp) :: Length_scale, Porosity_scale, Pressure_scale, Time_scale, Flux_scale, Volume_scale
   real(wp) :: Permeability_scale_mD, Permeability_scale_m2
   real(wp), dimension(:,:,:), allocatable :: poro_CInt, perm_CInt, Zvals_array
   real(wp), dimension(:), allocatable :: Z_layers


   contains

   subroutine read_reservoir_inputs(XML_filepath, input_string, common_files_string)

      implicit none
      !Main Input folder
      character(len=*), intent(in)   ::  XML_filepath, input_string, common_files_string

      !Input files
      character(1000)  ::  ceil_topo_file
      character(1000)  ::  base_topo_file
      character(1000)  ::  porosity_file
      character(1000)  ::  permeability_file
      character(1000)  ::  h_init_file
      character(1000)  ::  h_res_init_file
      character(1000)  ::  P_init_file
      character(1000)  :: Z_layers_file

      integer :: j_z, j_p!, io, j_y, iostat    


      write(*,'(A)') 'Reading in inputs from XML file:'
      write(*,'(A)') trim(XML_filepath)
      write(*,*)
      call read_input_XML(XML_filepath)
      
  
      ! Generate the relevant file paths based on the specified Input Folder
      call input_folders(input_string, common_files_string, &
      & ceil_topo_file, base_topo_file, porosity_file, permeability_file, &
      & h_init_file, h_res_init_file, P_init_file, Z_layers_file)

      ! -- Read in Inputs ----------------------------------------------------------------------------------------------

      !Allocate array sizes
      allocate(H0(0:nx-1,0:ny-1))
      allocate(B0(0:nx-1,0:ny-1))
      allocate(D0(0:nx-1,0:ny-1))
      allocate(permeability(0:nx-1,0:ny-1,0:nz-1))
      allocate(porosity(0:nx-1,0:ny-1,0:nz-1))
      allocate(h_init(0:nx-1,0:ny-1))
      allocate(h_res_init(0:nx-1,0:ny-1))
      allocate(P_init(0:nx-1,0:ny-1))
      allocate(Z_layers(0:nz-1))

      !Read in topography arrays
      write(*,'(A)') ' - ceiling topography array'
      call get_2D_array(H0, ceil_topo_file)

      write(*,'(A)') ' - basement topography array'
      call get_2D_array(B0, base_topo_file)

      ! Define the reservoir thickness from the caprock and basement topographies
      D0 = B0 - H0



      write(*,'(A)') ' - porosity array'
      call get_3D_array(porosity, porosity_file)

      write(*,'(A)') ' - permeability array'
      call get_3D_array(permeability, permeability_file)

      write(*,'(A)') ' - initial mobile CO2 thickness (h) array'
      call get_2D_array(h_init, h_init_file)

      write(*,'(A)') ' - initial residually trapped CO2 thickness (h_res) array'
      call get_2D_array(h_res_init, h_res_init_file)

      write(*,'(A)') ' - initial non-hydrostatic ambient pressure (P) array'
      call get_2D_array(P_init, P_init_file)

      write(*,'(A)') ' - Z layers'
      call get_1D_array(Z_layers, z_layers_file)


      write(*,'(A)') ''
      write(*,'(A)') 'Calculating Cumulative Integrals'

      allocate(poro_CInt(0:nx-1,0:ny-1,0:nz-1))
      allocate(perm_CInt(0:nx-1,0:ny-1,0:nz-1))
      allocate(Zvals_array(0:nx-1,0:ny-1,0:nz-1))


      do j_z = 0,nz-1
         Zvals_array(:,:,j_z) = Z_layers(j_z) * D0(:,:)
         ! Z_layers(j_z) = real(j_z,wp)/(nz-1)
      enddo

      call cumulative_integral(Porosity    ,Zvals_array,poro_CInt)
      call cumulative_integral(Permeability,Zvals_array,perm_CInt) 

      write(*,*) '-----'
      write(*,'(A,F15.4)') 'max porosity     = ', maxval(porosity)
      write(*,'(A,F15.4)') 'max permeability = ', maxval(permeability)
      write(*,'(A,F15.4)') 'max poro_CInt    = ', maxval(poro_CInt)
      write(*,'(A,F15.4)') 'max perm_CInt    = ', maxval(perm_CInt)
      write(*,'(A,F15.4)') 'max H0           = ', maxval(H0)
      write(*,'(A,F15.4)') 'max B0           = ', maxval(B0)
      write(*,'(A,F15.4)') 'max D0           = ', maxval(D0)
      write(*,'(A,F15.4)') 'max h_init       = ', maxval(h_init)
      write(*,'(A,F15.4)') 'max h_res_init   = ', maxval(h_res_init)
      write(*,'(A,F15.4)') 'max P_init       = ', maxval(P_init)
      write(*,*) '-----'


      !!! Nondimensionalise
      call nondimensionalise_inputs(H0,B0,D0,porosity,permeability,h_init,h_res_init,P_init,dx,dy,&
      & Length_scale, Porosity_scale, Pressure_scale, Time_scale, Flux_scale, &
      & Permeability_scale_mD, Permeability_scale_m2, M, Gamma_val, q_dissolve)


      !!! Nondimensionalise lengths !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      H0 = H0 / Length_scale
      B0 = B0 / Length_scale
      D0 = D0 / Length_scale

      dx = dx / Length_scale
      dy = dy / Length_scale
  
      !!! Nondimensionalise times !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Convert dimensional times from days to seconds, then nondimensionalise
      target_plot_times = target_plot_times * days_to_seconds / Time_scale
      t_flux_vals       = t_flux_vals       * days_to_seconds / Time_scale
  
      write(*,*) ''
      write(*,'(A)') 'Nondimensional plot times:'
      do j_z = 0,np-1
          write(*,'(A,I4,A,F15.4)') '[', j_z, ']: t = ', target_plot_times(j_z)
      end do
      write(*,*) ''
  
      write(*,'(A)') 'Nondimensional flux times:'
      do j_p = 0,n_flux_times-1
          write(*,'(A,I4,A,F15.4)') '[', j_p, ']: t = ', t_flux_vals(j_p)
      end do
      write(*,*) ''
  
      !!! Nondimensional flux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Q_flux_vals = Q_flux_vals / Flux_scale
  
      !!! Nondimensionalise porosity and permeability !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      porosity     = porosity     / Porosity_scale
      permeability = permeability / Permeability_scale_mD

      poro_CInt = poro_CInt / (Porosity_scale        * Length_scale)
      perm_CInt = perm_CInt / (Permeability_scale_mD * Length_scale)
      
      Zvals_array = Zvals_array / Length_scale

      !!! Nondimensionalise initial profiles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      h_init       = h_init     / Length_scale
      h_res_init   = h_res_init / Length_scale
      P_init       = P_init     / Pressure_scale


      write(*,*) '-----'
      write(*,'(A,F15.4)') 'max porosity     = ', maxval(porosity)
      write(*,'(A,F15.4)') 'max permeability = ', maxval(permeability)
      write(*,'(A,F15.4)') 'max poro_CInt    = ', maxval(poro_CInt)
      write(*,'(A,F15.4)') 'max perm_CInt    = ', maxval(perm_CInt)
      write(*,'(A,F15.4)') 'max H0           = ', maxval(H0)
      write(*,'(A,F15.4)') 'max B0           = ', maxval(B0)
      write(*,'(A,F15.4)') 'max D0           = ', maxval(D0)
      write(*,'(A,F15.4)') 'max h_init       = ', maxval(h_init)
      write(*,'(A,F15.4)') 'max h_res_init   = ', maxval(h_res_init)
      write(*,'(A,F15.4)') 'max P_init       = ', maxval(P_init)
      write(*,*) '-----'


   end subroutine read_reservoir_inputs


   subroutine input_folders(input_string, common_files_string, &
   & ceil_topo_file, base_topo_file, porosity_file, permeability_file, &
   & h_init_file, h_res_init_file, P_init_file, Z_layers_file)
      !This subroutine generates the relevant folder paths for each of the inputs for CO2GraVISim based on a specified
      !Input Folder path in input_string.


      implicit none

      !Main Input folder
      character(len=*), intent(in)    ::  input_string, common_files_string
      character(1000)                 ::  input_string_b, common_files_string_b
      integer                         ::  input_string_len, common_files_string_len

      !Input files
      character(1000), intent(out)  ::  ceil_topo_file
      character(1000), intent(out)  ::  base_topo_file
      character(1000), intent(out)  ::  porosity_file
      character(1000), intent(out)  ::  permeability_file

      character(1000), intent(out)  ::  h_init_file, h_res_init_file, P_init_file
      character(1000), intent(out)  ::  Z_layers_file


      input_string_len = len_trim(input_string)
      common_files_string_len = len_trim(common_files_string)


      if ( input_string_len==0 ) then
         !If the input is empty, use a default value
         input_string_b = './Input'
      else
         input_string_b = input_string
         !Check for a trailing '/', and remove it if present
         call trim_slashes(input_string_b)
      end if


      if ( common_files_string_len==0 ) then
         !Use the Input folder path
         common_files_string_b = input_string_b
      else
         common_files_string_b = common_files_string
         call trim_slashes(common_files_string_b)
      end if


      write(*,'(A)') ''
      write(*,'(A)') 'Locating input data:'

    !   call specific_or_common_file(flow_params_file       , input_string_b, common_files_string_b, 'flow_parameters'      )
    !   call specific_or_common_file(grid_params_file       , input_string_b, common_files_string_b, 'grid_parameters'      )
    !   call specific_or_common_file(BC_params_file         , input_string_b, common_files_string_b, 'boundary_conditions'  )
    !   call specific_or_common_file(target_plot_times_file , input_string_b, common_files_string_b, 'target_plot_times'    )
      call specific_or_common_file(ceil_topo_file         , input_string_b, common_files_string_b, 'ceil_topo'            )
      call specific_or_common_file(base_topo_file         , input_string_b, common_files_string_b, 'base_topo'            )
      call specific_or_common_file(porosity_file          , input_string_b, common_files_string_b, 'porosity'             )
      call specific_or_common_file(permeability_file      , input_string_b, common_files_string_b, 'permeability'         )
      ! call specific_or_common_file(inj_locs_file          , input_string_b, common_files_string_b, 'injection_locations'  )
      ! call specific_or_common_file(inj_prof_file          , input_string_b, common_files_string_b, 'injection_profile'    )
      call specific_or_common_file(h_init_file            , input_string_b, common_files_string_b, 'h_init'               )
      call specific_or_common_file(h_res_init_file        , input_string_b, common_files_string_b, 'h_res_init'           )
      call specific_or_common_file(P_init_file            , input_string_b, common_files_string_b, 'P_init'               )
      call specific_or_common_file(Z_layers_file          , input_string_b, common_files_string_b, 'z_layers'               )

   end subroutine input_folders
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine specific_or_common_file(full_string, input_string, common_files_string, filename)
      !This subroutine builds the filepath for each of the input files, checking whether
      !they should come from a specific Input folder or the Common Files folder, for the case
      !of a batch run

      implicit none

      character(len=*), intent(out) :: full_string
      character(len=*), intent(in) :: input_string, common_files_string, filename
      character(1000) :: tail_string
      logical :: file_exists

      tail_string = trim( "/" // filename // ".txt" )

      !Build the filepath for the target Input folder
      call combine_strings(full_string , input_string, tail_string )

      !Check that this file exists
      inquire(file=full_string, exist=file_exists)

      if (file_exists) then
         !This file was found in the specified Input folder
         write(*,'(A,A25,A)') ' - ' , trim(filename)//repeat(' ',25) , '[Input Folder]'
      else
         !This file was not found in the specified Input folder
         !Try the Common Files folder instead
         call combine_strings(full_string , common_files_string, tail_string )
         inquire(file=full_string, exist=file_exists)
         if (file_exists) then
            !This file was found in the Common Files folder
            write(*,'(A,A25,A)') ' - ' , trim(filename)//repeat(' ',25) , '[Common Files]'
         else
            !This file was not found in the Common Files folder, either!
            write(*,'(A,A,A)') '!! ', trim(filename) , ' not found in either Input or Common Files Folders !!'
         end if

      end if

   end subroutine specific_or_common_file
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine combine_strings(full_string,base_string,tail_string)
      !This subroutine combines two strings, dealing with allocated string lengths

      implicit none

      character(len=*), intent(out) :: full_string
      character(len=*), intent(in) :: base_string, tail_string
      integer :: full_len, base_len, tail_len

      base_len = len_trim(base_string)
      tail_len = len_trim(tail_string)
      full_len = len(full_string)

      if ( base_len + tail_len > full_len ) then
         !The combined string is too long for the full_string provided
         write(*,*) 'The combined string is too long!'
         write(*,*) base_string
         write(*,*) tail_string
         stop
      end if

      full_string(:base_len+tail_len) = base_string(:base_len) // tail_string(:tail_len)

      if ( base_len+tail_len < full_len ) then
         !Fill the remainder of the string with spaces
         full_string(base_len+tail_len+1:full_len) = repeat(' ', full_len - (base_len+tail_len))
      end if

   end subroutine combine_strings
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine trim_slashes(string)
      !This subroutine removes any slashes ('\' or '/') included at the end of a
      !file path string, so that they can be in a standard format

      implicit none
      character(len=*), intent(inout) :: string
      integer :: text_len

      text_len = len_trim(string)

      !Check for a trailing '/', and remove it if present
      if (string(text_len:text_len)=='/' .or. string(text_len:text_len)=='\') then
         string(text_len:text_len) = ' '
      end if

   end subroutine trim_slashes
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine get_1D_array(Array, filepath)

      implicit none

      real(wp), dimension(0:nz-1), intent(inout) :: Array
      character(1000), intent(in) :: filepath
      integer :: io, k

      open(newunit=io, file=filepath, status='old', action='read')
      do k=0,nz-1
         read(io,*) Array(k)
      enddo
      close(io)


   end subroutine get_1D_array
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine get_2D_array(Array, filepath)

      implicit none

      real(wp), dimension(0:nx-1,0:ny-1), intent(inout) :: Array
      character(1000), intent(in) :: filepath
      integer :: io, j

      open(newunit=io, file=filepath, status='old', action='read')
      do j=0,ny-1
         read(io,*) Array(:,j)
      enddo
      close(io)


   end subroutine get_2D_array
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine get_3D_array(Array, filepath)

      implicit none

      real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(inout) :: Array
      character(1000), intent(in) :: filepath
      integer :: io, j, k

      open(newunit=io, file=filepath, status='old', action='read')
      do k=0,nz-1
         do j=0,ny-1
            read(io,*) Array(:,j,k)
         enddo
      enddo
      close(io)


   end subroutine get_3D_array
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine cumulative_integral(F_array,Z,CInt_array)
      ! This subroutines calculates the cumulative integral \int_{H}^{H+z} f(x,y,w) dw
      ! that are used for porosity and permeability, where H is the ceiling topography
      ! and z increases downwards with gravity

      implicit none
      real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: F_array, Z
      real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(out) :: CInt_array
      real(wp), dimension(0:nx-1,0:ny-1) :: Sum, dz
      integer :: k

      !k = 0
      Sum(:,:)            = 0._wp
      CInt_array(:,:,0)   = 0._wp

      do k = 1,nz-1

         dz = Z(:,:,k) - Z(:,:,k-1)

         Sum = Sum + 0.5_wp * dz(:,:) * ( F_array(:,:,k-1) + F_array(:,:,k) )
         CInt_array(:,:,k) = Sum

      enddo

   end subroutine cumulative_integral
   !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
   subroutine calculate_length_scale(Length_scale, width_vals)

    implicit none
    real(wp), dimension(0:nx-1,0:ny-1), intent(in) :: width_vals
    real(wp), intent(out) :: Length_scale
    integer :: k

    ! Set the length scale to be the mean of the reservoir widths at the injection locations
    Length_scale = 0._wp
    do k = 0,n_inj_locs-1
        Length_scale = Length_scale + width_vals( Q_inj_locs(k,0) , Q_inj_locs(k,1) )
    end do

    Length_scale = Length_scale / n_inj_locs

    end subroutine calculate_length_scale
    !----------------------------------------------------------------------------------------------


    !----------------------------------------------------------------------------------------------
    subroutine calculate_flux_scale(Flux_scale)

    implicit none
    real(wp), intent(out) :: Flux_scale
    real(wp), dimension(0:n_flux_times) :: t_vals
    real(wp), dimension(0:n_flux_times-1) :: t_interval_vals, Q_sum_vals
    real(wp) :: total_flux, total_injecting_time
    integer :: k

    ! t_vals = [t_flux_1, t_flux_2, ..., t_flux_N, t_end]
    t_vals(0:n_flux_times-1) = t_flux_vals(0:n_flux_times-1)
    t_vals(n_flux_times) = target_plot_times(Ubound(target_plot_times,1))

    ! intervals between these times
    t_interval_vals = t_vals(1:n_flux_times) - t_vals(0:n_flux_times-1)

    ! Sum the flux values in each interval over all the wells
    Q_sum_vals = sum( Q_flux_vals, dim=2)

    ! Total time that any well is injecting
    total_injecting_time = sum( t_interval_vals , Q_sum_vals > 0._wp )

    if (total_injecting_time .le. 0._wp) then
      ! Total injecting time is zero - likely due to simulating an initial profile without injection.
      ! Let's pick a representative value for the Flux scale then
      Flux_scale = 1e0_wp !m^2 s^-1
      write(*,'(A)') '!!! Total injection time is zero !!!'
      write(*,'(A, G15.4)') 'Choosing a nominal value of ', Flux_scale

    else 
      ! Total flux (easy integration, since piecewise constant)
      total_flux = sum( Q_sum_vals * t_interval_vals )

      ! Flux scale as the average injected flux per well
      Flux_scale = abs(total_flux) / (n_inj_locs * total_injecting_time)
    end if

    end subroutine calculate_flux_scale
    !----------------------------------------------------------------------------------------------


    !----------------------------------------------------------------------------------------------
    subroutine calculate_perm_scale(Perm_scale_mD, Perm_scale_m2, Perm_array) !, ActNum_array)

    implicit none
    real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: Perm_array
    ! integer, dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: ActNum_array
    real(wp), intent(out) :: Perm_scale_mD, Perm_scale_m2
   !  real(wp) :: max_val, min_val

    ! Arithmetic mean
    Perm_scale_mD = sum( Perm_array ) / (nx*ny*nz) !size( Perm_array )

    !Convert from mD to m2
    Perm_scale_m2 = Perm_scale_mD * mD_to_m2

    end subroutine calculate_perm_scale
    !----------------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------------
    subroutine nondimensionalise_inputs(H0,B0,D0,porosity,permeability,h_init,h_res_init,P_init,dx,dy,&
        & Length_scale, Porosity_scale, Pressure_scale, Time_scale, Flux_scale, &
        & Permeability_scale_mD, Permeability_scale_m2, M, Gamma_val, q_dissolve)

        implicit none
        real(wp), dimension(0:nx-1,0:ny-1), intent(inout) :: H0, B0, D0, h_init, h_res_init, P_init
        real(wp), dimension(0:nx-1,0:ny-1,0:nz-1), intent(inout) :: porosity, permeability
        real(wp), intent(inout) :: dx, dy
        real(wp), intent(out) :: Length_scale, Porosity_scale, Pressure_scale, Time_scale, Flux_scale
        real(wp), intent(out) :: Permeability_scale_mD, Permeability_scale_m2
        real(wp), intent(out) :: M, Gamma_val, q_dissolve
        real(wp) :: u_b, u_Q
        integer :: k

        ! Calculate the relevant dimensional scales used to construct the nondimensional parameters
        ! and variables below

        call calculate_length_scale(Length_scale, D0)

        call calculate_flux_scale(Flux_scale)

        call calculate_perm_scale(Permeability_scale_mD, Permeability_scale_m2, Permeability)

        Porosity_scale = 1._wp 

        Volume_scale = Porosity_scale * Length_scale**3


        ! Calculate the remaining dimensional scales and nondimensional parameters

        !Buoyancy velocity scale
        u_b = ( (rho_a - rho_c) * g * Permeability_scale_m2 ) / mu_c
        !Injection flux velocity scale
        u_Q = Flux_scale / ( Length_scale**2 )
        !Advective time scale based on injection velocity
        Time_scale = ( Porosity_scale * Length_scale ) / u_Q
        !Pressure scale based in injection velocity
        Pressure_scale = ( mu_c / Permeability_scale_m2 ) * ( Flux_scale / Length_scale )
        !Viscosity ratio between CO2 and ambient
        M = mu_c / mu_a
        !Buoyancy-Injection velocity ratio
        Gamma_val = u_b / u_Q

        !Dimensional convective dissolution parameter
        q_dissolve = 0.12_wp * ( Porosity_scale * C_sat * D_mol / Length_scale ) &
        & * ( ( (rho_a_sat - rho_a) * g * perm_ratio * Permeability_scale_m2 * Length_scale )&
        &     / ( Porosity_scale * mu_a * D_mol ) )**(0.84) !Dimensional version

        q_dissolve = omega_conv * q_dissolve !Incorporate tuning parameter

        q_dissolve = q_dissolve / u_Q !Nondimensional version


        write(*,*)
        write(*,'(A)') '---- Dimensional Scales and Nondimensional Parameters ----'
        write(*,*)
        write(*,'(A,G15.4,A)') 'Time scale         = ', Time_scale,            ' [s]'
        write(*,'(A,G15.4,A)') 'Length scale       = ', Length_scale,          ' [m]'
        write(*,'(A,G15.4,A)') 'Flux scale         = ', Flux_scale,            ' [m^3 s^-1]'
        write(*,'(A,G15.4,A)') 'Pressure scale     = ', Pressure_scale,        ' [Pa]'
        write(*,'(A,G15.4,A)') 'Porosity scale     = ', Porosity_scale,        ' [-]'
        write(*,'(A,G15.4,A)') 'Permeability scale = ', Permeability_scale_m2, ' [m^2]'
        write(*,'(A,G15.4,A)') '                   = ', Permeability_scale_mD, ' [mD]'
        write(*,'(A,G15.4,A)') ''
        write(*,'(A,G15.4,A)') 'Buoyancy velocity  = ', u_b,                   ' [m s^-1]'
        write(*,'(A,G15.4,A)') 'Flux velocity      = ', u_Q,                   ' [m s^-1]'
        write(*,'(A,G15.4,A)') ''
        write(*,'(A,G15.4,A)') 'M                  = ', M,          ''
        write(*,'(A,G15.4,A)') 'Gamma_val          = ', Gamma_val,  ''
        write(*,'(A,G15.4,A)') 'q_dissolve         = ', q_dissolve, ''


    end subroutine nondimensionalise_inputs
    !----------------------------------------------------------------------------------------------


end module CO2GraVISim_input_parameters
