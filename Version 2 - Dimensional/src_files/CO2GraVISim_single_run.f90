program CO2GraVISim_single_run
   ! This program performs a single run of CO2GraVISim.
   ! Various input parameters and arrays are read into memory, the main solver program
   ! is called, and then the corresponding output files are written.

   use CO2GraVISim_global !wp, iter_max, tmax, dt_init, dt_min, errmax_t_step, errmax_P_iter, P_iter_omega, t0
   use CO2GraVISim_input_parameters  !H0, B0, D0, perm_h, poro_h, n_inj_locs, n_flux_times, Q_inj_locs, Q_flux_vals, t_flux_vals
   use CO2GraVISim_solver, only: solver_run
   ! use CO2GraVISim_solver_unconfined, only: unconfined_solver_run
   use CO2GraVISim_logging
   use omp_lib

   implicit none
   integer :: j_x, j_y, j_p, io, i_arg, arg_count
   real(wp), dimension(:,:,:), allocatable :: h_array, h_res_array, P_array, V_mobile_array, V_trapped_array
   real(wp), dimension(:,:),   allocatable :: Max_mobile_thickness, Max_pressure, Min_pressure
   real(wp), dimension(:),     allocatable :: plot_times
   real(wp) :: start_save, finish_save

   !Command-Line inputs
   character(1000)  :: input_folder, output_folder, batch_folder, common_files_folder, XML_filepath
   character(1000)  :: arg

   !Output files
   character(1000)  ::  output_times_file
   character(1000)  ::  param_output_file
   character(1000)  ::  dimensional_scales_file
   character(1000)  ::  inj_locs_output_file
   character(1000)  ::  Max_mobile_thickness_file
   character(1000)  ::  Max_pressure_file
   character(1000)  ::  Min_pressure_file

   character(30) :: Presence_row_format

   logical :: CL_Confined_flag, CL_Unconfined_flag
   logical :: Confined_flag
   logical :: XML_exist


   ! -- Read in Command Line inputs --------------------------------------------------------------------------------

   !Default values
   input_folder  = './Input'
   output_folder = './Output'
   XML_filepath  = ''
   batch_folder = ''
   common_files_folder = './Input' !If not specified via -batch, reuse the main Input folder
   CL_Confined_flag     = .False.
   CL_Unconfined_flag   = .False.
   Confined_flag        = .True.


   !Parse Command Line inputs - based on code by Jon D'Souza-Eva
   arg_count = command_argument_count()
   i_arg = 1
   do while( i_arg <= arg_count )
      call get_command_argument(i_arg, arg)

      select case(arg)

      case('-input')
         if (i_arg < arg_count) then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, input_folder)
            call trim_slashes(input_folder)
            if (XML_filepath == '' ) then
               !No XML filepath specified - set to the first one we can find in the Input folder
               XML_filepath = trim(input_folder) // '/*.xml'
            end if
         endif

      case('-output')
         if (i_arg < arg_count) then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, output_folder)
            call trim_slashes(output_folder)
         endif

      case('-batch')
         if (i_arg < arg_count) then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, batch_folder)
            call trim_slashes(batch_folder)
            !Set input and output folders relative to this
            input_folder = trim(batch_folder) // '/Input'
            output_folder = trim(batch_folder) // '/Output'
            common_files_folder = trim(batch_folder) // '/../Common_files'
            ! XML_filepath = trim(input_folder) // '/*.xml' 
            !Can't deal with wildcards here, they need to be resolved in the command line
         endif

      case('-xml')
         if (i_arg < arg_count) then
            i_arg = i_arg + 1
            call get_command_argument(i_arg, XML_filepath)
         endif

      case('-confined')
         !Set the command-line Confined_flag so that the unconfined options in the solver are used
         CL_Confined_flag = .True.
      
      case('-unconfined')
         !Set Confined_flag so that the unconfined options in the solver are used
         CL_Unconfined_flag = .True.
         ! Confined_flag = .False.

      case default
         !Unknown command-line argument encountered. List valid options.
         write(*,'(A,A)') 'Unknown argument: ', trim(arg)
         write(*,'(A)') ''
         write(*,'(A)') 'Possible arguments are:'

         write(*,'(A)') ' -xml "/Path to/Input XML File.xml" : &
                     &  Specify the path to the input XML file for this simulation'

         write(*,'(A)') ' -input "/Path to/target input folder/" : &
         & Specify the path to the target Input folder'
         
         write(*,'(A)') ' -output "/Path to/target output folder/" : &
         & Specify the path to the target Output folder'
         
         write(*,'(A)') ' -confined : &
         & Set the simulator to run in Confined mode'
         
         write(*,'(A)') ' -unconfined : &
         & Set the simulator to run in Unconfined mode'
         
         write(*,'(A)') ' -batch "Path to/ batch run folder n/" : &
         & Specify the main folder for a particular run in a batch. &
         & This run folder sits level with the Common Files folder, and contains Input and Output folders'
         
         write(*,'(A)') ''
         stop

      end select

      i_arg = i_arg + 1

   end do

   !Remove any trailing '/' or '\' that might have been included
   ! call trim_slashes(input_folder)
   ! call trim_slashes(common_files_folder)
   call trim_slashes(output_folder)

   ! write(*,'(A)') repeat('-',37)
   call log_input_line('-', '')
   call log_input_line('-', 'Filepaths')
   call log_input_line('-', '')
   write(*,'(3A)') 'XML input file:      [' , trim(XML_filepath)  , ']'
   write(*,'(3A)') 'Input Folder:        [' , trim(input_folder)  , ']'
   write(*,'(3A)') 'Common Files Folder: [' , trim(common_files_folder)  , ']'
   write(*,'(3A)') 'Output Folder:       [' , trim(output_folder) , ']'
   ! write(*,'(A)') repeat('-',37)
   call log_input_line('-', '')
   call log_input_line(' ', '')
   call log_input_line(' ', '')


   ! -- Read in Inputs ----------------------------------------------------------------------------------------------

   !Check that the target XML file exists
   inquire(file=trim(XML_filepath), exist=XML_exist)
   if (.not. XML_exist) then
       write(*,'(A)') '!! Unable to access XML input file !!'
       write(*,'(A)') trim(XML_filepath)
       write(*,*) ''
       stop
   end if

   call read_reservoir_inputs(XML_filepath,input_folder,common_files_folder) !from the -input_parameters- module

   ! -- Set Simulator mode (Confined/Unconfined) --------------------------------------------------------------------
    
   call set_simulator_mode(Confined_flag,CL_confined_flag,CL_unconfined_flag,XML_sim_mode)

   !Set some temporary values in the array of actual plot_times.
   !The final values will be close to the target_plot_times values, but
   !may differ slightly due to the time-stepping routine used
   allocate(plot_times(0:np-1))
   plot_times = target_plot_times

   write(*,'(A)') ''
   write(*,'(A)') 'Allocating output arrays'
   !Allocated array sizes
   allocate(h_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(h_res_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(P_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(V_mobile_array(0:nx-1, 0:ny-1, 0:np-1))
   allocate(V_trapped_array(0:nx-1, 0:ny-1, 0:np-1))

   allocate(Max_mobile_thickness(0:nx-1, 0:ny-1))
   allocate(Max_pressure(0:nx-1, 0:ny-1))
   allocate(Min_pressure(0:nx-1, 0:ny-1))


   ! -- Print Permeability and Porosity at the Injection Locations ---------------------------------------------------

   write(*,'(A)') 'Injection locations - properties in z:'
   write(*,'(A)') repeat('-', 5+5+8+25+25 + 4*3 + 2*2)
   write(*,'("| ", A5," | ",A5," | ",A13," | ",A13," | ",A15," | ",A15, " |")') &
   &                          'idx_x', 'idx_y', 'Min. porosity', 'Max. porosity', 'Min. perm. [mD]', 'Max. perm. [mD]'
   write(*,'(A)') repeat('-', 5+5+8+25+25 + 4*3 + 2*2)
   do j_x = 0, n_inj_locs-1
       write(*,'("| ", I5," | ", I5," | ", F13.3," | ", F13.3, " | ", F15.8," | ",F15.8, " |")') &
           &  Q_inj_locs(j_x,0) + 1 , Q_inj_locs(j_x,1) + 1 ,&
           &   minval(     Porosity( Q_inj_locs(j_x,0) , Q_inj_locs(j_x,1) , : ) ) * Porosity_scale        ,&
           &   maxval(     Porosity( Q_inj_locs(j_x,0) , Q_inj_locs(j_x,1) , : ) ) * Porosity_scale        ,&
           &   minval( Permeability( Q_inj_locs(j_x,0) , Q_inj_locs(j_x,1) , : ) ) * Permeability_scale_mD ,&
           &   maxval( Permeability( Q_inj_locs(j_x,0) , Q_inj_locs(j_x,1) , : ) ) * Permeability_scale_mD
   end do
   write(*,'(A)') repeat('-', 5+5+8+25+25 + 4*3 + 2*2)


   ! -- Run Solver ---------------------------------------------------------------------------------------------------

   call solver_run(h_array,h_res_array,P_array,V_mobile_array,V_trapped_array,&
   & Max_mobile_thickness,Max_pressure,Min_pressure,plot_times,output_folder,Confined_flag)

   write(*,'(A)') ''
   write(*,*) '-- Solver run complete --'
   write(*,'(A)') ''


   ! -- Redimensionalise Outputs -------------------------------------------------------------------------------------

   h_array              = h_array               * Length_scale
   h_res_array          = h_res_array           * Length_scale
   Max_mobile_thickness = Max_mobile_thickness  * Length_scale
   
   V_mobile_array       = V_mobile_array  * Volume_scale
   V_trapped_array      = V_trapped_array * Volume_scale

   P_array              = P_array      * Pressure_scale
   Max_pressure         = Max_pressure * Pressure_scale
   Min_pressure         = Min_pressure * Pressure_scale

   plot_times           = plot_times   * Time_scale * seconds_to_days !Express in days, rather than seconds


   ! -- Save Outputs -------------------------------------------------------------------------------------------------

   start_save = omp_get_wtime()

   write(*,'(A)') ''
   write(*,'(A)') 'Saving outputs:'

   !Save arrays that contain a different slice for each plot time
   write(*,'(A)') ' - mobile current thickness profiles'
   call save_profile_array(h_array, "/Current_Thickness/h")

   write(*,'(A)') ' - mobile current volume profiles'
   call save_profile_array(V_mobile_array, "/Current_Volume/V")

   write(*,'(A)') ' - trapped current thickness profiles'
   call save_profile_array(h_res_array, "/Current_Thickness/h_res")

   write(*,'(A)') ' - trapped current volume profiles'
   call save_profile_array(V_trapped_array, "/Current_Volume/V_res")

   write(*,'(A)') ' - ambient pressure profiles'
   call save_profile_array(P_array, "/Current_Pressure/P")


   write(*,'(A)') ' - maximum mobile thickness array'
   write(Max_mobile_thickness_file, '( A , "/Other/Max_mobile_thickness.txt")') trim(output_folder)
   open (newunit=io, file=Max_mobile_thickness_file, action='write')
   write(Presence_row_format, '(A, I0, A)') '(', nx, '(E16.8E3,1X))'
   do j_y = 0, ny-1
      write(io,Presence_row_format) Max_mobile_thickness(:,j_y)
   end do
   close(io)


   write(*,'(A)') ' - maximum pressure array'
   write(Max_pressure_file, '( A , "/Other/Max_pressure.txt")') trim(output_folder)
   open (newunit=io, file=Max_pressure_file, action='write')
   do j_y = 0, ny-1
      write(io,Presence_row_format) Max_pressure(:,j_y)
   end do
   close(io)


   write(*,'(A)') ' - minimum pressure array'
   write(Min_pressure_file, '( A , "/Other/Min_pressure.txt")') trim(output_folder)
   open (newunit=io, file=Min_pressure_file, action='write')
   do j_y = 0, ny-1
      write(io,Presence_row_format) Min_pressure(:,j_y)
   end do
   close(io)


   write(*,'(A)') ' - plot times'
   write(output_times_file, '( A , "/Other/plot_times.txt")') trim(output_folder)
   open (newunit=io, file=output_times_file, action='write')
   do j_p = 0, np-1
      write(io, *) plot_times(j_p)
   end do
   close(io)


   write(*,'(A)') ' - run parameters'
   write(param_output_file, '( A , "/Other/parameter_values.txt")') trim(output_folder)
   open (newunit=io, file=param_output_file, action='write')
   write(io, '(I5)') nx
   write(io, '(I5)') ny
   write(io, '(I5)') nz
   write(io, '(F10.3)') dx * Length_scale
   write(io, '(F10.3)') dy * Length_scale
   write(io, '(F10.3)') M
   write(io, '(F10.3)') Gamma_val
   write(io, '(F10.3)') s_c_r
   write(io, '(F10.3)') s_a_i
   write(io, '(F10.3)') C_sat
   write(io, '(F10.3)') q_dissolve
   close (io)


   write(*,'(A)') ' - dimensional scales'
   write(dimensional_scales_file, '( A , "/Other/dimensional_scales.txt")') trim(output_folder)
   open (newunit=io, file=dimensional_scales_file, action='write')
   write(io, '(A)') '-- Length scale [m]'
   write(io, '(E16.8E3)') Length_scale
   write(io, '(A)') '-- Time scale [s]'
   write(io, '(E16.8E3)') Time_scale
   write(io, '(A)') '-- Time scale [days]'
   write(io, '(E16.8E3)') Time_scale * seconds_to_days
   write(io, '(A)') '-- Flux scale [m^3 / s]'
   write(io, '(E16.8E3)') Flux_scale
   write(io, '(A)') '-- Pressure scale [Pa]'
   write(io, '(E16.8E3)') Pressure_scale
   write(io, '(A)') '-- Porosity scale [-]'
   write(io, '(E16.8E3)') Porosity_scale
   write(io, '(A)') '-- Permeability scale [milliDarcy]'
   write(io, '(E16.8E3)') Permeability_scale_mD
   write(io, '(A)') '-- Permeability scale [m^2]'
   write(io, '(E16.8E3)') Permeability_scale_m2

   close (io)


   write(*,'(A)') ' - injection locations'
   write(inj_locs_output_file, '( A , "/Other/injection_locations.txt")') trim(output_folder)
   open (newunit=io, file=inj_locs_output_file, action='write')
   do j_p = 0, n_inj_locs-1
      write(io, *) Q_inj_locs(j_p,0), Q_inj_locs(j_p,1)
   end do
   close(io)


   finish_save = omp_get_wtime()

   write(*,'(A,G0.3,A)') 'Time taken to save data = ', finish_save - start_save, ' seconds.' 


   write(*,'(A)') ''
   write(*,*) 'Run complete.'
   write(*,'(A)') ''


   ! Deallocate arrays now that they are no longer needed
   ! (this is mainly to placate Valgrind, but might have other benefits)
   deallocate(plot_times)
   deallocate(h_array)
   deallocate(h_res_array)
   deallocate(P_array)
   deallocate(V_mobile_array)
   deallocate(V_trapped_array)
   deallocate(Max_mobile_thickness)
   deallocate(Max_pressure)
   deallocate(Min_pressure)

   deallocate(Q_inj_locs)
   deallocate(Q_flux_vals)
   deallocate(t_flux_vals)

   deallocate(H0)
   deallocate(B0)
   deallocate(D0)
   ! deallocate(perm_h)
   ! deallocate(poro_h)
   deallocate(permeability)
   deallocate(porosity)
   deallocate(target_plot_times)
   deallocate(poro_CInt)
   deallocate(perm_CInt)
   deallocate(Zvals_array)
   deallocate(Z_layers)
   
contains

! ---------------------------------------------------------------------------------------------   
subroutine save_profile_array(Array,target_folder)

      implicit none
      real(wp), dimension(0:nx-1,0:ny-1,0:np-1), intent(in) :: Array
      character(len=*), intent(in) :: target_folder
      character(1000) :: f_output_file, row_format
      real(wp), dimension(0:nx-1,0:ny-1) :: A_p

      write(row_format,'(A,I0,A)') '(', nx, '(E16.8E3, 1X))'

      ! !$OMP PARALLEL SHARED(output_folder,target_folder,Array,ny,np) PRIVATE(A_p,j_p,j_y,f_output_file,io) NUM_THREADS(4)
      ! !$OMP DO SCHEDULE(DYNAMIC)
      do j_p = 0, np-1
         
         A_p(:,:) = Array(:,:,j_p)

         write(f_output_file, '( A , A, I2.2, ".txt")') trim(output_folder), trim(target_folder), j_p

         open (newunit=io, file=f_output_file, action='write')
         do j_y = 0, ny-1
            ! write(io, *) A_p(:,j_y)
            write(io, row_format) A_p(:,j_y)
         end do
         close(io)

      end do
      ! !$OMP END DO
      ! !$OMP END PARALLEL

   end subroutine save_profile_array
   ! ---------------------------------------------------------------------------------------------


   ! ---------------------------------------------------------------------------------------------
   subroutine log_input_line(line_character, text_string)
        
      implicit none
      character(1), intent(in) :: line_character
      character(*), intent(in) :: text_string
      integer :: total_length, text_length, left_length, right_length
      

      total_length = 70
      text_length = len_trim(text_string)

      if (text_length .eq. 0) then
          !Blank title, used to print a line of hyphens with no gaps
          write(*,'(A)') repeat(line_character, total_length)
      else
          ! Put half (after rounding) of the hyphens on the left, once we've accounted for the title string
          ! and the padding spaces
          left_length = ( total_length - (text_length+2) )/2 !Integer division should round as I want
          !Put the remainder of the hyphens on the right
          right_length = total_length - text_length - left_length - 2
  
          !Print this full string to STDOUT
          write(*,'(A,1X,A,1X,A)') repeat(line_character, left_length) , trim(text_string) , repeat(line_character, right_length) 
          
      end if

   end subroutine log_input_line
   ! ---------------------------------------------------------------------------------------------


   ! ---------------------------------------------------------------------------------------------
   subroutine set_simulator_mode(Confined_flag, CL_confined_flag, CL_unconfined_flag, XML_sim_mode)
    
      implicit none
      logical, intent(out) :: Confined_flag
      logical, intent(in) :: CL_confined_flag, CL_unconfined_flag
      character(1), intent(in) :: XML_sim_mode
      
      if (CL_confined_flag .and. CL_unconfined_flag) then
          !Both -confined and -unconfined specified at the command line!
          write(*,'(A)') ''
          write(*,'(A)') '!!! Both -confined and -unconfined have been specified at the command line! !!!'
          write(*,'(A)') 'Please choose one or the other.'
          write(*,'(A)') ''
          stop
          
      end if
      
          
      !Now either -confined has been set, -unconfined has been set, or neither.
      !   1) If -confined, I want Confined_flag = .True.
      !   2) If -unconfined, I want Confined_flag = .False.
      !   3) If neither, I want to set Confined_flag based on XML_sim_mode
          
      if (CL_confined_flag) then
          !Confined mode requested at command line - this takes precedence over any XML input value
          Confined_flag = .True.
          
          write(*,'(A)') ''
          write(*,'(A)') ' ========================================================='
          write(*,'(A)') ' | -- Simulator Mode: Confined [set via command line] -- |'
          write(*,'(A)') ' ========================================================='
          write(*,'(A)') ''
          
      else if (CL_unconfined_flag) then
          !Unconfined mode requested at command line - this takes precedence over any XML input value
          Confined_flag = .False.
          
          write(*,'(A)') ''
          write(*,'(A)') ' ==========================================================='
          write(*,'(A)') ' | -- Simulator Mode: Unconfined [set via command line] -- |'
          write(*,'(A)') ' ==========================================================='
          write(*,'(A)') ''
          
      else
          
          !Mode not set at the command line - now defer to the XML input value
              
          select case (XML_sim_mode)
                  
          case ('C')
              !XML input requested Confined
              Confined_flag = .True.
              
              write(*,'(A)') ''
              write(*,'(A)') ' ======================================================'
              write(*,'(A)') ' | -- Simulator Mode: Confined [set via XML Input] -- |'
              write(*,'(A)') ' ======================================================'   
              write(*,'(A)') ''
              
          case ('U')
              !XML input requested Unconfined
              Confined_flag = .False.
              
              write(*,'(A)') ''
              write(*,'(A)') ' ========================================================'
              write(*,'(A)') ' | -- Simulator Mode: Unconfined [set via XML Input] -- |'
              write(*,'(A)') ' ========================================================' 
              write(*,'(A)') ''
              
          case ('')
              !Nothing requested from XML input - default to Unconfined
              Confined_flag = .False.
              
              write(*,'(A)') ''
              write(*,'(A)') ' ==================================================================='
              write(*,'(A)') ' | -- Simulator Mode: Unconfined [not specified; using default] -- |'
              write(*,'(A)') ' ===================================================================' 
              write(*,'(A)') ''
                  
          end select
          
          
      end if
          
      
      end subroutine set_simulator_mode
      ! ---------------------------------------------------------------------------------------------
      

end program CO2GraVISim_single_run
