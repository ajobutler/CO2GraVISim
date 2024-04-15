module CO2GraVISim_input_parameters
    !This module is used to load in various input parameters that are specific to an 
    !individual run, such as reservoir topography and permeabililty arrays, so that 
    !other modules and subroutines can access them.

    use CO2GraVISim_global
    use CO2GraVISim_injection_profiles

    implicit none

    integer :: nx, ny, np
    real(wp), dimension(:,:),   allocatable :: H0, B0, D0
    real(wp), dimension(:,:),   allocatable :: perm_h, poro_h
    real(wp), dimension(:),     allocatable :: target_plot_times
    integer, dimension(0:3) :: h_BC_params, P_BC_params
    real(wp) :: dx, dy
    real(wp) :: M, Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve


    contains 

    subroutine read_reservoir_inputs(input_string)

        implicit none
        !Main Input folder
        character(len=*), intent(in)   ::  input_string

        !Input files
        character(1000)  ::  flow_params_file       
        character(1000)  ::  grid_params_file        
        character(1000)  ::  BC_params_file  
        character(1000)  ::  target_plot_times_file  
        character(1000)  ::  ceil_topo_file          
        character(1000)  ::  base_topo_file       
        character(1000)  ::  poro_h_file             
        character(1000)  ::  perm_h_file 
        character(1000)  ::  inj_locs_file 
        character(1000)  ::  inj_prof_file 

        integer :: io, j_y


        ! Generate the relevant file paths based on the specified Input Folder
        call input_folders(input_string,flow_params_file,grid_params_file,BC_params_file, target_plot_times_file, &
        & ceil_topo_file, base_topo_file, poro_h_file, perm_h_file,inj_locs_file,inj_prof_file)

        ! -- Read in Inputs ----------------------------------------------------------------------------------------------

        write(*,*) 'Reading in flow parameters'
        !Read in parameters
        open(newunit=io, file=flow_params_file, status='old', action='read')
        read(io,*) !Skip the first line, as it's formatting instructions
        read(io,*) M
        read(io,*) Gamma_val
        read(io,*) s_c_r
        read(io,*) s_a_i
        read(io,*) C_sat
        read(io,*) q_dissolve
        close(io)


        write(*,*) 'Reading in grid parameters'
        !Read in parameters
        open(newunit=io, file=grid_params_file, status='old', action='read')
        read(io,*) !Skip the first line, as it's formatting instructions
        read(io,*) nx
        read(io,*) ny
        read(io,*) dx
        read(io,*) dy
        close(io)

        write(*,*) 'Reading in boundary condition parameters'
        open(newunit=io, file=BC_params_file, status='old', action='read')
        !Skip the first 6 lines (I need a better way of dealing with comments)
        read(io,*)
        read(io,*)
        read(io,*)
        read(io,*)
        read(io,*)
        read(io,*)
        read(io,*) h_BC_params
        read(io,*) P_BC_params

        !Check if a Neumann BC has been chosen for the pressure on all four sides of the domain. 
        !This prevents ambient fluid from leaving the domain as CO2 is injected, and as a result 
        !the computation falters.
        if ( all(P_BC_params == 2) ) then
            write(*,*) ''
            write(*,*) '!!! Neumann BC chosen for the pressure on all four sides. This makes injection difficult! !!!'
            write(*,*) ''
        end if


        write(*,*) 'Reading in plot times'
        !Read in plot times
        open (newunit=io, file=target_plot_times_file, status='old', action='read')
        read(io,*)     !Skip the first line, as it's formatting instructions
        read(io, *) np !Read in the number of plot times
        allocate (target_plot_times(0:np-1))
        do j_y = 0, np-1
           read(io, *) target_plot_times(j_y)
        end do
        close (io)


        write(*,*) 'Allocating arrays'
        !Allocate array sizes
        allocate(H0(0:nx-1,0:ny-1))
        allocate(B0(0:nx-1,0:ny-1))
        allocate(perm_h(0:nx-1,0:ny-1))
        allocate(poro_h(0:nx-1,0:ny-1))



        
        !Read in input arrays
        write(*,*) 'Reading in ceiling topography array'
        open (newunit=io, file=ceil_topo_file, status='old', action='read')
        do j_y = 0, ny-1
            read (io, *) H0(:,j_y)
        end do
        close(io)

        write(*,*) 'Reading in basement topography array'
        open (newunit=io, file=base_topo_file, status='old', action='read')
        do j_y = 0, ny-1
            read (io, *) B0(:,j_y)
        end do
        close(io)

        write(*,*) 'Reading in porosity array'
        open (newunit=io, file=poro_h_file, status='old', action='read')
        do j_y = 0, ny-1
            read (io, *) poro_h(:,j_y)
        end do
        close(io)

        write(*,*) 'Reading in permeability array'
        open (newunit=io, file=perm_h_file, status='old', action='read')
        do j_y = 0, ny-1
            read (io, *) perm_h(:,j_y)
        end do
        close(io)


        ! Define the reservoir thickness from the caprock and basement topographies
        D0 = B0 - H0


        !Injection information
        call read_injection_data(inj_locs_file,inj_prof_file)


    end subroutine read_reservoir_inputs


    subroutine input_folders(input_string,flow_params_file,grid_params_file,BC_params_file, target_plot_times_file, &
        & ceil_topo_file, base_topo_file, poro_h_file, perm_h_file,inj_locs_file,inj_prof_file)
        !This subroutine generates the relevant folder paths for each of the inputs for CO2GraVISim based on a specified
        !Input Folder path in input_string.


        implicit none

        !Main Input folder
        character(len=*), intent(in)    ::  input_string
        character(1000)                 ::  input_string_b
        integer                         ::  input_string_len

        !Input files
        character(1000), intent(out)  ::  flow_params_file 
        character(1000), intent(out)  ::  grid_params_file 
        character(1000), intent(out)  ::  BC_params_file    
        character(1000), intent(out)  ::  target_plot_times_file    
        character(1000), intent(out)  ::  ceil_topo_file    
        character(1000), intent(out)  ::  base_topo_file   
        character(1000), intent(out)  ::  poro_h_file      
        character(1000), intent(out)  ::  perm_h_file   

        character(1000), intent(out)  ::  inj_locs_file
        character(1000), intent(out)  ::  inj_prof_file


        input_string_len = len_trim(input_string)

        
        if ( input_string_len==0 ) then
            !If the input is empty, use a default value
            input_string_b = './Input'
        else
            input_string_b = input_string
            !Check for a trailing '/', and remove it if present
            call trim_slashes(input_string_b)
        end if


        !Append this folder string to the front of the relevant file paths
        call combine_strings(flow_params_file       , input_string_b, "/flow_parameters.txt"    )
        call combine_strings(grid_params_file       , input_string_b, "/grid_parameters.txt"    )
        call combine_strings(BC_params_file         , input_string_b, "/boundary_conditions.txt")
        call combine_strings(target_plot_times_file , input_string_b, "/target_plot_times.txt"  )
        call combine_strings(ceil_topo_file         , input_string_b, "/ceil_topo.txt"          )
        call combine_strings(base_topo_file         , input_string_b, "/base_topo.txt"          )
        call combine_strings(poro_h_file            , input_string_b, "/porosity.txt"           )
        call combine_strings(perm_h_file            , input_string_b, "/permeability.txt"       )
        call combine_strings(inj_locs_file          , input_string_b, "/injection_locations.txt")
        call combine_strings(inj_prof_file          , input_string_b, "/injection_profile.txt"  )

    end subroutine input_folders


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


end module CO2GraVISim_input_parameters