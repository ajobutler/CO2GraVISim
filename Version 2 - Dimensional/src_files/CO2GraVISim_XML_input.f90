    ! CO2GraVISim_XML_input.f90

    module CO2GraVISim_XML_input

    use CO2GraVISim_global
    use FoX_dom


    implicit none

    !General FoX variables
    type(Node), pointer :: myDoc
    type(DOMException) :: Exptn
    integer :: i, j

    ! !Dataset Filepath
    ! character(1000) :: dataset_folder

    !Simulation Mode
    character(1) :: XML_sim_mode

    !Dimensional flow parameters
    real(wp) :: rho_c, rho_a, rho_a_sat, mu_c, mu_a, g, D_mol

    !Nondimensional flow parameters
    real(wp) :: s_c_r, s_a_i, C_sat, perm_ratio, omega_conv
    real(wp) :: krn_mobile, krw_residual


    !Grid parameters
    integer :: nx, ny, nz
    real(wp) :: dx, dy

    ! !Injection interval parameters
    ! integer :: z_idx_ceil, z_idx_base

    !Boundary conditions
    integer, dimension(0:3) :: h_BC_params, P_BC_params

    !Injection locations
    integer :: n_inj_locs
    integer, dimension(:,:), allocatable :: Q_inj_locs
    character(10), dimension(:), allocatable :: inj_loc_tags

    !Injection schedule
    integer :: n_flux_times
    real(wp), dimension(:), allocatable :: t_flux_vals
    real(wp), dimension(:,:), allocatable :: Q_flux_vals

    !Plot times
    integer :: np
    real(wp), dimension(:), allocatable :: target_plot_times

    !Simulation parameters
    integer :: itermax_I, P_iter_check_I
    real(wp) :: tmax_I, dt_init_I, dt_min_I, errmax_t_step_I, errmax_P_iter_I, P_iter_omega_I

    character(1000) :: row_fmt

    Private

    Public :: rho_c, rho_a, rho_a_sat, mu_c, mu_a, g, D_mol, s_c_r, s_a_i, krn_mobile, krw_residual, &
        & C_sat, perm_ratio, omega_conv, &
        & nx, ny, nz, dx, dy, &
        & h_BC_params, P_BC_params, &
        & n_inj_locs, Q_inj_locs, inj_loc_tags, n_flux_times, t_flux_vals, Q_flux_vals, &
        & np, target_plot_times, &
        & XML_sim_mode, &
        & read_input_XML



    contains


    subroutine read_input_XML(XML_filepath)

    implicit none
    character(1000), intent(in) :: XML_filepath

    ! Load in the XML document
    myDoc => parseFile( trim(XML_filepath), ex = Exptn )
    call Dom_error_check(Exptn)


    ! Header for log of information read from the XML input file
    call section_title('')
    call section_title('Input Values from XML File')
    call section_title('')
    write(*,'(A)') ''
    write(*,'(A)') ''


    ! ! Dataset Filepath
    ! call section_title('Dataset Filepath')
    ! call get_filepath(myDoc,dataset_folder)
    ! write(*,'(A)') '- filepath: '
    ! write(*,'(A)') trim(dataset_folder)

    ! Simulation Mode
    call section_title('Simulation Mode')
    call get_simulation_mode(myDoc,XML_sim_mode)
    
    ! Grid parameters
    call section_title('Grid parameters')
    call get_grid_parameters(myDoc,nx,ny,nz,dx,dy)


    ! Boundary Conditions
    call section_title('Boundary Conditions')
    call get_boundary_conditions(myDoc,h_BC_params,P_BC_params)


    ! Flow parameters
    call section_title('Dimensional Flow parameters')
    call get_dim_flow_parameters(myDoc,rho_c,rho_a,rho_a_sat,mu_c,mu_a,s_c_r,s_a_i,krn_mobile,krw_residual,&
        & C_sat,g,D_mol,perm_ratio,omega_conv)
    

    ! ! Injection Interval parameters
    ! call section_title('Injection Interval parameters')
    ! call get_injection_interval_parameters(myDoc,z_idx_ceil,z_idx_base)



    ! Injection Locations
    call section_title('Injection Locations')
    call get_inj_locs(myDoc,n_inj_locs,inj_loc_tags,Q_inj_locs,nx,ny)


    !Injection Schedule
    call section_title('Injection Schedule')
    call get_flux_info(myDoc,n_inj_locs,inj_loc_tags,n_flux_times,t_flux_vals,Q_flux_vals)


    !Plot times
    call section_title('Plot Times')
    call get_plot_times(myDoc,np,target_plot_times)


    !Simulation parameters
    call section_title('Simulation Parameters')
    call get_simulation_parameters(myDoc,itermax_I, tmax_I, dt_init_I, dt_min_I, errmax_t_step_I, &
        & errmax_P_iter_I, P_iter_omega_I, P_iter_check_I)


    call section_title('')
    call section_title('')

    ! Clear up all allocated memory
    call destroy(myDoc)


    end subroutine read_input_XML
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine get_grid_parameters(doc,nx,ny,nz,dx,dy)
        ![get_grid_parameters] looks for the <grid_parameters> block and extracts the grid parameters used
        !in the simulation

    implicit none

    type(node), pointer, intent(in) :: doc
    type(node), pointer :: p_0, p_1
    type(nodelist), pointer :: parameterlist, c_0
    integer, intent(out) :: nx, ny, nz
    real(wp), intent(out) :: dx, dy
    integer :: i, j, iostat
    character(1000) :: node_path
    type(node), pointer :: node_pointers(:)
    logical :: nx_check, ny_check, nz_check, dx_check, dy_check

    nx_check = .False.
    ny_check = .False.
    nz_check = .False.
    dx_check = .False.
    dy_check = .False.

    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "grid_parameters")

    if ( getlength(parameterlist) == 0 ) then
        write(*,'(A)') '<grid_parameters> not specified!'
        stop
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,'(A)') 'Too many <grid_parameters> blocks specified!'
        stop
    endif

    !look at the sub-blocks of the first <grid_parameters> block, and extract the relevant values
    p_0 => item(parameterlist, 0)
    c_0 => getchildnodes(p_0)

    allocate(node_pointers(2))

    do i = 0, getlength(c_0)-1

        p_1 => item(c_0, i)

        node_pointers = [p_0, p_1]
        node_path = build_node_path(node_pointers)

        select case ( getlocalname(p_1) )

        case ('nx')
            call extractdatacontent(p_1, nx, iostat=iostat)
            nx_check = .True.
        case ('ny')
            call extractdatacontent(p_1, ny, iostat=iostat)
            ny_check = .True.
        case ('nz')
            call extractdatacontent(p_1, nz, iostat=iostat)
            nz_check = .True.
        case ('dx')
            call extractdatacontent(p_1, dx, iostat=iostat)
            dx_check = .True.
        case ('dy')
            call extractdatacontent(p_1, dy, iostat=iostat)
            dy_check = .True.
        case default 
            iostat = 0
        end select

        call error_check(iostat, node_path)

    end do

    !Check the presence and validity of these values

    !-- nx -------------------------------------------------------------------------------------------
    if(.not. nx_check) then
        ! nx hasn't be set
        write(*,'(A)') "!! A value hasn't been specified for nx! !!"
        stop
    else
        if( nx < 0 ) then
            !The specified value is invalid
            write(*,'(A,I0,A)') "!! The value specified for nx is invalid : ", nx, " - need (nx>0) !!"
            stop
        else
            !The specified value is valid
            write(*,'(A)') '- nx: '
            write(*,'(I0)') nx
        endif
    endif


    !-- ny -------------------------------------------------------------------------------------------
    if(.not. ny_check) then
        ! ny hasn't be set
        write(*,'(A)') "!! A value hasn't been specified for ny! !!"
        stop
    else
        if( ny < 0 ) then
            !The specified value is invalid
            write(*,'(A,I0,A)') "!! The value specified for ny is invalid : ", ny, " - need (ny>0) !!"
            stop
        else
            !The specified value is valid
            write(*,'(A)') '- ny: '
            write(*,'(I0)') ny
        endif
    endif


    !-- nz -------------------------------------------------------------------------------------------
    if(.not. nz_check) then
        ! nz hasn't be set
        write(*,'(A)') "!! A value hasn't been specified for nz! !!"
        stop
    else
        if( nz < 0 ) then
            !The specified value is invalid
            write(*,'(A,I0,A)') "!! The value specified for nz is invalid : ", nz, " - need (nz>0) !!"
            stop
        else
            !The specified value is valid
            write(*,'(A)') '- nz: '
            write(*,'(I0)') nz
        endif
    endif


    !-- dx -------------------------------------------------------------------------------------------
    if(.not. dx_check) then
        ! dx hasn't be set
        write(*,'(A)') "!! A value hasn't been specified for dx! !!"
        stop
    else
        if( dx < 0 ) then
            !The specified value is invalid
            write(*,'(A,F8.3,A)') "!! The value specified for dx is invalid : ", dx, " - need (dx>0) !!"
            stop
        else
            !The specified value is valid
            write(*,'(A)') '- dx: '
            write(*,'(F8.3)') dx
        endif
    endif


    !-- dy -------------------------------------------------------------------------------------------
    if(.not. dy_check) then
        ! dy hasn't be set
        write(*,'(A)') "!! A value hasn't been specified for dy! !!"
        stop
    else
        if( dy < 0 ) then
            !The specified value is invalid
            write(*,'(A,F8.3,A)') "!! The value specified for dy is invalid : ", dy, " - need (dy>0) !!"
            stop
        else
            !The specified value is valid
            write(*,'(A)') '- dy: '
            write(*,'(F8.3)') dy
        endif
    endif



    end subroutine get_grid_parameters
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine get_simulation_mode(doc,XML_sim_mode)

    implicit none
    type(node), pointer, intent(in) :: doc
    character(1), intent(out) :: XML_sim_mode
    character(1) :: mode_input
    type(node), pointer :: p_0
    type(nodelist), pointer :: parameterlist, c_0


    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "Simulation_mode")

    if ( getlength(parameterlist) == 0 ) then
        !write(*,'(A)') '<Simulation_mode> not specified!'
        !Set mode to an empty value - later this will either be overwritten if a command-line flag
        !has been used, or set to a default value
        XML_sim_mode = ''
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,'(A)') 'Too many <Simulation_mode> blocks specified!'
        stop

    else

        !look at the sub-blocks of the first <Simulation_mode> block, and extract the relevant values
        p_0 => item(parameterlist, 0)
        c_0 => getchildnodes(p_0)

        !Extract the type tag associated with this <Simulation_mode> block
        if (hasattribute(p_0, "type")) then
            !store the Simulation Mode type tag value
            mode_input = trim( getattribute(p_0, "type") )


            if (mode_input == 'C' .or. mode_input == 'c') then
                !Confined mode specified
                XML_sim_mode = 'C'
            else if (mode_input == 'U' .or. mode_input == 'u') then
                !Unconfined mode specified
                XML_sim_mode = 'U'
            else
                !Invalid type value specified
                write(*,'(A)') ''
                write(*,'(A,A,A)') '!!! Invalid type specified for <Simulation_mode>: ', trim(mode_input), ' !!!'
                write(*,'(A)') ''
                stop
            end if

        else
            !A <Simulation_mode> block appears, but no type tag is specified
            write(*,'(A)') ''
            write(*,*) '!! No type specified for <Simulation_mode>! !!'
            write(*,'(A)') ''
            stop
        end if

    endif


    !Print set values to STDOUT
    write(*,'(A)') '- Simulation Mode: '
    select case (XML_sim_mode)
    case ('C')
        write(*,'(A)') 'Confined'
    case ('U')
        write(*,'(A)') 'Unconfined'
    case ('')
        write(*,'(A)') 'Not set'
        case default
        write(*,'(A,A)') 'Simulation mode type error! : ', XML_sim_mode
        stop
    end select
    write(*,'(A)') ''



    end subroutine get_simulation_mode
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine get_boundary_conditions(doc,h_bc_params,p_bc_params)
    ![get_boundary_conditions] looks for the <boundary_conditions> block and extracts the
    !corresponding boundary condition flags h_bc_params and p_bc_params for the current thickness h and
    !ambient pressure p. these are specified within the sub-blocks <current_thickness> and
    !<ambient_pressure>, respectively. the bc flag vectors each consist of the four values
    !corresponding to the north, east, south, and west boundaries of the rectangular domain.

    implicit none

    type(node), pointer, intent(in) :: doc
    type(node), pointer :: p_0, p_1, p_2
    type(nodelist), pointer :: parameterlist, c_0, c_1
    integer, dimension(0:3), intent(out) :: h_bc_params, p_bc_params
    integer :: i, j, iostat
    character(1000) :: node_path
    type(node), pointer :: node_pointers(:)
    logical :: h_bc_check, p_bc_check = .FALSE.


    !find all of the '<boundary_conditions>...</boundary_conditions>' blocks
    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "boundary_conditions")

    !deal with the wrong number of such blocks
    if ( getlength(parameterlist) == 0 ) then
        write(*,'(A)') '<boundary_conditions> not specified!'
        stop
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,'(A)') 'Too many <boundary_conditions> blocks specified!'
        stop
    endif

    !take the first of these blocks (i.e. ignore repeats if there are any)
    p_0 => item(parameterlist, 0)
    !take all of the sub-blocks within this block
    c_0 => getchildnodes(p_0)

    !initialise h_bc_params and p_bc_params with invalid values, to catch the case where the bcs haven't
    !been specified (valid values are 1 or 2, for dirichlet or neumann bcs respectively)
    h_bc_params = [-1,-1,-1,-1]
    p_bc_params = [-1,-1,-1,-1]


    allocate(node_pointers(3))


    do i = 0, getlength(c_0)-1
        !look within the sub-blocks for 'current_thickness' and 'ambient_pressure',
        !and extract the relevant values

        p_1 => item(c_0, i)

        if (getlocalname(p_1) =="current_thickness") then
            c_1 => getchildnodes(p_1)
            do j=0, getlength(c_1)-1

                p_2 => item(c_1, j)

                node_pointers = [p_0, p_1, p_2]
                node_path = build_node_path(node_pointers)

                select case ( getlocalname(p_2) )

                case ('north')
                    call extractdatacontent(p_2, h_bc_params(0), iostat=iostat)
                case ('east')
                    call extractdatacontent(p_2, h_bc_params(1), iostat=iostat)
                case ('south')
                    call extractdatacontent(p_2, h_bc_params(2), iostat=iostat)
                case ('west')
                    call extractdatacontent(p_2, h_bc_params(3), iostat=iostat)
                case default 
                    iostat = 0
                end select

                call error_check(iostat, node_path)

                !Flag that boundary conditions for h have been set
                h_bc_check = .TRUE.

            enddo

        elseif (getlocalname(p_1) =="ambient_pressure") then

            c_1 => getchildnodes(p_1)
            do j=0, getlength(c_1)-1

                p_2 => item(c_1, j)

                node_pointers = [p_0, p_1, p_2]
                node_path = build_node_path(node_pointers)

                select case ( getlocalname(p_2) )

                case ('north')
                    call extractdatacontent(p_2, p_bc_params(0), iostat=iostat)
                case ('east')
                    call extractdatacontent(p_2, p_bc_params(1), iostat=iostat)
                case ('south')
                    call extractdatacontent(p_2, p_bc_params(2), iostat=iostat)
                case ('west')
                    call extractdatacontent(p_2, p_bc_params(3), iostat=iostat)
                case default 
                    iostat = 0
                end select

                call error_check(iostat, node_path)

                !Flag that boundary conditions for P have been set
                p_bc_check = .TRUE.

            enddo

        endif
    enddo
    
    !Check that h and P boundary conditions have been set
    if (.not. h_bc_check) then
        write(*,*) ''
        write(*,*) '!! No Boundary Conditions for Current_Thickness specified! !!'
        write(*,*) ''
        stop
    elseif (.not. p_bc_check) then
        write(*,*) ''
        write(*,*) '!! No Boundary Conditions for Ambient_Pressure specified! !!'
        write(*,*) ''
        stop
    endif

    !Check if a Neumann BC has been chosen for the pressure on all four sides of the domain.
    !This prevents ambient fluid from leaving the domain as CO2 is injected, and as a result
    !the computation falters.
    if ( all(P_BC_params == 2) ) then
        write(*,*) ''
        write(*,*) '!!! Neumann BC chosen for the pressure on all four sides. This makes injection difficult if Confined! !!!'
        write(*,*) ''
    end if


    !Print set values to STDOUT
    write(*,'(A)') '- h_BC_params: '
    write(*,'(4(I1,1X))') h_BC_params
    write(*,'(A)') '- P_BC_params: '
    write(*,'(4(I1,1X))') P_BC_params
    write(*,'(A)') ''

    end subroutine get_boundary_conditions
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine get_dim_flow_parameters(doc,rho_c,rho_a,rho_a_sat,mu_c,mu_a,&
        &                               s_c_r,s_a_i,krn_m,krw_r,c_sat,g,D_mol,perm_ratio,omega_conv)
    !get_flow_parameters(doc,m,gamma_val,s_c_r,s_a_i,c_sat,q_dissolve)
    ![get_dim_flow_parameters] looks for the <flow_parameters> block and extracts the values for the fixed
    !non-dimensional flow parameters used in the simulation.

    implicit none

    type(node), pointer, intent(in) :: doc
    type(node), pointer :: p_0, p_1
    type(nodelist), pointer :: parameterlist, c_0
    real(wp), intent(out) :: rho_c, rho_a, rho_a_sat, mu_c, mu_a, s_c_r, s_a_i, krn_m, krw_r 
    real(wp), intent(out) :: c_sat, g, D_mol, perm_ratio, omega_conv
    integer :: i, iostat
    character(1000) :: node_path
    type(node), pointer :: node_pointers(:)
    logical :: check_vals(14) = .FALSE.


    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "flow_parameters")

    if ( getlength(parameterlist) == 0 ) then
        write(*,'(A)') '<flow_parameters> not specified!'
        stop
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,'(A)') 'Too many <flow_parameters> blocks specified!'
        stop
    endif

    !look at the sub-blocks of the first <flow_parameters> block, and extract the relevant values
    p_0 => item(parameterlist, 0)
    c_0 => getchildnodes(p_0)

    allocate(node_pointers(2))

    do i = 0, getlength(c_0)-1

        p_1 => item(c_0, i)

        node_pointers = [p_0, p_1]
        node_path = build_node_path(node_pointers)

        select case ( getlocalname(p_1) )

        case ('rho_c')
            call extractdatacontent(p_1, rho_c, iostat=iostat)
            check_vals(1) = .TRUE.
        case ('rho_a_unsat')
            call extractdatacontent(p_1, rho_a, iostat=iostat)
            check_vals(2) = .TRUE.
        case ('rho_a_sat')
            call extractdatacontent(p_1, rho_a_sat, iostat=iostat)
            check_vals(3) = .TRUE.
        case ('mu_c')
            call extractdatacontent(p_1, mu_c, iostat=iostat)
            check_vals(4) = .TRUE.
        case ('mu_a')
            call extractdatacontent(p_1, mu_a, iostat=iostat)
            check_vals(5) = .TRUE.
        case ('s_c_r')
            call extractdatacontent(p_1, s_c_r, iostat=iostat)
            check_vals(6) = .TRUE.
        case ('s_a_i')
            call extractdatacontent(p_1, s_a_i, iostat=iostat)
            check_vals(7) = .TRUE.
        case ('krn_mobile')
            call extractdatacontent(p_1, krn_m, iostat=iostat)
            check_vals(8) = .TRUE.
        case ('krw_residual')
            call extractdatacontent(p_1, krw_r, iostat=iostat)
            check_vals(9) = .TRUE.
        case ('C_sat')
            call extractdatacontent(p_1, c_sat, iostat=iostat)
            check_vals(10) = .TRUE.
        case ('g')
            call extractdatacontent(p_1, g, iostat=iostat)
            check_vals(11) = .TRUE.
        case ('D_mol')
            call extractdatacontent(p_1, D_mol, iostat=iostat)
            check_vals(12) = .TRUE.
        case ('perm_ratio')
            call extractdatacontent(p_1, perm_ratio, iostat=iostat)
            check_vals(13) = .TRUE.
        case ('omega_conv')
            call extractdatacontent(p_1, omega_conv, iostat=iostat)
            check_vals(14) = .TRUE.
        case default 
            iostat = 0
        end select

        call error_check(iostat, node_path)

    end do



    !Check the values that have been assigned. Give an error if invalid or not specified and required.
    !If optional and not required, set to a default value.
    !Finally, print the set values (if valid) to STDOUT
    ! Format:                        parameter, check_flag, valid_flag, required_flag, param_str, unit_str, valid_str, default_val (optional)
    ! 1) Required variables
    call flow_params_check_and_print(rho_c,        check_vals(1),  (0.<rho_c),                 &
    &                                                   .True.,  'rho_c',        'kg m^-3',  'rho_c>0'      )
    
    call flow_params_check_and_print(rho_a,        check_vals(2),  (0.<rho_a),                 &
    &                                                   .True.,  'rho_a',        'kg m^-3',  'rho_a>0'      )
    
    call flow_params_check_and_print(rho_a_sat,    check_vals(3),  (0.<rho_a_sat),             &
    &                                                   .True.,  'rho_a_sat',    'kg m^-3',  'rho_a_sat>0'  )
    
    call flow_params_check_and_print(mu_c,         check_vals(4),  (0.<mu_c),                  &
    &                                                   .True.,  'mu_c',         'Pa s',     'mu_c>0'       )
    
    call flow_params_check_and_print(mu_a,         check_vals(5),  (0.<mu_a),                  &
    &                                                   .True.,  'mu_a',         'Pa s',     'mu_a>0'       )
    
    call flow_params_check_and_print(s_c_r,        check_vals(6),  (0.<=s_c_r .and. s_c_r<1.), &
    &                                                   .True.,  's_c_r',        '-',        '0<=s_c_r<1'   )
    
    call flow_params_check_and_print(s_a_i,        check_vals(7),  (0.<=s_a_i .and. s_a_i<1.), &
    &                                                   .True.,  's_a_i',        '-',        '0<=s_a_i<1'   )
    
    ! 2) Optional variables
    call flow_params_check_and_print(krn_mobile,   check_vals(8),  (0.<krn_mobile),            &
    &                                                   .False., 'krn_mobile',   '-',        'krn_mobile>0',   default_val=1.0_wp)

    call flow_params_check_and_print(krw_residual, check_vals(9),  (0.<krw_residual),          &
    &                                                   .False., 'krw_residual', '-',        'krw_residual>0', default_val=1.0_wp)

    call flow_params_check_and_print(C_sat,        check_vals(10), (0.<=C_sat .and. C_sat<1.), &
    &                                                   .False., 'C_sat',        '-',        '0=<C_sat<1',     default_val=0.04_wp)

    call flow_params_check_and_print(g,            check_vals(11), (0.<g),                     &
    &                                                   .False., 'g',            'm s^-2',   'g>0',            default_val=9.81_wp)

    call flow_params_check_and_print(D_mol,        check_vals(12), (0.<D_mol),                 &
    &                                                   .False., 'D_mol',        'm^2 s^-1', 'D_mol>0',        default_val=2e-9_wp)

    call flow_params_check_and_print(perm_ratio,   check_vals(13), (0.<perm_ratio),            &
    &                                                   .False., 'perm_ratio',   '-',        'perm_ratio>0',   default_val=0.1_wp)

    call flow_params_check_and_print(omega_conv,   check_vals(14), (0.<=omega_conv),           &
    &                                                   .False., 'omega_conv',   '-',        'omega_conv>=0',  default_val=0.0_wp)


    write(*,'(A)') ''


    end subroutine get_dim_flow_parameters
    ! ---------------------------------------------------------------------------------------------


    ! ! ! ! ---------------------------------------------------------------------------------------------
    ! ! ! subroutine get_injection_interval_parameters(doc,z_idx_ceil,z_idx_base)
    ! ! ! ![get_injection_interval_parameters] looks for the <injection_interval_parameters> block and
    ! ! ! !extracts the grid parameters used in the simulation

    ! ! ! implicit none

    ! ! ! type(node), pointer, intent(in) :: doc
    ! ! ! type(node), pointer :: p_0, p_1
    ! ! ! type(nodelist), pointer :: parameterlist, c_0
    ! ! ! integer, intent(out) :: z_idx_ceil, z_idx_base
    ! ! ! integer :: i, j, iostat
    ! ! ! character(1000) :: node_path
    ! ! ! type(node), pointer :: node_pointers(:)

    ! ! ! parameterlist => getelementsbytagnamens(doc, &
    ! ! !     "http://www.xml-cml.org/schema", "injection_interval_parameters")

    ! ! ! if ( getlength(parameterlist) == 0 ) then
    ! ! !     write(*,'(A)') '<injection_interval_parameters> not specified!'
    ! ! ! elseif ( getlength(parameterlist) > 1 ) then
    ! ! !     write(*,'(A)') 'Too many <injection_interval_parameters> blocks specified!'
    ! ! ! endif

    ! ! ! !look at the sub-blocks of the first <grid_parameters> block, and extract the relevant values
    ! ! ! p_0 => item(parameterlist, 0)
    ! ! ! c_0 => getchildnodes(p_0)

    ! ! ! allocate(node_pointers(2))  

    ! ! ! do i = 0, getlength(c_0)-1

    ! ! !     p_1 => item(c_0, i)

    ! ! !     node_pointers = [p_0, p_1]
    ! ! !     node_path = build_node_path(node_pointers) 

    ! ! !     select case ( getlocalname(p_1) )

    ! ! !     case ('z_idx_ceil')
    ! ! !         call extractdatacontent(p_1, z_idx_ceil, iostat=iostat)
    ! ! !     case ('z_idx_base')
    ! ! !         call extractdatacontent(p_1, z_idx_base, iostat=iostat)
    ! ! !     case default 
    ! ! !         iostat = 0
    ! ! !     end select
        
    ! ! !     call error_check(iostat, node_path)

    ! ! ! end do


    ! ! ! end subroutine get_injection_interval_parameters
    ! ! ! ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine get_inj_locs(doc,n_inj_locs,inj_loc_tags,q_inj_locs,nx,ny)
    ![get_inj_locs] looks for the <injection_locations> block and extracts the number of
    !injection locations n_inj_locs, the id tags used to refer to each for reference in
    ![get_flux_info], and the grid coordinates of each injection location.

    implicit none

    type(node), pointer, intent(in) :: doc
    type(node), pointer :: p_0, p_1, p_2
    type(nodelist), pointer :: parameterlist, c_0, c_1
    integer, dimension(:,:), allocatable, intent(out) :: q_inj_locs
    character(10), dimension(:), allocatable, intent(out) :: inj_loc_tags
    integer, intent(out) :: n_inj_locs
    integer, intent(in) :: nx, ny
    integer :: i, j, loc_count, iostat, idx_x, idx_y
    character(1000) :: node_path, id_tag_val
    type(node), pointer :: node_pointers(:)
    logical :: idx_x_check, idx_y_check

    idx_x_check = .FALSE.
    idx_y_check = .FALSE.


    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "injection_locations")

    if ( getlength(parameterlist) == 0 ) then
        write(*,'(A)') '<injection_locations> not specified!'
        stop
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,'(A)') 'Too many <injection_locations> blocks specified!'
        stop
    endif

    !look at the sub-blocks of the first <injection_locations> block
    p_0 => item(parameterlist, 0)
    c_0 => getchildnodes(p_0)

    allocate(node_pointers(3))

    !determine the number of '<location>' tags - getlength(c_0) includes some blank ones
    n_inj_locs = getlength( getelementsbytagnamens(p_0, &
        "http://www.xml-cml.org/schema", "location") )


    !set the size of the array of injection location indices, and the default value
    allocate(q_inj_locs(0:n_inj_locs-1, 0:1))
    q_inj_locs(:,:) = 0

    !set the size of the array of injection location tags, with an initial value
    allocate(inj_loc_tags(0:n_inj_locs-1))
    inj_loc_tags(:) = ''

    !need a separate counter to keep track of the number of <location> blocks
    !encountered, since c_0 includes some 'blank' tags
    loc_count = 0

    do i = 0,getlength(c_0)-1

        p_1 => item(c_0, i)

        if (getlocalname(p_1) =="location") then

            idx_x_check = .FALSE.
            idx_y_check = .FALSE.

            if (hasattribute(p_1, "id")) then
                !store the location id tag value, for reference in get_flux_info
                id_tag_val = trim(getattribute(p_1, "id"))

                if (len_trim(id_tag_val) .lt. 1) then
                    !Empty id tag
                    write(*,'(A,I0,A)') '!! The id tag for injection location ', loc_count+1, ' is empty! !!'
                    stop
                endif

                inj_loc_tags(loc_count) = id_tag_val

                !look at the sub-blocks within each <location> block, and extract the
                !relevant values
                c_1 => getchildnodes(p_1)

                do j = 0, getlength(c_1)-1

                    p_2 => item(c_1, j)

                    node_pointers = [p_0, p_1, p_2]
                    node_path = build_node_path(node_pointers)

                    select case ( getlocalname(p_2) )

                    case ('i')
                        call extractdatacontent( p_2, q_inj_locs(loc_count,0), iostat=iostat )
                        idx_x_check = .TRUE.
                    case('j')
                        call extractdatacontent( p_2, q_inj_locs(loc_count,1), iostat=iostat )
                        idx_y_check = .TRUE.
                    case default 
                        iostat = 0
                    end select
                    
                    call error_check(iostat, node_path)

                enddo

            else
                !encountered a <location> block without an id tag
                write(*,'(A)') ''
                write(*,'(A,I0,A)') '<location> block ', loc_count+1,' is missing an id tag!'
                write(*,'(A)') ''
                stop
            endif

            !Check that both i and j index values were set for this injection location
            if (.not. idx_x_check) then
                !Didn't set an i index value for this location
                write(*,'(A)') ''
                write(*,'(A,I0,A)') '<location> block ', loc_count+1,' is missing an i index value!'
                write(*,'(A)') ''
                stop
            elseif (.not. idx_y_check) then
                !Didn't set an j index value for this location
                write(*,'(A)') ''
                write(*,'(A,I0,A)') '<location> block ', loc_count+1,' is missing a j index value!'
                write(*,'(A)') ''
                stop
            end if


            loc_count = loc_count + 1

        endif

    end do

    if (loc_count .le. 0) then
        write(*,*) '!! No injection locations defined! !!'
        stop
    end if

    ! Check for duplicates among the injection location tags
    call check_string_duplicates(inj_loc_tags, n_inj_locs, 'Injection Location Tags')
    

    !! Check injection indices

    do i = 0,n_inj_locs-1

        idx_x = q_inj_locs(i,0)
        idx_y = q_inj_locs(i,1)

        ! 1) Check injection indices are sensible
        if ( (idx_x .lt. 0) .or. (idx_x .gt. nx-1) .or. (idx_y .lt. 0) .or. (idx_y .gt. ny-1) ) then
            write(*,'(A)') ''
            write(*,'(A,A,A)') '!!! Invalid injection location indices for ', trim(inj_loc_tags(i)), ' !!!'
            write(*,'(A,I0,A,I0,A,I0,A,I0)') ' (i,j) = (', idx_x, ',', idx_y, ') - need 0<= i <=', nx-1,' and 0<=j<=', ny-1
            write(*,'(A)') ''
            stop
        end if

        ! 2) Check injection indices are unique
        do j = i+1,n_inj_locs-1

            if( (q_inj_locs(j,0) == idx_x) .and. (q_inj_locs(j,1) == idx_y) ) then
                !Duplicate pairs of indices found
                write(*,'(A)') ''
                write(*,'(5A)') '!!! Duplicate indices found for ', trim(inj_loc_tags(i)), ' and ', trim(inj_loc_tags(j)), ' !!!'
                write(*,'(A,I0,A,I0,A)') ' (i,j) = (', idx_x, ',', idx_y, ') '
                write(*,'(A)') ''
                stop
            end if
        end do

    end do


    !Print set values to STDOUT
    write(*,'(A)') '- n_inj_locs: '
    write(*,*) n_inj_locs
    write(*,'(A)') '- inj_loc_tags: '
    write(*,*) inj_loc_tags
    write(*,'(A)') '- Q_inj_locs: '
    do i = 0,n_inj_locs-1
        write(*,'(2(I5,1X))') Q_inj_locs(i,0), Q_inj_locs(i,1)
    enddo
    write(*,*)


    end subroutine get_inj_locs
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine get_flux_info(doc,n_inj_locs,inj_loc_tags,n_flux_times,t_flux_vals,q_flux_vals)
    ![get_flux_info] looks for the <injection_schedule> block, and extracts the number of times
    !when the injection fluxes are changed, the corresponding times, and the updated values for
    !each of the injection wells.


    implicit none

    type(node), pointer, intent(in) :: doc
    type(node), pointer :: p_0, p_1, p_2, p_3
    type(nodelist), pointer :: parameterlist, c_0, c_1, c_2
    integer, intent(in) :: n_inj_locs
    character(10), dimension(0:n_inj_locs-1), intent(in) :: inj_loc_tags
    integer, intent(out) :: n_flux_times
    real(wp), dimension(:), allocatable, intent(out) :: t_flux_vals
    real(wp), dimension(:,:), allocatable, intent(out) :: q_flux_vals
    logical, dimension(:,:), allocatable :: flux_val_is_set
    character(100) :: str, loc_tag
    integer :: i, j, k, l, schedule_count, iostat, idx_sort
    character(1000) :: node_path
    type(node), pointer :: node_pointers(:)
    logical :: time_check, rate_check
    real(wp) :: rate_val
    integer :: current_tag_idx
    real(wp), dimension(:), allocatable :: t_flux_vals_temp
    real(wp), dimension(:,:), allocatable :: q_flux_vals_temp
    logical, dimension(:,:), allocatable :: flux_val_is_set_temp
    integer, dimension(:), allocatable :: t_flux_vals_ordering


    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "injection_schedule")

    if ( getlength(parameterlist) == 0 ) then
        write(*,*) '<injection_schedule> not specified!'
        stop
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,*) 'Too many <injection_schedule> blocks specified!'
        stop
    endif

    !look at the sub-blocks of the first <injection_schedule> block
    p_0 => item(parameterlist, 0)
    c_0 => getchildnodes(p_0)

    allocate(node_pointers(4))

    !determine the number of 'schedule' tags - getlength(c_0) includes some blank ones
    n_flux_times = getlength( getelementsbytagnamens(p_0, &
        "http://www.xml-cml.org/schema", "schedule") )


    allocate(t_flux_vals(0:n_flux_times-1))
    t_flux_vals(:) = 0.

    allocate(q_flux_vals(0:n_flux_times-1,0:n_inj_locs-1))
    !set the default value to -1, so that i can track the entries that haven't been specified in the loop
    !below and then set them to their preceeding value
    q_flux_vals = -1.

    allocate(flux_val_is_set(0:n_flux_times-1,0:n_inj_locs-1))
    !This array is used to track which entries in q_flux_vals have been updated from the XML input, and
    !which are missing and then need to take assumed values afterwards
    flux_val_is_set = .FALSE.


    !need a separate counter to keep track of the number of <schedule> blocks
    !encountered, since c_0 includes some 'blank' tags
    schedule_count = 0

    do i = 0,getlength(c_0)-1

        p_1 => item(c_0, i)

        if (getlocalname(p_1) =="schedule") then

            time_check = .FALSE.

            !look at each of the sub-blocks within this <schedule> block
            c_1 => getchildnodes(p_1)

            do j = 0, getlength(c_1)-1

                p_2 => item(c_1, j)

                node_pointers = [p_0, p_1, p_2, p_2]
                node_path = build_node_path(node_pointers)

                select case ( getlocalname(p_2) )

                case('time')
                    !extract the time value corresponding to this flux change, and
                    !store it in t_flux_vals
                    call extractdatacontent( p_2, t_flux_vals(schedule_count), iostat=iostat )
                    call error_check(iostat, node_path)
                    time_check = .TRUE.

                case('location')
                    !store the id tag value for the current location block, so that i can search for it
                    !among the values stored by [get_injection_locations]
                    loc_tag = trim(getattribute(p_2, "id"))

                    !loop over each of the id tags stored by [get_injection_locations] until we have a match,
                    !then extract the corresponding flux value and store it in q_flux_vals in the appropriate
                    !entry
                    ! ---  there is definitely a better way to do this bit !
                    current_tag_idx = -1
                    do k = 0, n_inj_locs-1
                        if ( trim(inj_loc_tags(k)) == trim(loc_tag) ) then
                            current_tag_idx = k
                            exit !Step out of this do loop, we've found the index we're after
                        endif
                    enddo


                    if (current_tag_idx == -1) then
                        ! Current loc_tag value wasn't found among the inj_loc_tags values previously set!!
                        write(*,'(A)') ''
                        write(*,'(A)') ' !! Encountered an id tag that does not match one of the injection location tags &
                        &                       set via the <injection_locations> block !!'
                        write(*,'(A,I0,A,A)') ' <schedule> block ', schedule_count+1, ', location id ', trim(loc_tag)
                        write(*,'(A)') ''
                        stop
                    else
                        ! Current value matches one of the established injection locations
                        ! Now process the associated rate value for this location

                        rate_check = .FALSE.

                        c_2 => getchildnodes(p_2)

                        do l=0, getlength(c_2)-1
                            p_3 => item(c_2, l)

                            node_pointers = [p_0, p_1, p_2, p_3]
                            node_path = build_node_path(node_pointers)

                            if ( getlocalname(p_3) == "rate" ) then
                                !Extract the rate value, and check for standard errors
                                call extractdatacontent( p_3, rate_val, iostat=iostat )
                                call error_check(iostat, node_path)

                                !Check that the specified rate value is sensible
                                !(i.e. non-negative here)
                                if (rate_val .lt. 0) then
                                    write(*,'(A)') ''
                                    write(*,'(A,I0,A,A,A)') '!! Invalid rate value specified in <schedule> block ', &
                                        & schedule_count+1, ', location id ', trim(loc_tag), ' !!'
                                    write(*,'(A,F12.5)') 'rate_val = ', rate_val
                                    write(*,'(A)') ''
                                    stop
                                endif

                                !Record that we have read in a valid rate value
                                rate_check = .TRUE.
                                !Add this to the q_flux_vals array
                                q_flux_vals(schedule_count,current_tag_idx) = rate_val
                                !Record that we have entered a value here
                                flux_val_is_set(schedule_count,current_tag_idx) = .TRUE.
                            endif
                        enddo

                        if (.not. rate_check) then
                            !Rate value wasn't specified
                            write(*,'(A)') ''
                            write(*,'(A,I0,A,A,A)') '!! No rate value specified in <schedule> block ', &
                                & schedule_count+1, ', location id ', trim(loc_tag), ' !!'
                            write(*,'(A)') ''
                            stop
                        endif

                    end if

                end select

            enddo

            schedule_count = schedule_count + 1

            if (.not. time_check) then
                write(*,'(A)') ''
                write(*,'(A,I0,A)') '!! Time not specified for <schedule> block ', schedule_count, ' !!'
                write(*,'(A)') ''
                stop
            end if

        endif

    end do

    !finally, we want to make sure that t_flux_vals, and the corresponding rows of q_flux_vals, are
    !in ascending order.
    allocate(t_flux_vals_temp(0:n_flux_times-1))
    t_flux_vals_temp = t_flux_vals !Set to t_flux_vals, which we want to sort

    allocate(q_flux_vals_temp(0:n_flux_times-1,0:n_inj_locs-1))
    q_flux_vals_temp = q_flux_vals !Set to q_flux_vals, just to initialise with some reasonable values that we will then overwrite

    allocate(flux_val_is_set_temp(0:n_flux_times-1,0:n_inj_locs-1))
    flux_val_is_set_temp = flux_val_is_set !Similarly, initialise with values that we will then overwrite

    allocate(t_flux_vals_ordering(0:n_flux_times-1))
    t_flux_vals_ordering(:) = 0

    !Sort t_flux_vals, and keep the ordering for use in sorting q_flux_vals
    call sort_real_array(t_flux_vals_temp, n_flux_times, t_flux_vals_ordering)

    !Update t_flux_vals to this sorted ordering
    t_flux_vals = t_flux_vals_temp
    !Reorder to rows of q_flux_vals and flux_vals_check to match this sorted ordering
    do i = 0,n_flux_times-1
        idx_sort = t_flux_vals_ordering(i)
        q_flux_vals_temp(i,:)    = q_flux_vals( idx_sort , : )
        flux_val_is_set_temp(i,:) = flux_val_is_set( idx_sort , : )
    end do
    q_flux_vals = q_flux_vals_temp
    flux_val_is_set = flux_val_is_set_temp


    !deal with flux values that weren't changed at each given flux time
    !these values are currently set to -1, so that I can track them
    !For the first row, replace -1 with 0
    do j = 0, n_inj_locs-1
        if ( .not. flux_val_is_set(0,j) ) then
            q_flux_vals(0,j) = 0.
        end if
    enddo
    !For the remaining rows, set to the value in the previous row, i.e. assuming no change
    !in flux at that time for that location
    do i = 1, n_flux_times-1
        do j = 0, n_inj_locs-1
            if ( .not. flux_val_is_set(i,j) ) then
                q_flux_vals(i,j) = q_flux_vals(i-1,j)
            end if
        enddo
    end do


    !Print set values to STDOUT
    write(*,'(A)') '- n_flux_times: '
    write(*,*) n_flux_times
    write(*,'(A)') '- [t_flux_vals(i), Q_flux_vals(i,:)] : '

    write(row_fmt,'(A,I0,A)') '(A,F16.5,', n_inj_locs,'(F8.5),A)'
    do i = 0,n_flux_times-1
        write(*,row_fmt) '[', t_flux_vals(i), Q_flux_vals(i,:) , ']'
    enddo
    write(*,*)


    end subroutine get_flux_info
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine get_plot_times(doc,np,target_plot_times)
    ![get_plot_times] looks for the <plot_times> block, and extracts the nondimensional times at
    !which the simulator should record the results of the simulation.

    implicit none

    type(node), pointer, intent(in) :: doc
    type(node), pointer :: p_0, p_1
    type(nodelist), pointer :: parameterlist, c_0
    real(wp), dimension(:), allocatable, intent(out) :: target_plot_times
    integer, intent(out) :: np
    integer :: i, time_count, iostat
    character(1000) :: node_path
    type(node), pointer :: node_pointers(:)


    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "plot_times")

    if ( getlength(parameterlist) == 0 ) then
        write(*,*) '<plot_times> not specified!'
        stop
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,*) 'Too many <plot_times> blocks specified!'
        stop
    endif

    p_0 => item(parameterlist, 0)
    c_0 => getchildnodes(p_0)

    allocate(node_pointers(2))

    !determine the number of 'time' tags - getlength(c_0) includes some blank ones
    np = getlength( getelementsbytagnamens(p_0, &
        "http://www.xml-cml.org/schema", "time") )


    allocate(target_plot_times(0:np-1))
    target_plot_times(:) = 0

    time_count = 0 !need a separate counter for this, since c_0 includes some 'blank' tags
    do i = 0,getlength(c_0)-1

        p_1 => item(c_0, i)

        node_pointers = [p_0, p_1]
        node_path = build_node_path(node_pointers)

        if (getlocalname(p_1) =="time") then

            call extractdatacontent( p_1, target_plot_times(time_count), iostat=iostat )
            call error_check(iostat, node_path)

            time_count = time_count + 1
        endif

    end do

    !Check if any plot times have been set
    if (np .le. 0) then
        write(*,'(A)') ''
        write(*,'(A)') '!! No plot times specified! !!'
        write(*,'(A)') ''
        stop
    end if

    !Sort target_plot_times so that they are in increasing order
    call sort_real_array(target_plot_times, np)


    !Check for negative values
    if(target_plot_times(0) .lt. 0.) then
        write(*,'(A)') ''
        write(*,'(A)') '!! Specified plot times contain a negative value! !!'
        write(*,'(F12.5)') target_plot_times(0)
        write(*,'(A)') ''
        stop
    end if


    !Print set values to STDOUT
    write(*,'(a)') '- np: '
    write(*,*) np
    write(*,'(a)') '- target_plot_times: '
    write(row_fmt,'(a,i0,a)') '(a,i' , ceiling( log10( real(np-1,8) ) ) , ',a,f16.5)'
    do i = 0,np-1
        write(*,row_fmt) '[', i, ']: ', target_plot_times(i)
    enddo
    write(*,*)


    end subroutine get_plot_times
    ! ---------------------------------------------------------------------------------------------

    
    ! ---------------------------------------------------------------------------------------------
    subroutine get_simulation_parameters(doc, itermax_I, tmax_I, dt_init_I, dt_min_I, errmax_t_step_I, &
        & errmax_P_iter_I, P_iter_omega_I, P_iter_check_I)
    ![get_simulation_parameters] looks for the <simulator_parameters> block, and extracts any user-specified values
    !for parameters that are used for the running of the main solver.

    implicit none

    type(node), pointer, intent(in) :: doc
    type(node), pointer :: p_0, p_1
    type(nodelist), pointer :: parameterlist, c_0
    integer, intent(out) :: itermax_I, P_iter_check_I
    real(wp), intent(out) :: tmax_I, dt_init_I, dt_min_I, errmax_t_step_I, errmax_P_iter_I, P_iter_omega_I
    integer :: i, iostat
    character(1000) :: node_path
    type(node), pointer :: node_pointers(:)
    logical :: check_vals(8) = .FALSE.

    !Placeholder values
    !These should all be positive, so use -1 to indicate no assigned value
    itermax_I       = -1
    P_iter_check_I  = -1

    tmax_I          = -1._wp
    dt_init_I       = -1._wp
    dt_min_I        = -1._wp
    errmax_t_step_I = -1._wp
    errmax_P_iter_I = -1._wp
    P_iter_omega_I  = -1._wp



    parameterlist => getelementsbytagnamens(doc, &
        "http://www.xml-cml.org/schema", "simulator_parameters")

    if ( getlength(parameterlist) == 0 ) then
        !write(*,*) '<simulator_parameters> not specified!'
        ! No <simulator_parameters> block in the XML document - use default values
        check_vals = .FALSE.
    elseif ( getlength(parameterlist) > 1 ) then
        write(*,*) 'Too many <simulator_parameters> blocks specified!'
        stop
    else
        !One and only one <simulator_parameters> block in the XML document - parse it
        p_0 => item(parameterlist, 0)
        c_0 => getchildnodes(p_0)

        allocate(node_pointers(2))

        do i = 0, getlength(c_0)-1

            p_1 => item(c_0, i)

            node_pointers = [p_0, p_1]
            node_path = build_node_path(node_pointers)

            select case ( getlocalname(p_1) )

            case ('max_iterations')
                call extractdatacontent(p_1, itermax_I, iostat=iostat)
                check_vals(1) = .TRUE.
            case ('max_time')
                call extractdatacontent(p_1, tmax_I, iostat=iostat)
                check_vals(2) = .TRUE.
            case ('dt_initial')
                call extractdatacontent(p_1, dt_init_I, iostat=iostat)
                check_vals(3) = .TRUE.
            case ('dt_minimum')
                call extractdatacontent(p_1, dt_min_I, iostat=iostat)
                check_vals(4) = .TRUE.
            case ('t_error_threshold')
                call extractdatacontent(p_1, errmax_t_step_I, iostat=iostat)
                check_vals(5) = .TRUE.
            case ('P_error_threshold')
                call extractdatacontent(p_1, errmax_P_iter_I, iostat=iostat)
                check_vals(6) = .TRUE.
            case ('P_iter_omega')
                call extractdatacontent(p_1, P_iter_omega_I, iostat=iostat)
                check_vals(7) = .TRUE.
            case ('P_iter_check')
                call extractdatacontent(p_1, P_iter_check_I, iostat=iostat)
                check_vals(8) = .TRUE.

                case default
                iostat = 0
            end select

            call error_check(iostat, node_path)

        end do

    endif


    !Assign specified values if they are allowable, or use the preset defaults
    !Finally, print set values to STDOUT

    !Format:: parameter, input_val, set_flag, valid_flag, param_str

    !Integer parameter
    call sim_params_check_and_print_int(itermax,        itermax_I,       check_vals(1), &
    &                                                           (itermax_I > 0),        'itermax',       'itermax>0')
    !Real parameters
    call sim_params_check_and_print_real(tmax,          tmax_I,          check_vals(2), &
    &                                                           (tmax_I > 0),           'tmax',          'tmax>0')

    call sim_params_check_and_print_real(dt_init,       dt_init_I,       check_vals(3), &
    &                                                           (dt_init_I > 0),        'dt_init',       'dt_init>0')

    call sim_params_check_and_print_real(dt_min,        dt_min_I,        check_vals(4), &
    &                                                           (dt_min_I > 0),         'dt_min',        'dt_min>0')

    call sim_params_check_and_print_real(errmax_t_step, errmax_t_step_I, check_vals(5), &
    &                                                           (errmax_t_step_I > 0),  'errmax_t_step', 'errmax_t_step>0')

    call sim_params_check_and_print_real(errmax_P_iter, errmax_P_iter_I, check_vals(6), &
    &                                                           (errmax_P_iter_I > 0),  'errmax_P_iter', 'errmax_P_iter>0')

    call sim_params_check_and_print_real(P_iter_omega,  P_iter_omega_I,  check_vals(7), &
    &                                                           (P_iter_omega_I > 0),   'P_iter_omega',  'P_iter_omega>0')
    
    !Integer parameter
    call sim_params_check_and_print_int(P_iter_check,   P_iter_check_I,  check_vals(8), &
    &                                                           (P_iter_check_I > 0),   'P_iter_check',  'P_iter_check>0')

    write(*,*)

    end subroutine get_simulation_parameters
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    function build_node_path(node_pointers) result(node_path)
        ! This function builds out a string describing the name path in the XML document
        ! for the provided nodes node_path = [p_0, p_1, ..., p_n]

    implicit none
    type(node), pointer :: node_pointers(:)
    type(node), pointer :: p
    character(1000), allocatable :: node_path
    character(1000), allocatable :: node_name
    integer :: count
    
    node_path = ""

    do count = 1, size(node_pointers)

        p => node_pointers(count)
        node_name = getlocalname(p)

        node_path = trim(node_path) // '/' // trim(node_name)

    end do

    end function build_node_path
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine error_check(iostat_val, node_string)

    implicit none
    integer, intent(in) :: iostat_val
    character(*), intent(in) :: node_string
    character(100) :: error_type
    logical :: error_status
    integer :: total_length = 70

    select case ( iostat_val )

        case(0)
            error_status = .False.

        case (-1)
            error_type = 'Too Few Elements'
            error_status = .True.
        case (1)
            error_type = 'Too Many Elements'
            error_status = .True.
        case (2)
            error_type = 'Format/Data-Type Issue'
            error_status = .True.
        case default
            error_type = 'Uncommon Error Type'
            error_status = .True.
    end select


    if (error_status) then
        write(*,'(A)') ''
        write(*,'(A)') repeat('!',total_length)
        write(*,'(A,A)') '   An error has been encountered at ', trim(node_string)
        write(*,'(A,A)') '   Error type: ', trim(error_type)
        write(*,'(A)') repeat('!',total_length)
        write(*,'(A)') ''
        stop
    end if


    end subroutine error_check
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    subroutine remove_quotes(string)

    implicit none
    character(*), intent(inout) :: string
    character(len=len(string)) :: trimmed_string
    integer :: str_length, idx_start, idx_end

    str_length = len_trim(string)

    !Check first character for a quotation mark, and move the start index accordingly
    if ((string(1:1) == '"') .or. (string(1:1) == "'")) then
        idx_start = 2
    else
        idx_start = 1
    end if

    !Check last character for a quotation mark, and move the end index accordingly
    if ((string(str_length:str_length) == '"') .or. (string(str_length:str_length) == "'")) then
        idx_end = str_length-1
    else
        idx_end = str_length
    end if

    !Build the 'trimmed' string, and then set the original string to this for output
    trimmed_string = string(idx_start:idx_end)
    string = trimmed_string


    end subroutine remove_quotes
    ! ---------------------------------------------------------------------------------------------
    
    
    ! ---------------------------------------------------------------------------------------------
    subroutine section_title(title_str)
    ! This subroutine prints a string of a fixed length, with a given title string in the middle and
    ! a space on each side along with hyphens filling the rest of the space.
    ! This is for me to use to display section titles for each of the input types read in above.
    ! e.g.
    !           ------------------- Injection Interval parameters --------------------
    
    implicit none
    character(*), intent(in) :: title_str
    integer :: total_length
    integer :: title_length, left_length, right_length

    
    total_length = 70
    title_length = len_trim(title_str)
    
    if (title_length .eq. 0) then
        
        !Blank title, used to print a line of hyphens with no gaps
        write(*,'(A)') repeat('-', total_length)
        
    else
    
        ! Put half (after rounding) of the hyphens on the left, once we've accounted for the title string
        ! and the padding spaces
        left_length = ( total_length - (title_length+2) )/2 !Integer division should round as I want
        !Put the remainder of the hyphens on the right
        right_length = total_length - title_length - left_length - 2

        !Print this full string to STDOUT
        write(*,'(A,1X,A,1X,A)') repeat('-', left_length) , trim(title_str) , repeat('-', right_length) 
        
    end if

    
    end subroutine section_title
    ! ---------------------------------------------------------------------------------------------

    
    ! ---------------------------------------------------------------------------------------------
    subroutine check_string_duplicates(str_array, N, description_string)

        implicit none
        integer, intent(in) :: N
        character(*), intent(in) :: description_string
        character(*), dimension(0:N-1), intent(in) :: str_array
        character(1000), dimension(0:N-1) :: sorted_array
        integer, dimension(0:N-1) :: ordering
        logical, dimension(N) :: duplicate_flags
        integer :: i
    
        duplicate_flags = .FALSE.
    
        sorted_array = str_array
    
        call sort_string_array(sorted_array,N,ordering)
    
        do i = 0,N-2
            if (sorted_array(i) == sorted_array(i+1)) then
                duplicate_flags(i:i+1) = .TRUE. !Flag both as duplicates
            end if
        end do
    
        if (any(duplicate_flags)) then
            !Duplicates found
            write(*,'(A,A,A)') 'Duplicates found for ', description_string, ':'
            do i = 0,N-1
                if (duplicate_flags(i)) then
                    write(*,'(A)') trim(sorted_array(i))
                end if
            end do
            stop
        end if
    
    
        end subroutine
        ! ---------------------------------------------------------------------------------------------
    
    
        ! ---------------------------------------------------------------------------------------------
        subroutine sort_string_array(array,N,ordering)
        ! This function sorts a 1xN array of strings.
        ! The new ordering of elements can be returned as an optional output, 'ordering'
    
        implicit none
        integer, intent(in) :: N
        character(*), dimension(0:N-1), intent(inout) :: array
        integer, dimension(0:N-1), optional, intent(out) :: ordering
        integer, dimension(0:N-1) :: ordering_temp
        integer :: i, j, idx_temp
        character(1000) :: temp
        logical :: swap_check
    
        !Initial ordering
        ordering_temp = [(i, i=0,N-1)]
    
        !Bubble sort
        do i = 0, N-2
    
            swap_check = .FALSE.
    
            do j = 0, N-2
                if (array(j) > array(j+1)) then
                    !Swap the elements of the array
                    temp        = array(j)
                    array(j)    = array(j+1)
                    array(j+1)  = temp
    
                    !Swap the corresponding indices in the ordering array
                    idx_temp            = ordering_temp(j)
                    ordering_temp(j)    = ordering_temp(j+1)
                    ordering_temp(j+1)  = idx_temp
    
                    swap_check = .TRUE.
                end if
            end do
    
            if (.not. swap_check) then
                exit !No swaps were made, so the array is already sorted
            end if
    
        end do
    
        !Save the new ordering of elements, if it has been requested as an output
        if (present(ordering)) then
            ordering = ordering_temp
        end if
    
    
        end subroutine sort_string_array
        ! ---------------------------------------------------------------------------------------------
    
    
        ! ---------------------------------------------------------------------------------------------
        subroutine sort_real_array(array,N,ordering)
        ! This function sorts a 1xN array of real values.
        ! The new ordering of elements can be returned as an optional output, 'ordering'
    
        implicit none
        integer, intent(in) :: N
        real(wp), dimension(0:N-1), intent(inout) :: array
        integer, dimension(0:N-1), optional, intent(out) :: ordering
        integer, dimension(0:N-1) :: ordering_temp
        integer :: i, j, idx_temp
        real(wp) :: temp
        logical :: swap_check
    
        !Initial ordering
        ordering_temp = [(i, i=0,N-1)]
    
        !Bubble sort
        do i = 0, N-2
    
            swap_check = .FALSE.
    
            do j = 0, N-2
                if (array(j) > array(j+1)) then
                    !Swap the elements of the array
                    temp        = array(j)
                    array(j)    = array(j+1)
                    array(j+1)  = temp
    
                    !Swap the corresponding indices in the ordering array
                    idx_temp            = ordering_temp(j)
                    ordering_temp(j)    = ordering_temp(j+1)
                    ordering_temp(j+1)  = idx_temp
    
                    swap_check = .TRUE.
                end if
            end do
    
            if (.not. swap_check) then
                exit !No swaps were made, so the array is already sorted
            end if
    
        end do
    
        !Save the new ordering of elements, if it has been requested as an output
        if (present(ordering)) then
            ordering = ordering_temp
        end if
    
    
        end subroutine sort_real_array
        ! ---------------------------------------------------------------------------------------------
    
    
        ! ---------------------------------------------------------------------------------------------
        subroutine flow_params_check_and_print(p,set_flag,valid_flag,required_flag,param_str,units_str,valid_str,default_val)
    
        implicit none
        real(wp), intent(inout) :: p
        logical, intent(in) :: set_flag, valid_flag, required_flag
        character(*), intent(in) :: param_str, units_str, valid_str
        real(wp), optional, intent(in) :: default_val
        character(100) :: label_fmt, err_val_fmt, err_req_fmt
    
    
        label_fmt      = "('- ', A, ' [', A, ']:' )"
        err_val_fmt = "(' !!! Invalid value specified for ', A, ' : ', E12.4E3, ' - need ( ', A, ' ) !!!')"
        err_req_fmt = "(' !!! A value is required for ', A, ' : ', E12.4E3, ' !!!')"
    
        if (set_flag) then
            !A value has been specified in the XML document
    
            if (valid_flag) then
                !Valid value
                write(*,label_fmt) param_str, units_str
                if (required_flag) then
                    write(*,'(G20.8)') p
                else
                    write(*,'(G20.8, 3X, A)') p , '(User Specified)'
                end if
            else
                !Invalid value
                write(*,err_val_fmt) param_str, p, valid_str
                stop
            end if
    
        else
            !No value specified in XML document
    
            if (required_flag) then
                !We require a value for this variable!
                write(*,err_req_fmt) param_str, p
                stop
            else
                !An input value isn't required for this variable - set to the default value
                p = default_val
                write(*,label_fmt) param_str, units_str
                write(*,'(G20.8, 3X, A)') p , '(Default)'
            end if
    
        end if
    
        end subroutine flow_params_check_and_print
        ! ---------------------------------------------------------------------------------------------
        

    ! ---------------------------------------------------------------------------------------------
    subroutine sim_params_check_and_print_real(p, input_val, set_flag, valid_flag, param_str, valid_str)

        implicit none
        real(wp), intent(inout) :: p
        real(wp), intent(inout) :: input_val
        logical, intent(in) :: set_flag, valid_flag
        character(*), intent(in) :: param_str, valid_str
        character(100) :: row_fmt, err_fmt
    
        !Logging statement formats for displaying the set/default simulation parameters below
        row_fmt = "('- ', A, ':')"
        err_fmt = "(' !!! Invalid value specified for ', A, ' : ', G20.8, ' - need ( ', A, ' ) !!!')"
    
        if (set_flag) then
            !A value has been specified in the XML document
            if (valid_flag) then
                !Valid value, update global parameter
                p = input_val
                write(*,row_fmt) param_str
                write(*,'(G20.8, 3X, A)') p , '(User Specified)'
            else
                !Invalid value
                write(*,err_fmt) param_str, input_val, valid_str
                stop
            end if
        else
            !No value specified in XML document - stick to default
            write(*,row_fmt) param_str
            write(*,'(G20.8, 3X, A)') p , '(Default)'
        end if
    
    
        end subroutine sim_params_check_and_print_real
        ! ---------------------------------------------------------------------------------------------
    
    
        ! ---------------------------------------------------------------------------------------------
        subroutine sim_params_check_and_print_int(p, input_val, set_flag, valid_flag, param_str, valid_str)
    
        implicit none
        integer, intent(inout) :: p
        integer, intent(inout) :: input_val
        logical, intent(in) :: set_flag, valid_flag
        character(*), intent(in) :: param_str, valid_str
        character(100) :: row_fmt, err_fmt
    
        !Logging statement formats for displaying the set/default simulation parameters below
        row_fmt = "('- ', A, ':')"
        err_fmt  = "(' !!! Invalid value specified for ', A, ' : ', I20, ' - need ( ', A, ' ) !!!')"
    
    
        if (set_flag) then
            !A value has been specified in the XML document
            if (valid_flag) then
                !Valid value, update global parameter
                p = input_val
                write(*,row_fmt) param_str
                write(*,'(I20, 3X, A)') p , '(User Specified)'
            else
                !Invalid value
                write(*,err_fmt) param_str, input_val, valid_str
                stop
            end if
        else
            !No value specified in XML document - stick to default
            write(*,row_fmt) param_str
            write(*,'(I20, 3X, A)') p , '(Default)'
        end if
    
    
        end subroutine sim_params_check_and_print_int
        ! ---------------------------------------------------------------------------------------------
    
    
        ! ---------------------------------------------------------------------------------------------
    subroutine Dom_error_check(E)

        implicit none
        type(DOMException), intent(inout) :: E
        integer :: error_code
    
        if (inException(E)) then
            write(*,'(A)') '!! An error has occured while attempting to read in the XML file !!'
            
            error_code = getExceptionCode(E)
            
            select case(error_code)
                
            !DOM error codes (0<e<200)
            case(1)     
                write(*,'(A)') 'DOM error: INDEX_SIZE_ERR '
            case(2)     
                write(*,'(A)') 'DOM error: DOMSTRING_SIZE_ERR '
            case(3)     
                write(*,'(A)') 'DOM error: HIERARCHY_REQUEST_ERR '
            case(4)     
                write(*,'(A)') 'DOM error: WRONG_DOCUMENT_ERR '
            case(5)     
                write(*,'(A)') 'DOM error: INVALID_CHARACTER_ERR '
            case(6)     
                write(*,'(A)') 'DOM error: NO_DATA_ALLOWED_ERR '
            case(7)     
                write(*,'(A)') 'DOM error: NO_MODIFICATION_ALLOWED_ERR '
            case(8)     
                write(*,'(A)') 'DOM error: NOT_FOUND_ERR '
            case(9)     
                write(*,'(A)') 'DOM error: NOT_SUPPORTED_ERR '
            case(10)    
                write(*,'(A)') 'DOM error: INUSE_ATTRIBUTE_ERR '
            case(11)    
                write(*,'(A)') 'DOM error: INVALID_STATE_ERR '
            case(12)    
                write(*,'(A)') 'DOM error: SYNTAX_ERR '
            case(13)    
                write(*,'(A)') 'DOM error: INVALID_MODIFICATION_ERR '
            case(14)    
                write(*,'(A)') 'DOM error: NAMESPACE_ERR '
            case(15)    
                write(*,'(A)') 'DOM error: INVALID_ACCESS_ERR '
            case(16)    
                write(*,'(A)') 'DOM error: VALIDATION_ERR '
            case(17)    
                write(*,'(A)') 'DOM error: TYPE_MISMATCH_ERR '
            case(51)    
                write(*,'(A)') 'DOM error: INVALID_EXPRESSION_ERR '
            case(52)    
                write(*,'(A)') 'DOM error: TYPE_ERR '
            case(81)    
                write(*,'(A)') 'DOM error: PARSE_ERR '
                write(*,'(A)') 'This is likely due to incomplete or incorrect brackets.'
            case(82)    
                write(*,'(A)') 'DOM error: SERIALIZE_ERR '
                
            !FoX error codes (e>200) 
            case(201) 
                write(*,'(A)') 'FoX error: FoX_INVALID_NODE '
            case(202) 
                write(*,'(A)') 'FoX error: FoX_INVALID_CHARACTER '
            case(203) 
                write(*,'(A)') 'FoX error: FoX_NO_SUCH_ENTITY '
            case(204) 
                write(*,'(A)') 'FoX error: FoX_INVALID_PI_DATA '
            case(205) 
                write(*,'(A)') 'FoX error: FoX_INVALID_CDATA_SECTION '
            case(206) 
                write(*,'(A)') 'FoX error: FoX_HIERARCHY_REQUEST_ERR '
            case(207) 
                write(*,'(A)') 'FoX error: FoX_INVALID_PUBLIC_ID '
            case(208) 
                write(*,'(A)') 'FoX error: FoX_INVALID_SYSTEM_ID '
            case(209) 
                write(*,'(A)') 'FoX error: FoX_INVALID_COMMENT '
            case(210) 
                write(*,'(A)') 'FoX error: FoX_NODE_IS_NULL '
            case(211) 
                write(*,'(A)') 'FoX error: FoX_INVALID_ENTITY '
            case(212) 
                write(*,'(A)') 'FoX error: FoX_INVALID_URI '
            case(213) 
                write(*,'(A)') 'FoX error: FoX_IMPL_IS_NULL '
            case(214) 
                write(*,'(A)') 'FoX error: FoX_MAP_IS_NULL '
            case(215) 
                write(*,'(A)') 'FoX error: FoX_LIST_IS_NULL '
            case(999) 
                write(*,'(A)') 'FoX error: FoX_INTERNAL_ERROR '
                
                case default
                write(*,'(A,I0)') 'Unknown error : ', error_code
    
            end select
            
    
            stop
            
        end if
    
        end  subroutine Dom_error_check
        ! ---------------------------------------------------------------------------------------------
    
    end module CO2GraVISim_XML_input