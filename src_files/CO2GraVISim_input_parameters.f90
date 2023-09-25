module CO2GraVISim_input_parameters
    !This module is used to load in various input parameters that are specific to an 
    !individual run, such as reservoir topography and permeabililty arrays, so that 
    !other modules and subroutines can access them.

    use CO2GraVISim_global

    implicit none

    integer :: nx, ny, np
    real(wp), dimension(:,:),   allocatable :: H0, B0, D0
    real(wp), dimension(:,:),   allocatable :: perm_h, poro_h
    integer, dimension(0:3) :: h_BC_params, P_BC_params
    real(wp) :: dx, dy
    real(wp) :: M, Gamma_val, s_c_r, s_a_i, C_sat, q_dissolve

    !Input files
    character(50)  ::  flow_params_file         = "./Input/flow_parameters.txt"
    character(50)  ::  grid_params_file         = "./Input/grid_parameters.txt"
    character(50)  ::  BC_params_file           = "./Input/boundary_conditions.txt"
    character(50)  ::  ceil_topo_file           = "./Input/ceil_topo.txt"
    character(50)  ::  base_topo_file           = "./Input/base_topo.txt"
    character(50)  ::  poro_h_file              = "./Input/porosity.txt"
    character(50)  ::  perm_h_file              = "./Input/permeability.txt"


    contains 
    subroutine read_flow_parameters

        implicit none
        integer :: io

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

    end subroutine read_flow_parameters


    subroutine read_reservoir_inputs

        implicit none
        integer :: io, j_y

        ! -- Read in Inputs ----------------------------------------------------------------------------------------------

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


    end subroutine read_reservoir_inputs

end module CO2GraVISim_input_parameters