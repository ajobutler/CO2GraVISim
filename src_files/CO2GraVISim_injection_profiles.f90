module CO2GraVISim_injection_profiles
    !This module loads the specified injection information for use by the main solver.
    !This includes the number and locations of the injection points (n_inj_locs, Q_inj_locs),
    !the specified times at which the injection fluxes change (t_flux_vals), and the corresponding
    !values of the flux for each injection location (Q_flux_vals).

    use CO2GraVISim_global !wp

    implicit none

    integer :: n_inj_locs, n_flux_times
    integer,  dimension(:,:), allocatable :: Q_inj_locs
    real(wp), dimension(:,:), allocatable :: Q_flux_vals
    real(wp), dimension(:),   allocatable :: t_flux_vals


    !Injection files are now passed to the subroutine below as inputs
    !when called in CO2GraVISim_input_parameters
    ! ! !Input files
    ! ! character(50)  ::  inj_locs_file = "./Input/injection_locations.txt"
    ! ! character(50)  ::  inj_prof_file = "./Input/injection_profile.txt"


    contains 
    subroutine read_injection_data(inj_locs_file,inj_prof_file)

        implicit none
        integer :: io, j
        real(wp), dimension(:,:), allocatable :: Q_flux_vals_array

        character(len=*), intent(in)    ::  inj_locs_file
        character(len=*), intent(in)    ::  inj_prof_file

        write(*,*) 'Reading in injection locations'
        !Read in parameters
        open(newunit=io, file=inj_locs_file, status='old', action='read')

        read(io,*) !Skip the first two lines, as they're formatting instructions
        read(io,*)
        read(io,*) n_inj_locs

        !allocate the size of the injection locations array, now that we know it
        allocate(Q_inj_locs(0:n_inj_locs-1, 0:1))

        !read the injection locations into this array
        do j = 0,n_inj_locs-1
            read(io,*) Q_inj_locs(j,:)
        end do

        close(io)



        write(*,*) 'Reading in injection flux information'
        !Read in parameters
        open(newunit=io, file=inj_prof_file, status='old', action='read')

        read(io,*) !Skip the first three lines, as they're formatting instructions
        read(io,*)
        read(io,*)
        read(io,*) n_flux_times

        !allocate arrays now that we know how large they should be
        allocate(t_flux_vals(0:n_flux_times-1))
        allocate(Q_flux_vals(0:n_flux_times-1,0:n_inj_locs-1))
        allocate(Q_flux_vals_array(0:n_flux_times-1,0:n_inj_locs))

        !Read in the full flux data
        do j = 0,n_flux_times-1
            read(io,*) Q_flux_vals_array(j,:)
        end do

        close(io)

        !Separate out the first column as the times at which the fluxes are updated,
        !and the remaining columns as the corresponding flux values for the specified
        !injection locations
        t_flux_vals(:)      = Q_flux_vals_array(:,0)
        Q_flux_vals(:,:)    = Q_flux_vals_array(:,1:n_inj_locs)



    end subroutine read_injection_data



end module CO2GraVISim_injection_profiles