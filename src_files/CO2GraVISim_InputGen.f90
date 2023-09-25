program CO2GraVISim_InputGen
   !This program can be used to generate permeability, porosity, and topography arrays 
   !for use with CO2GraVISim.
   !
   !There are four choices for the porosity distribution: Uniform, Low-porosity block,
   !High-porosity channel, and Checkerboard, specified via user input.
   !The permeability array is then generated from this via the Kozeny-Carman relation.
   !
   !The user can the generate sloped ceilingtopography, specifying the nondimensional 
   !slope in x and y separately. The basement topography is then offset from this by
   !a uniform depth of 1.


   use CO2GraVISim_global

   implicit none 
   real(wp), parameter :: PI = 3.14159265358979
   real(wp), dimension(:,:), allocatable :: ceil_topo, base_topo, poro, perm!, Q
   real(wp), dimension(:), allocatable :: x_vect, y_vect
   integer :: nx, ny, i, j, io, ierror1, ierror2, Poro_type
   real(wp) :: dx, dy, perm_scale

   !For the checkerboard pattern
   integer :: tile_nx, tile_ny                  !Cell dimensions
   integer :: num_tiles_x, num_tiles_y          !Number of points in x,y per tile
   real(wp) :: phi_a, phi_b, phi_sum, phi_diff  !Porosity values for Checkboard pattern

   !For the sloped topography
   real(wp) :: slope_angle_x, slope_angle_y   !Slope angles for reservoir in x and y directions

   !Input file
   character(50)  ::  grid_params_file = "./Input/grid_parameters.txt"

   !Output files
   character(50)  ::  ceil_topo_file   = "./Input/ceil_topo.txt"
   character(50)  ::  base_topo_file   = "./Input/base_topo.txt"
   character(50)  ::  poro_file        = "./Input/porosity.txt"
   character(50)  ::  perm_file        = "./Input/permeability.txt"

   write(*,*) 'Reading in grid parameters'
   !Read in parameters
   open(newunit=io, file=grid_params_file, status='old', action='read')
   read(io,*) !Skip the first line, as it's formatting instructions
   read(io,*) nx
   read(io,*) ny
   read(io,*) dx
   read(io,*) dy
   close(io)

   !Allocate arrays now we know the size of them
   allocate(ceil_topo(0:nx-1,0:ny-1))
   allocate(base_topo(0:nx-1,0:ny-1))
   allocate(poro(0:nx-1,0:ny-1))
   allocate(perm(0:nx-1,0:ny-1))

   allocate(x_vect(0:nx-1))
   allocate(y_vect(0:ny-1))

   !Construct the nondimensional spatial grid vectors for use in defining the arrays below
   do i = 0,nx-1
      x_vect(i) = ( i - nint((nx-1)/2._wp) ) * dx
   end do

   do j = 0,ny-1
      y_vect(j) = ( j - nint((ny-1)/2._wp) ) * dy
   end do


   ! -- Porosity -----------------------------------------------------------------------------------------

   !Ask the user to choose one of the preset porosity distribtion:
   ! 1)Uniform 2)Low-perm_rectangle 3)High-perm channel 4)Checkerboard
   write(*,*) 'Choose a Porosity distribution:'
   write(*,*) '1)Uniform 2)Low-porosity rectangle 3)High-porosity channel 4)Checkerboard'

   read(*,'(i10)',iostat=ierror1) Poro_type



   Select Case (Poro_type)
    case (1)
      write(*,*) 'You chose Uniform (1)'

      poro = 1._wp


    case (2)
      write(*,*) 'You chose Low-Porosity Rectangle (2)'

      write(*,*) 'Choose a value for (phi_min / phi_max):'
      read *, phi_b

      !Build a low-porosity rectangle a prescribed nondimensional distance away from the origin
      do j = 0,ny-1
         do i= 0,nx-1

            if (x_vect(i)>=-3 .and. x_vect(i)<=3 .and. y_vect(j)>=2 .and. y_vect(j)<=5 ) then
               poro(i,j) = phi_b
            else
               poro(i,j) = 1._wp
            end if

         end do
      end do



    case (3)
      write(*,*) 'You chose High-Porosity Channel (3)'

      write(*,*) 'Choose a value for (phi_min / phi_max):'
      read *, phi_b

      !Build a high-porosity channel a prescribed nondimensional distance away from the origin
      do i= 0,nx-1

         if (x_vect(i)>=2 .and. x_vect(i)<=5 ) then
            poro(i,:) = 1._wp
         else
            poro(i,:) = phi_b
         end if

      end do


    case (4)
      write(*,*) 'You chose Checkerboard (4)'

      print*, "Choose the no. of points in x and y of a checkerboard tile:"
      write(*,*) 'grid dimensions are (Nx,Ny) = ', nx, ny
      !Read something from stdin, and then output it in stdout
      read(*,'(i10)',iostat=ierror2) tile_nx, tile_ny
      !!Check the type of the answer given is appropriate
      if ( ierror2 /= 0 ) then
         write(*,*) ierror2
         write(*,*) 'An error occured'
         write(*,*) 'You entered ', tile_nx, tile_ny
         write(*,*) 'These should be integers. Please try again.'
         stop  !Error has occurred, so stop the entire script
      endif
      write(*,*) 'Checkerboard values are ', tile_nx, tile_ny

      write(*,*) 'Choose a value for (phi_min / phi_max):'
      read *, phi_b

      !Number of tiles that will span the x and y directions
      num_tiles_x = ceiling( float(nx-1) / float(tile_nx) )
      num_tiles_y = ceiling( float(ny-1) / float(tile_ny) )

      write(*,*) num_tiles_x, num_tiles_y


      !Permeability values, starting with phi_a in the top left
      phi_a      = 1.0_wp
      ! phi_b      = 0.8_wp
      phi_sum    = phi_a + phi_b
      phi_diff   = phi_a - phi_b

      write(*,*) 'Checkerboard parameters:', tile_nx, tile_ny, num_tiles_x, num_tiles_y

      !Start by constructing a checkerboard pattern of +/- 1 values
      ! Build up the first column
      do i = 0,num_tiles_x - 2
         ! All but the last row, first column
         poro(i * tile_nx : (i + 1) * tile_nx, 0:tile_ny) = (-1.0) ** (i)
      end do
      ! Last row, first column
      poro((num_tiles_x - 1) * tile_nx : nx-1, 0:tile_ny) = (-1.0) ** (num_tiles_x - 1)

      ! Build the remaining columns by alternating the sign from the first column. 
      ! First do this for the 2nd to penultimate columns
      do j = 1, num_tiles_y - 2
         poro(:, j * tile_ny : (j + 1) * tile_ny) = (-1.0) ** (j) * poro(0:nx-1, 0:tile_ny)
      end do
      ! Now do this for the final column
      poro(0:nx-1, (num_tiles_y - 1) * tile_ny : ny-1) = &
      &                          ( (-1.0) ** (num_tiles_y - 1) ) * poro(0:nx-1, 0:tile_ny)
      

      ! Switch from +/- 1 values to ka/kb values
      poro = 0.5*( phi_sum + phi_diff * poro )

    case default
      write(*,*) 'You chose incorrectly!'
      stop

   end select


   ! -- Permeability -------------------------------------------------------------------------------------------------

   !! Build the permeability from the porosity using the Kozeny-Carman relation
   !Scaling porosity by a nominal dimensional value of 15% just so that this formula works (i.e. 1-phi isn't zero)
   perm = ( (0.15*poro)**3 ) * ( 1._wp - (0.15*poro) )**(-2)

   !Rescale so that the mean value is 1, for nondimensional form
   perm_scale = sum(perm) / (nx*ny)
   perm = perm / perm_scale


   ! -- Reservoir Topography -----------------------------------------------------------------------------------------

   write(*,*) 'Choose a topographic slope in x:'
   read *, slope_angle_x
   write(*,*) 'Choose a topographic slope in y:'
   read *, slope_angle_y

   !Ceiling topography
   ceil_topo = 0._wp
   !Slope of specified angles in the x and y directions
   do i=0,nx-1
      do j=0,ny-1
         ceil_topo(i,j) = tan(slope_angle_x *pi/180._wp) * dx * (i - (nx-1)/2) &
         &              + tan(slope_angle_y *pi/180._wp) * dy * (j - (ny-1)/2)
      end do
   end do

   !Basement Topography
   base_topo = ceil_topo + 1._wp
   ! base_topo = max(base_topo,ceil_topo) !To avoid base_topo passing through ceil_topo, max because we're measuring downward!



   !-- Output Data --------------------------------------------------------------------------------

   open(newunit=io, file=ceil_topo_file)
   do j=0,ny-1
      write(io,*) ceil_topo(0:nx-1,j)
   end do

   open(newunit=io, file=base_topo_file)
   do j=0,ny-1
      write(io,*) base_topo(0:nx-1,j)
   end do

   open(newunit=io, file=poro_file)
   do j=0,ny-1
      write(io,*) poro(0:nx-1,j)
   end do

   open(newunit=io, file=perm_file)
   do j=0,ny-1
      write(io,*) perm(0:nx-1,j)
   end do

   write(*,*) 'Porosity, Permeability, Topography, and Injection array built'

end program CO2GraVISim_InputGen
