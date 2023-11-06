program EMWavePropagation_2D
  implicit none
  
  ! Simulation parameters
  integer, parameter :: nx = 50         ! Number of spatial grid points in the x-direction
  integer, parameter :: ny = 50         ! Number of spatial grid points in the y-direction
  integer, parameter :: nt = 50       ! Number of time steps
  real(8), parameter :: c0 = 3.0e8       ! Speed of light (m/s)
  real(8), parameter :: dx = 0.01        ! Spatial grid spacing in the x-direction (m)
  real(8), parameter :: dy = 0.01        ! Spatial grid spacing in the y-direction (m)
  real(8), parameter :: dt = dx / (2.0 * c0) ! Time step (s) using the CFL condition
  real(8), dimension(0:nx-1, 0:ny-1) :: E, H
  real(8), parameter :: pulse_frequency = 1.0e12 ! Frequency of the source (Hz)
  real(8), parameter :: pulse_width = 2.0e-12    ! Width of the Gaussian pulse (s)
  integer :: t, i, j
  real(8), dimension(0:nx-1) :: pml_coeff
  real(8), parameter :: pml_decay_factor = 0.5  ! PML decay factor (0 < pml_decay_factor < 1)

  ! Constants for the Gaussian pulse source
  real(8) :: pulse_center_x, pulse_center_y
  pulse_center_x = 0.5 * real(nx) * dx
  pulse_center_y = 0.5 * real(ny) * dy

  ! Initialize the fields and PML coefficients
  E = 0.0
  H = 0.0
  pml_coeff = 0.6

  ! Time-stepping loop
  do t = 1, nt
    ! Update electric field
    do i = 1, nx-1
      do j = 1, ny-1
        E(i, j) = E(i, j) - c0 * dt * ((H(i, j) - H(i-1, j)) / dx + (H(i, j) - H(i, j-1)) / dy)
      end do
    end do

    ! Update magnetic field
    do i = 0, nx-2
      do j = 0, ny-2
        H(i, j) = H(i, j) - c0 * dt * ((E(i+1, j) - E(i, j)) / dx + (E(i, j+1) - E(i, j)) / dy)
      end do
    end do
    
    ! Apply a Gaussian pulse source at the center
    do i = 0, nx-1
      do j = 0, ny-1
        E(i, j) = E(i, j) + exp(-((real(i)*dx - pulse_center_x)**2 / (pulse_width**2) + &
                            (real(j)*dy - pulse_center_y)**2 / (pulse_width**2))) * &
                            sin(2.0 * 3.14159265359 * pulse_frequency * t * dt)
      end do
    end do
    ! Print some data at center every 10 time steps
    if (mod(t, 10) == 0) then
      print *, "Time step:", t
      print *, "Electric field at center:", E(nx/2, ny/2)
    end if
    ! Save the electric field data at each time step
    write(11) E
  end do
end program EMWavePropagation_2D
