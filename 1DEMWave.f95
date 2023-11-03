program ElectromagneticWavePropagation
  implicit none
  real(8) :: maxValue
  real(8) :: roundedVal
  real(8) :: curValue
  character(15) :: outputLine
  integer, parameter :: nx = 30  ! Number of grid points
  integer, parameter :: nt = 500 ! Number of time steps
  real(8), parameter :: c0 = 3.0e8 ! Speed of light (m/s)
  real(8), parameter :: dx = 0.01  ! Spatial grid spacing (m)
  real(8), parameter :: dt = 1.0e-12 ! Time step (s)
  real(8), dimension(nx) :: Ez, Hy ! Electric Field , Magnetic Field
  integer :: i, t

  ! Initialize arrays
  Ez = 0.0
  Hy = 0.0

  ! Main time-stepping loop
  do t = 1, nt
    ! Update Magnetic Field
    do i = 1, nx - 1
      Hy(i) = Hy(i) + (Ez(i + 1) - Ez(i)) * (dt / dx / c0)
    end do

    ! Update Electric field
    do i = 2, nx
      Ez(i) = Ez(i) + (Hy(i) - Hy(i - 1)) * (dt * c0 / dx)
    end do

    ! Source excitation
    Ez(1) = exp(-0.5 * ((real(t) - 30.0) / 10.0)**2)

    ! Print the Ez field
    if (mod(t, 10) == 0) then
      write(*, *) 'Time step:', t, ' Ez at point 1:', Ez(1)
    end if
  end do
end program ElectromagneticWavePropagation
