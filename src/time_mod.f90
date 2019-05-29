
module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none

  integer(i4b)                           :: n_t                ! Number of x-values
  real(dp),    allocatable, dimension(:) :: x_t                ! Grid of relevant x-values
  real(dp),    allocatable, dimension(:) :: a_t                ! Grid of relevant a-values

  integer(i4b)                           :: n_eta              ! Number of eta grid poins
  real(dp),    allocatable, dimension(:) :: x_eta              ! Grid points for eta
  real(dp),    allocatable, dimension(:) :: eta, eta2          ! Eta and eta'' at each grid point

contains

  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n0, n1, n2
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init, eta_init, eps, h1, hmin, yp1, yp2

    ! Define two epochs, 1) during and 2) after recombination. A third epoch before recombination is added 
    ! for milestone 3. 

    n0          = 300                        ! Number of grid points before recombination            
    n1          = 200                        ! Number of grid points during recombination
    n2          = 300                        ! Number of grid points after recombination
    n_t         = n0 + n1 + n2               ! Total number of grid points
    z_start_rec = 1630.4d0                   ! Redshift of start of recombination
    z_end_rec   = 614.2d0                    ! Redshift of end of recombination
    z_0         = 0.d0                       ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)   ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)     ! x of end of recombination
    x_0         = 0.d0                       ! x today
    
    n_eta       = 1000                       ! Number of eta grid points (for spline)
    a_init      = 1.d-10                     ! Start value of a for eta evaluation
    x_eta1      = log(a_init)                ! Start value of x for eta evaluation
    x_eta2      = 0.d0                       ! End value of x for eta evaluation
    eta_init    = c*a_init/(H_0*sqrt(Omega_r)) ! Start value of eta for eta evalutation



    ! Define variables for ODE-solver and spline 
    eps = 1.d-10
    hmin = 0
    yp1 = 1.d30 
    yp2 = 1.d30  



    ! ------------------------ Filling in x- and a-grids -------------------------

    allocate(x_t(n_t))
    allocate(a_t(n_t))
    

    ! Before recombination
    dx = (x_start_rec - log(1.d-8))/(n0-1)

    do i = 1, n0
       x_t(i) =  log(1.d-8) + (i-1)*dx
    end do 


    ! During recombination
    dx = (x_end_rec- x_start_rec)/n1

    do i = 1, n1
       x_t(n0+i) = x_start_rec + i*dx
    end do 


    ! After recombination 
    dx = (x_0 - x_end_rec)/n2

    do i = 1, n2
       x_t(n0+n1+i) = x_end_rec + i*dx
    end do


    ! Calculating a-grid 
    a_t = exp(x_t)




    ! ----------- Computing and splining the conformal time at each eta time step ---------------

    allocate(x_eta(n_eta))
    allocate(eta(n_eta))
    allocate(eta2(n_eta))

    
    ! Filling in an independent x-grid for the eta values
    dx = (x_eta2 - x_eta1)/n_eta

    do i=1, n_eta 
       x_eta(i) = x_eta1 + (i-1)*dx
     end do 
    

    ! Setting initial value for eta and integrating eta
    eta(1) = eta_init 
    h1 = abs(dx)

    do i=2, n_eta
       eta(i) = eta(i-1)
       call odeint(eta(i:i), x_eta(i-1), x_eta(i), eps, h1, hmin, derivs, bsstep, output)
    end do


    ! Calling spline function
    call spline(x_eta, eta, yp1, yp2, eta2)


  end subroutine initialize_time_mod



  
  ! --------------------------------- Subroutines ------------------------------------
  
  ! Subroutine calculating the derivative of eta with respect to x
  subroutine  derivs(x, eta, detadx)
    implicit none 
    
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: eta
    real(dp), dimension(:), intent(out) :: detadx
    
    detadx = c/get_H_p(x)
  
  end subroutine derivs


  ! Output function for odeint
  subroutine  output(x, eta)
    implicit none 
    
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: eta
    
  end subroutine output



  ! --------------------------------- Functions ------------------------------------

  ! Function computing H at a given x
  function get_H(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H

    get_H = H_0 * sqrt( (Omega_m + Omega_b)*exp(-3*x) +(Omega_r+Omega_nu)*exp(-4*x) + Omega_lambda)

  end function get_H


  ! Function computing H' = a*H at given x
  function get_H_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_H_p

    get_H_p = exp(x) * get_H(x)

  end function get_H_p


  ! Function computing dH'/dx at given x
  function get_dH_p(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p

    get_dH_p =  -((H_0**2)/(get_H_p(x))) * (0.5*(Omega_b + Omega_m)*exp(-x) + Omega_r*exp(-2.d0*x) - Omega_lambda*exp(2.d0*x))

  end function get_dH_p


  ! Function computing eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none

    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta

    get_eta = splint(x_eta, eta, eta2, x_in)

  end function get_eta

  

  ! ---------------------------- Write to file -------------------------------

  subroutine write_to_file_time_mod
    implicit none 

    real(dp)     :: rho_crit, z, rho_m, rho_b, rho_r, rho_lambda
    integer(i4b) :: i

    ! Writing to files 
    open(1, file = '../Milestone1/data/hubble_x.dat', status = 'replace')      ! Hubble-values with respect to x
    open(2, file = '../Milestone1/data/hubble_z.dat', status = 'replace')      ! Hubble-values with respect to z
    open(3, file = '../Milestone1/data/eta.dat', status = 'replace')           ! Conformal time-values
    open(4, file = '../Milestone1/data/omegas.dat', status = 'replace')        ! Density parameters
    open(5, file = '../Milestone1/data/splined_eta.dat', status = 'replace')   ! Splined conformal time-values

    do i = 1, n_eta

       z = 1/exp(x_eta(i)) - 1
       
       write(1,*) x_eta(i), get_H(x_eta(i))
       write(2,*) z, get_H(x_eta(i))
       write(3,*) x_eta(i), eta(i)

       ! Calculating density parameters 
       rho_crit = rho_c * get_H(x_eta(i))**2 / H_0**2
       rho_m = Omega_m * rho_c * exp(-3*x_eta(i))
       rho_b = Omega_b * rho_c * exp(-3*x_eta(i))
       rho_r = Omega_r * rho_c * exp(-4*x_eta(i))
       rho_lambda = Omega_lambda * rho_c 

       write(4,"(5(E17.8))") x_eta(i), rho_m/rho_crit, rho_b/rho_crit, rho_r/rho_crit, rho_lambda/rho_crit
 
    end do

    do i = 1, n_t
       write(5, "(5(E17.8))") x_t(i), get_eta(x_t(i))
    end do

    close(1)
    close(2)
    close(3)
    close(4)
    close(5)

  end subroutine write_to_file_time_mod


end module time_mod


