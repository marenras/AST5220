
module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                             ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec                         ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22              ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: log_tau, log_tau2, log_tau22  ! Splined (log of) tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2                     ! Splined electron density, n_e
  real(dp), allocatable, dimension(:), private :: log_n_e, log_n_e2             ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: g, g2, g22                    ! Splined visibility function
  real(dp), allocatable, dimension(:), private :: X_e                           ! Fractional electron density, n_e / n_H

contains

  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b) :: i, j, k
    real(dp)     :: saha_limit, y, T_b, n_b, dydx, xmin, xmax, dx, f, n_e0, X_e0, xstart, xstop, eps, h1, hmin,  yp1, yp2, const, z
    logical(lgt) :: use_saha


    ! Define variables for ODE-solver and spline 
    eps = 1.d-8
    hmin = 1.d-5
    yp1 = 1.d30 
    yp2 = 1.d30  


    ! Defining grid-parameters
    saha_limit = 0.99d0       ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)  ! Start grids at a = 10^-10
    xstop      = 0.d0         ! Stop  grids at a = 1
    n          = 1000         ! Number of grid points between xstart and xstop


    ! Allocating arrays
    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(log_tau(n))
    allocate(log_tau2(n))
    allocate(log_tau22(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(log_n_e(n))
    allocate(log_n_e2(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))



    ! Filling in grid with x-values 
    dx = (xstop-xstart)/(n-1.d0)
    
    do i=1, n+1
       x_rec(i) = xstart + (i-1.d0)*dx
    end do


    ! ----------------- Computing X_e and n_e at all grid times -------------------

    use_saha = .true.
    
    h1 = abs(dx)
    
    do i = 1, n
       ! Computing the baryon density 
       n_b = (Omega_b * rho_c)/(m_H * exp(3*x_rec(i)))

       if (use_saha) then
          ! Use the Saha equation
          T_b = T_0 / exp(x_rec(i))
          const = 1.d0/n_b * ((m_e * T_b * k_b)/(2.d0*pi* hbar**2))**(1.5d0) * exp(-epsilon_0/(k_b*T_b))
  
          !X_e(i) = (- const + sqrt(const**2 + 4.d0*const)) / 2.d0     ! Standard quadratic formula
          X_e(i) = 2.d0/(1.d0 + sqrt(1.d0 + 4.d0/const))               ! Stable quadratic formula
          

          if (X_e(i) < saha_limit) use_saha = .false.
       else
          ! Use the Peebles equation
          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), eps, h1, hmin, derivs_X_e, rkqs, output)

       end if

       ! Computing electron density 
       n_e(i) = X_e(i) * n_b

    end do


    ! Computing splined (log of) electron density function
    log_n_e = log(n_e)
    call spline(x_rec, log_n_e, yp1, yp2, n_e2) 




    ! -------------- Computing the optical depth at all grid points -----------------
    
    ! Initial condition
    tau(n) = 0.d0

    ! Integration
    do i = n-1, 1, -1
       tau(i) = tau(i+1)
       call odeint(tau(i:i), x_rec(i+1), x_rec(i), eps*1.d-3, h1, hmin, derivs_tau, rkqs, output)
    end do 


    ! Setting tau(n) equal to the calculated value for tau(n-1) to avoid log(0)
    tau(n) = tau(n-1)

    log_tau = log(tau)


    ! Computing splined (log of) optical depth and 
    ! splined second derivative of (log of) optical depth

    call spline(x_rec, log_tau, yp1, yp2, log_tau2)
    call spline(x_rec, log_tau2, yp1, yp2, log_tau22) 




    ! -------------------- Computing the visibility function -----------------------

    ! Computing splined visibility function and
    ! splined second derivative of visibility function

    do i = 1, n
       g(i) = - get_dtau(x_rec(i))*exp(-get_tau(x_rec(i))) 
    end do

    call spline(x_rec, g, yp1, yp2, g2)
    call spline(x_rec, g2, yp1, yp2, g22) 



  end subroutine initialize_rec_mod




  
  ! ----------------------------- Subroutines ---------------------------------

  ! Subroutine calculating the derivative of X_e
  subroutine derivs_X_e(x, X_e, dX_edx)
    implicit none

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: X_e
    real(dp), dimension(:), intent(out) :: dX_edx

    real(dp) :: phi_2, alpha_2, beta, beta_2, n_1s, lambda_alpha, lambda_2s_1s, C_r, T_b, n_H

    T_b = T_0 / exp(x)
    n_H = (Omega_b * rho_c)/(m_H * exp(3.d0*x))
    phi_2 = 0.448d0 * log(epsilon_0/(T_b*k_b))
    alpha_2 = (64.d0*pi)/(sqrt(27.d0*pi)) * (alpha/m_e)**2.d0 * sqrt(epsilon_0/(T_b*k_b)) * phi_2 * hbar**2 / c
    beta = alpha_2 * ((m_e*T_b*k_b)/(2.d0*pi*hbar**2))**(1.5d0) * exp(-epsilon_0/(T_b*k_b)) 
    !beta_2 = beta * exp(3.d0*epsilon_0/(4.d0*T_b*k_b))                                                 ! Original beta_2
    beta_2 = alpha_2 * ((m_e*T_b*k_b)/(2.d0*pi*hbar**2))**(1.5d0) * exp(-epsilon_0/(4.d0*T_b*k_b))      ! More stable beta_2
    n_1s = (1-X_e(1))*n_H
    lambda_alpha = get_H(x) * (3.d0*epsilon_0/(hbar*c))**3.d0/((8.d0*pi)**2 * n_1s)
    lambda_2s_1s = 8.227d0 
    C_r = (lambda_2s_1s + lambda_alpha)/(lambda_2s_1s + lambda_alpha + beta_2) 


    dX_edx = (C_r/get_H(x)) * (beta*(1.d0-X_e(1)) - n_H*alpha_2*X_e(1)**2)


  end subroutine derivs_X_e
  

  ! Subroutine computing the derivative of tau
  subroutine derivs_tau(x, tau, dtaudx)
    implicit none 

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: tau
    real(dp), dimension(:), intent(out) :: dtaudx

    dtaudx = - (get_n_e(x)*sigma_T*exp(x)*c)/get_H_p(x)

  end subroutine derivs_tau


  
  ! ----------------------------- Functions ---------------------------------

  ! Function for computing n_e at arbitrary x
  function get_n_e(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
   
    get_n_e = exp(splint(x_rec, log_n_e, log_n_e2, x))

  end function get_n_e



  ! Function for computing tau at arbitrary x
  function get_tau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_tau

    get_tau = exp(splint(x_rec, log_tau, log_tau2, x))

  end function get_tau


  ! Function for computing the derivative of tau at arbitrary x
  function get_dtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dtau

    get_dtau = get_tau(x) * splint_deriv(x_rec, log_tau, log_tau2, x)

  end function get_dtau


  ! Function for computing the second derivative of tau at arbitrary x
  function get_ddtau(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau
    
    get_ddtau = get_tau(x)*splint(x_rec, log_tau2, log_tau22, x) + get_dtau(x)**2.d0/get_tau(x)

  end function get_ddtau


  ! Function for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_g

    get_g = splint(x_rec, g, g2, x)

  end function get_g


  ! Function for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_dg

    get_dg = splint_deriv(x_rec, g, g2, x)

  end function get_dg


  ! Function for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none

    real(dp), intent(in) :: x
    real(dp)             :: get_ddg

    get_ddg = splint(x_rec, g2, g22, x)

  end function get_ddg



  ! ----------------------------- Write to file ----------------------------

  subroutine write_to_file_rec_mod
    implicit none 

    integer(i4b) :: i
    real(dp)     :: z

    ! Writing to files 
    open(1, file = '../Milestone2/data/X_e.dat', status = 'replace')    
    open(2, file = '../Milestone2/data/tau.dat', status = 'replace') 
    open(3, file = '../Milestone2/data/g.dat',   status = 'replace') 
  
    do i = 1, n
       z = 1.d0/exp(x_rec(i)) - 1                  ! Redshift

       write(1,"(11(E17.8E3))") z, X_e(i)
       write(2,"(11(E17.8E3))") x_rec(i), get_tau(x_rec(i)), get_dtau(x_rec(i)), get_ddtau(x_rec(i)) 
       write(3,"(11(E17.8E3))") x_rec(i), get_g(x_rec(i)), get_dg(x_rec(i)), get_ddg(x_rec(i))

    end do

    close(1)
    close(2)
    close(3)



  end subroutine write_to_file_rec_mod


end module rec_mod


