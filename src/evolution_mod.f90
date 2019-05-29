
module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  use spline_2D_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private :: a_init   = 1.d-8
  real(dp),     parameter, private :: x_init   = log(a_init)
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100
  integer(i4b), parameter, private :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current
  integer(i4b), private :: npar = 6+lmax_int
  real(dp), allocatable, dimension(:)     :: dydx

  ! Variables for hires source function
  integer(i4b), parameter            :: n_hires = 5000



contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k_hires, x_hires, S)
    implicit none

    real(dp), pointer, dimension(:),   intent(out) :: k_hires, x_hires 
    real(dp), pointer, dimension(:,:), intent(out) :: S
    real(dp), allocatable,     dimension(:,:)      :: S_lores
    real(dp), allocatable, dimension(:,:,:,:)      :: S_coeff

    integer(i4b) :: i, j, k
    real(dp)     :: g, dg, ddg, tau, dt, ddt, H_p, dH_p, ddHH_p, Pi, dPi, ddPi, ck, x_0


    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    !
    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    !   2) Then spline this function with a 2D spline
    !   3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays
    
    
    allocate(x_hires(n_hires))
    allocate(k_hires(n_hires))
   
    x_0 = 0.d0

    do i=1, n_hires
       x_hires(i) = x_init + (x_0 - x_init)*(i-1.d0)/(n_hires-1.d0)
       k_hires(i) = k_min  + (k_max - k_min)*(i-1.d0)/(n_hires-1.d0)
    end do
  
    
    allocate(S_lores(n_t, n_k))        ! Low resolution source function
    allocate(S(n_hires, n_hires))      ! High resolution source function
    allocate(S_coeff(4, 4, n_t, n_k))


    ! --------------- Calculating the source function for all values of k and x ----------------

    do k=1, n_k
       k_current = ks(k)
       ck = c*k_current

       do j=1, n_t
          g    = get_g(x_t(j))
          dg   = get_dg(x_t(j))
          ddg  = get_ddg(x_t(j))
          tau  = get_tau(x_t(j))
          dt   = get_dtau(x_t(j))
          ddt  = get_ddtau(x_t(j))
          H_p  = get_H_p(x_t(j))
          dH_p = get_dH_p(x_t(j))
          Pi   = Theta(j,2,k)
          dPi  = dTheta(j,2,k)
         
          ddPi = (2.d0*ck)/(5.d0*H_p) * (- dH_p/H_p*Theta(j,1,k) + dTheta(j,1,k)) & 
                 + 3.d0/10.d0 * (ddt*Pi + dt*dPi) & 
                 - (3.d0*ck)/(5.d0*H_p) * (-dH_p/H_p*Theta(j,3,k) + dTheta(j,3,k))


          ddHH_p = H_0**2/2.d0 * ((Omega_b + Omega_m)/exp(x_t(j)) + & 
                  4.d0*Omega_r/exp(2.d0*x_t(j)) + 4.d0*Omega_lambda*exp(2.d0*x_t(j)))

          S_lores(j,k) = g*(Theta(j,0,k) + Psi(j,k) + 1.d0/4.d0*Pi) + exp(-tau)*(dPsi(j,k) - dPhi(j,k)) - & 
                        1.d0/ck * (dH_p*g*v_b(j,k) + H_p*dg*v_b(j,k) + H_p*g*dv_b(j,k)) + &
                        3.d0/(4.d0*ck**2.d0) * (g*Pi*ddHH_p + 3.d0*H_p*dH_p*(dg*Pi + g*dPi) + & 
                        H_p**2.d0*(ddg*Pi + 2.d0*dg*dPi + g*ddPi))

       end do
    end do 
    


    ! Splining the source function
    call splie2_full_precomp(x_t, ks, S_lores, S_coeff)

    ! Storing the splined source funtion in array
    do k=1, n_hires
      do j=1, n_hires
          S(j,k) = splin2_full_precomp(x_t, ks, S_coeff, x_hires(j), k_hires(k))
       end do
    end do


  end subroutine get_hires_source_function




  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none

    integer(i4b) :: l, i, k

    ! Initializing k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do k=1, n_k
       ks(k) = k_min + (k_max - k_min) * ((k-1.d0)/(n_k-1.d0))**2.d0
    end do 

    ! Allocate arrays for perturbation quantities
    allocate(Theta(1:n_t, 0:lmax_int, n_k))
    allocate(delta(1:n_t, n_k))
    allocate(delta_b(1:n_t, n_k))
    allocate(v(1:n_t, n_k))
    allocate(v_b(1:n_t, n_k))
    allocate(Phi(1:n_t, n_k))
    allocate(Psi(1:n_t, n_k))
    allocate(dPhi(1:n_t, n_k))
    allocate(dPsi(1:n_t, n_k))
    allocate(dv_b(1:n_t, n_k))
    allocate(dTheta(1:n_t, 0:lmax_int, n_k))



    ! ---------- Setting up initial conditions for the Boltzmann and Einstein equations ----------

    Theta(:,:,:)   = 0.d0
    dTheta(:,:,:)  = 0.d0
    dPhi(:,:)      = 0.d0
    dPsi(:,:)      = 0.d0

    Phi(1,:)     = 1.d0
    Psi(1,:)     = - Phi(1,:)
    delta(1,:)   = 3.d0/2.d0 * Phi(1,:)
    delta_b(1,:) = delta(1,:)
    Theta(1,0,:) = 1.d0/2.d0 * Phi(1,:)
       
    do i = 1, n_k
       v(1,i)       = (c*ks(i))/(2.d0*get_H_p(x_init))*Phi(1,i)
       v_b(1,i)     = v(1,i)
       Theta(1,1,i) = - (c*ks(i))/(6.d0*get_H_p(x_init)) * Phi(1,i)
       Theta(1,2,i) = - (20.d0*c*ks(i))/(45.d0*get_H_p(x_init)*get_dtau(x_init)) * Theta(1,1,i) 
       do l = 3, lmax_int
          Theta(1,l,i) = -l/(2.d0*l+1.d0) * (c*ks(i))/(get_H_p(x_init)*get_dtau(x_init)) * Theta(1,l-1,i)
       end do
    end do

  end subroutine initialize_perturbation_eqns



  subroutine integrate_perturbation_eqns
    implicit none

    integer(i4b) :: i, j, k, l, j_tc
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, H_p, dt, t1, t2
    real(dp)     :: a, ckH, dtau

    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx


    eps    = 1.d-8
    hmin   = 0.d0
    h1     = 1.d-5

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))

    y_tight_coupling = 0.d0
    y                = 0.d0
    dydx             = 0.d0


    
    ! -------------------- Integrating perturbation equations -----------------------
    
    write(*,*) '--------------------------------------------------------------'
    write(*,*) 'Integrating perturbation equations'
    
    ! Propagate each k-mode independently
    do k = 1, n_k
       
       write(*,*) 'k =', k

       k_current = ks(k)  ! Store k_current as a global module variable

       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(1,k)
       y_tight_coupling(2) = delta_b(1,k)
       y_tight_coupling(3) = v(1,k)
       y_tight_coupling(4) = v_b(1,k)
       y_tight_coupling(5) = Phi(1,k)
       y_tight_coupling(6) = Theta(1,0,k)
       y_tight_coupling(7) = Theta(1,1,k)
       
       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       ! Integrating from x_init until the end of tight coupling, using
       ! the tight coupling equations
       
       j = 2

       do while (x_t(j) < x_tc)
          a = exp(x_t(j))
          ckH = c*k_current/get_H_p(x_t(j))
          dtau = get_dtau(x_t(j))
          
          ! Integrate next step 
          call odeint(y_tight_coupling, x_t(j-1), x_t(j), eps, h1, hmin, derivs_ytc, bsstep, output)
          
          ! Save current values 
          delta(j,k)    = y_tight_coupling(1)
          delta_b(j,k)  = y_tight_coupling(2)
          v(j,k)        = y_tight_coupling(3)
          v_b(j,k)      = y_tight_coupling(4) 
          Phi(j,k)      = y_tight_coupling(5)
          Theta(j,0,k)  = y_tight_coupling(6)
          Theta(j,1,k)  = y_tight_coupling(7)
          Theta(j,2,k)  = - 20.d0*ckH/(45.d0*dtau) * Theta(j,1,k)
          Psi(j,k)      = -Phi(j,k) - 12.d0*(H_0/(c*k_current*a))**2.d0 * Omega_r * Theta(j,2,k) 

          do l = 3, lmax_int
             Theta(j,l,k) = - l/(2.d0*l + 1.d0) * ckH/dtau * Theta(j,l-1,k)
          end do

          ! Store derivatives required for C_l estimation
          call derivs_ytc(x_t(j), y_tight_coupling, dydx)
          
          dv_b(j,k)     = dydx(4) 
          dPhi(j,k)     = dydx(5)
          dTheta(j,0,k) = dydx(6)
          dTheta(j,1,k) = dydx(7) 
          dTheta(j,2,k) = 2.d0/5.d0 * ckH * Theta(j,1,k) - 3.d0/5.d0 * ckH * Theta(j,3,k) &
                          + dtau * 0.9d0 * Theta(j,2,k) 
          dPsi(j,k)     = - dPhi(j,k) + 24.d0 * (H_0/(c*k_current*a))**2.d0 * Omega_r * Theta(j,2,k) &
                          - 12.d0 * (H_0/(c*k_current*a))**2.d0 * Omega_r * dTheta(j,2,k) 
        
          do l = 3, lmax_int - 1
             dTheta(j,l,k) =  l/(2.d0*l + 1.d0) * ckH * Theta(j,l-1,k) &
                    - (l + 1.d0)/(2.d0*l + 1.d0) * ckH * Theta(j,l+1,k) + dtau * Theta(j,l,k) 
          end do

          dTheta(j,lmax_int,k) = ckH * Theta(j,lmax_int-1,k) - c * & 
                           (l + 1.d0)/(get_H_p(x_t(j))*get_eta(x_t(j))) * Theta(j,lmax_int,k) &
                           + dtau * Theta(j,lmax_int,k)

          ! Update j 
          j = j+1
          
       end do


       ! Save the last j-value 
       j_tc = j 


       ! Setting up variables for integration from the end of tight coupling 
       ! until today
       y(1:7) = y_tight_coupling(1:7)
       y(8)   = Theta(j_tc-1,2,k)
       do l = 3, lmax_int
          y(6+l) = Theta(j_tc-1,l,k)
       end do

       
       ! Integrating equations from tight coupling to today

       do i = j_tc, n_t
          a = exp(x_t(i))
          ckH = c*k_current/get_H_p(x_t(i))
          dtau = get_dtau(x_t(i))
                    
          call odeint(y, x_t(i-1), x_t(i), eps, h1, hmin, derivs_y, bsstep, output)

          ! Task: Store variables at time step i in global variables
          delta(i,k)   = y(1)
          delta_b(i,k) = y(2)
          v(i,k)       = y(3)
          v_b(i,k)     = y(4)
          Phi(i,k)     = y(5)
          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l)
          end do
          Psi(i,k)     = -Phi(i,k) - 12.d0*(H_0/(c*k_current*a))**2.d0 * Omega_r * Theta(i,2,k) 

          ! Task: Store derivatives that are required for C_l estimation
          call derivs_y(x_t(i), y, dydx) 

          dPhi(i,k)     = dydx(5)
          dv_b(i,k)     = dydx(4)
          do l= 0, lmax_int
             dTheta(i,l,k) = dydx(6+l)
          end do
          dPsi(i,k)     =  - dPhi(j,k) + 24.d0 * (H_0/(c*k_current*a))**2.d0 * Omega_r * Theta(j,2,k) &
                          - 12.d0 * (H_0/(c*k_current*a))**2.d0 * Omega_r * dTheta(j,2,k) 


       end do

    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)
    
    write(*,*) 'Finished integrating perturbation equations'


  end subroutine integrate_perturbation_eqns




  ! ----------------------------- Subroutines ---------------------------------

  ! Subroutine for calculating the derivatives in the tight-coupling epoch 
  subroutine derivs_ytc(x, y_tc, dydx)
    implicit none 
    
    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y_tc
    real(dp), dimension(:), intent(out) :: dydx

    real(dp) :: delta, delta_b, v, v_b, Phi, Psi, Theta0, Theta1, Theta2
    real(dp) :: ddelta, ddelta_b, dv, dv_b, dPhi, dTheta0, dTheta1 
    real(dp) :: R, q, a, ckH, dtau, ddtau

    a     = exp(x)
    ckH   = c*k_current/get_H_p(x) 
    dtau  = get_dtau(x)
    ddtau = get_ddtau(x)

    ! Extract values from the tight-coupling vector
    delta    = y_tc(1) 
    delta_b  = y_tc(2)
    v        = y_tc(3)
    v_b      = y_tc(4)
    Phi      = y_tc(5)
    Theta0   = y_tc(6) 
    Theta1   = y_tc(7)

    ! Calculate the derivatives 
    R = (4.d0*Omega_r)/(3.d0*Omega_b*a) 
   
    Theta2 = - (20.d0*ckH)/(45.d0*dtau) * Theta1

    Psi =  -Phi - 12.d0*(H_0/(c*k_current*a))**2.d0 * Omega_r * Theta2

    dPhi = Psi - ckH**2.d0 / 3.d0 * Phi + (H_0/get_H_p(x))**2.d0 / 2.d0 & 
           * (Omega_m/a*delta + Omega_b/a*delta_b + 4.d0*Omega_r*Theta0/a**2.d0)

    dTheta0 = - ckH * Theta1 - dPhi

    ddelta = ckH * v - 3.d0*dPhi

    dv = - v - ckH * Psi 

    ddelta_b = ckH * v_b - 3.d0*dPhi

    q = (-((1.d0-2.d0*R)*dtau + (1.d0+R)*ddtau) * (3.d0*Theta1 + v_b) - ckH * Psi & 
        + (1.d0 - get_dH_p(x)/get_H_p(x)) * ckH * (- Theta0 + 2.d0*Theta2) &
        - ckH*dTheta0) / ((1.d0 + R)*dtau + get_dH_p(x)/get_H_p(x) - 1.d0)
 
    dv_b = 1.d0/(1.d0 + R) * (-v_b - ckH * Psi + R*(q + ckH * (-Theta0 + 2.d0*Theta2) - ckH*Psi))

    dTheta1 = 1.d0/3.d0 * (q - dv_b) 

    
    ! Filling in array with all the derivatives 
    dydx(1) = ddelta 
    dydx(2) = ddelta_b 
    dydx(3) = dv 
    dydx(4) = dv_b 
    dydx(5) = dPhi 
    dydx(6) = dTheta0 
    dydx(7) = dTheta1 

  end subroutine derivs_ytc


  ! Subroutine for calculating the derivatives 
  subroutine derivs_y(x, y, dydx) 
    implicit none 

    real(dp),               intent(in)  :: x
    real(dp), dimension(:), intent(in)  :: y
    real(dp), dimension(:), intent(out) :: dydx

    real(dp), dimension(0:6) :: Theta
    real(dp), dimension(0:6) :: dTheta

    real(dp) :: delta, delta_b, v, v_b, Phi, Psi
    real(dp) :: ddelta, ddelta_b, dv, dv_b, dPhi
    real(dp) :: R, a, ckH, dtau

    integer(i4b) :: l

    a     = exp(x)
    ckH   = c*k_current/get_H_p(x) 
    dtau  = get_dtau(x)

    ! Extract values from y_tc vector
    delta    = y(1)
    delta_b  = y(2)
    v        = y(3)
    v_b      = y(4)
    Phi      = y(5)
    Theta(0) = y(6)
    Theta(1) = y(7)
    Theta(2) = y(8)
    Theta(3) = y(9)
    Theta(4) = y(10)
    Theta(5) = y(11)
    Theta(6) = y(12)

    R = (4.d0*Omega_r)/(3.d0*Omega_b*a) 

    Psi =  -Phi - 12.d0*(H_0/(c*k_current*a))**2.d0 * Omega_r * Theta(2)

    dPhi =  Psi - ckH**2.d0 / 3.d0 * Phi + (H_0/get_H_p(x))**2.d0 / 2.d0 & 
           * (Omega_m/a*delta + Omega_b/a*delta_b + 4.d0*Omega_r*Theta(0)/a**2.d0)

    dTheta(0) = - ckH * Theta(1) - dPhi 

    ddelta = ckH*v - 3.d0*dPhi 

    ddelta_b = ckh*v_b - 3.d0*dPhi 

    dv = - v - ckH * Psi 

    dv_b = - v_b - ckH * Psi + dtau*R*(3.d0*Theta(1) + v_b)
    
    dTheta(1) = ckH/3.d0 * Theta(0) - 2.d0/3.d0*ckH*Theta(2) + ckH/3.d0 * Psi + & 
                dtau*(Theta(1) + v_b/3.d0) 

    dTheta(2) = 2.d0/5.d0 * ckH * Theta(1) - 3.d0/5.d0 * ckH * Theta(3) + dtau*0.9d0*Theta(2)

    do l = 3, lmax_int - 1
       dTheta(l) =  l/(2.d0*l + 1.d0) * ckH * Theta(l-1) &
                    - (l + 1.d0)/(2.d0*l + 1.d0) * ckH * Theta(l+1) + dtau * Theta(l) 
    end do
    
    dTheta(lmax_int) = ckH * Theta(lmax_int-1) - c * & 
                           (l + 1.d0)/(get_H_p(x)*get_eta(x)) * Theta(lmax_int) &
                           + dtau * Theta(lmax_int)

    ! Making array with all the derivatives
    dydx(1)  = ddelta
    dydx(2)  = ddelta_b
    dydx(3)  = dv
    dydx(4)  = dv_b 
    dydx(5)  = dPhi
    dydx(6)  = dTheta(0)
    dydx(7)  = dTheta(1)
    dydx(8)  = dTheta(2)
    dydx(9)  = dTheta(3)
    dydx(10) = dTheta(4)
    dydx(11) = dTheta(5)
    dydx(12) = dTheta(6)


  end subroutine derivs_y


  ! Subroutine returning the time at which tight coupling ends.
  ! In this project, we define this as either when dtau < 10 or
  ! c*k/(H_p*dt) > 0.1 or x > x(start of recombination)
  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    integer(i4b)          :: i,n
    real(dp)              :: x, x_start_rec

    n =1d4
    x_start_rec =  -log(1.d0 + 1630.4d0)
 
    do i=0,n
        x = x_init +i*(0.d0-x_init)/n
        if (x < x_start_rec .and. &
            abs(c*k/(get_H_p(x)*get_dtau(x))) <= 0.1d0 .and.& 
            abs(get_dtau(x)) > 10.d0) then 
            get_tight_coupling_time = x
        end if
    end do


  end function get_tight_coupling_time



  ! --------------------------- Write to file -----------------------------

  subroutine write_to_file_evolution_mod
    implicit none

    integer(i4b), dimension(6) :: k_values
    integer(i4b) :: k, j
    
    write(*,*) 'Writing to file'
    
    open(1, file = '../Milestone3/data/perturbations_1.dat', status = 'replace') 
    open(2, file = '../Milestone3/data/perturbations_2.dat', status = 'replace') 
    open(3, file = '../Milestone3/data/perturbations_3.dat', status = 'replace') 
    open(4, file = '../Milestone3/data/perturbations_4.dat', status = 'replace') 
    open(5, file = '../Milestone3/data/perturbations_5.dat', status = 'replace') 
    open(6, file = '../Milestone3/data/perturbations_6.dat', status = 'replace') 
    open(7, file = '../Milestone3/data/k_values.dat', status = 'replace')

    ! k_values = int( 1.d0/6.d0 * (/ n_k, 2*n_k, 3*n_k, 4*n_k, 5*n_k, 6*n_k /))
    k_values = (/ 1, 10, 30, 50, 80, 100 /)

    do k=1,6
       do j=1, n_t
          write(k,"(9(E17.8E3))") x_t(j), Phi(j,k_values(k)), Psi(j,k_values(k)), delta(j,k_values(k)), &
               delta_b(j,k_values(k)), v(j,k_values(k)), v_b(j,k_values(k)), Theta(j,0,k_values(k)), Theta(j,1,k_values(k))
       end do
       
       write(7, "(1(E17.8E3))") ks(k_values(k))*c/H_0

       close(k)

    end do


  end subroutine write_to_file_evolution_mod

end module evolution_mod


