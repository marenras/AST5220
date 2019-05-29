module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none
    integer(i4b) :: l_num
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: integrand1
    real(dp),     allocatable, dimension(:,:)     :: integrand2
    real(dp),     pointer, dimension(:)           :: x_hires, k_hires, cls_hires, ls_hires
    real(dp),     allocatable, dimension(:,:)     :: Theta_transfer 

contains

  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, k, x_num, n_spline, j_loc
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e, sum_integral, dk
    real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
    real(dp),     pointer,     dimension(:)       :: x
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     pointer,     dimension(:,:)     :: S, S2
    real(dp),     allocatable, dimension(:)       :: z_spline, j_l_spline, j_l_spline2
    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, cls, cls2, ls_dp

    real(dp)           :: t1, t2, integral
    logical(lgt)       :: exist
    character(len=128) :: filename
    real(dp), allocatable, dimension(:) :: y, y2

    ! -------------------- Set up which l's to compute ----------------------
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    

    ! --------------- Get source function from evolution_mod  ---------------   
   
    
    allocate(S(n_hires, n_hires))
    allocate(k_hires(n_hires))
    allocate(x_hires(n_hires))
    
   
    call get_hires_source_function(k_hires, x_hires, S)


    ! ---------------- Initialize spherical bessel functions for each l ---------------

    n_spline = 5400
    allocate(z_spline(n_spline))    ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))
    allocate(j_l2(n_spline, l_num))

    do i=1, n_spline
       z_spline(i) = 3500.d0/(n_spline-1.d0) * (i-1.d0)
    end do 

    ! Calculate Bessel functions 
    do l=1, l_num
       j_l(1,l) = 0.d0
       do i=2, n_spline
          call sphbes(ls(l), z_spline(i), j_l(i,l))
       end do
    end do
    
    ! Spline Bessel functions for each l
    do l=1, l_num
       call spline(z_spline, j_l(:,l), 1.d30, 1.d30, j_l2(:,l))
    end do 
    


   
    ! -------------------- Computing the C_l's for each l -------------------------

    write(*,*) '--------------------------------------------------------------'
    write(*,*) 'Integrating transfer function and power spectrum'

    allocate(integrand1(n_hires))
    allocate(integrand2(l_num, n_hires))
    allocate(Theta_transfer(l_num, n_hires))
    allocate(cls(l_num))
    allocate(cls2(l_num))

    ! Defining step sizes
    dx = x_hires(2)-x_hires(1)
    dk = k_hires(2)-k_hires(1)

    j_loc = locate_dp(k_hires,340.d0*H_0/c)

    do l = 1, l_num
       write(*,*) 'l =', l

       ! Computing the transfer function, Theta_l(k)
       do k=1, n_hires
          sum_integral = 0.d0
          
          do i=1, n_hires
             integrand1(i) = S(i,k)*splint(z_spline, j_l(:,l), j_l2(:,l), k_hires(k)*(get_eta(0.d0) - get_eta(x_hires(i))))  
             !sum_integral = sum_integral + integrand1(i)
          end do


          do i = 1, n_hires
             sum_integral = sum_integral + 0.5d0*(integrand1(i-1) + integrand1(i))*dx
          end do
          

          !Theta_transfer(l,k) = dx * sum_integral
          Theta_transfer(l,k) = sum_integral

       end do 


       ! Integrating P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
       do k=1, n_hires
          integrand2(l,k) = (c*k_hires(k)/H_0)**(n_s-1.d0) * Theta_transfer(l,k)**2.d0 / k_hires(k)
       end do
       
       sum_integral = 0.d0
       do k=2, n_hires
          sum_integral = sum_integral + 0.5d0*(integrand2(l,k-1)+integrand2(l,k))*dk
       end do

       ! Task: Store C_l in an array. Optionally output to file
       cls(l) = sum_integral * ls(l)*(ls(l)+1.d0)/(2.d0*pi)
       
    end do


    ! Splining C_l's found above, and output smooth C_l curve for each integer l

    ! Spline requires double presision, making new array with ls with double presision
    allocate(ls_dp(l_num))

    do j=1, l_num
       ls_dp(j) = ls(j)
    end do 

    call spline(ls_dp, cls, 1.d30, 1.d30, cls2)


    ! Making array with all values of l
    allocate(ls_hires(ls(l_num)))

    do j=1, ls(l_num)
       ls_hires(j) = j 
    end do

    ! Making high resolution array with C_l 
    allocate(cls_hires(ls(l_num)))

    do l=1, ls(l_num)
       cls_hires(l) = splint(ls_dp, cls, cls2, ls_hires(l))
    end do
    
    write(*,*) 'Finished integrating tranfer function and power spectrum'
    write(*,*) '--------------------------------------------------------------'


  end subroutine compute_cls
  



  subroutine write_to_file_cl_mod
    implicit none 

    integer(i4b), dimension(6) :: l_values
    integer(i4b) :: k, l

    write(*,*) 'Writing to files'
    
    open(1, file = '../Milestone4/data/theta_integrand_1.dat', status = 'replace')
    open(2, file = '../Milestone4/data/theta_integrand_2.dat', status = 'replace')
    open(3, file = '../Milestone4/data/theta_integrand_3.dat', status = 'replace')
    open(4, file = '../Milestone4/data/theta_integrand_4.dat', status = 'replace')
    open(5, file = '../Milestone4/data/theta_integrand_5.dat', status = 'replace')
    open(6, file = '../Milestone4/data/theta_integrand_6.dat', status = 'replace')
    open(7, file = '../Milestone4/data/l_values.dat', status = 'replace')
    open(8, file = '../Milestone4/data/powerspectras/CMB_spectrum_best_fit.dat', status = 'replace')

    l_values = (/ 1, 5, 10, 20, 30, 44 /)
    
    ! Writing transfer function and integrand to file for six different l-values
    do l=1, 6
       do k=1, n_hires
          write(l, '(3(E17.8E3))')  k_hires(k)*c/H_0, Theta_transfer(l_values(l), k), integrand2(l_values(l),k) / (c*k_hires(k)/H_0)**(n_s-1.d0)
       end do
       write(7, '(1(I17))') ls(l_values(l))
    end do

    ! Writing c_l values to file 
    do l=1, ls(l_num)
       write(8, '(2(E17.8E3))') ls_hires(l), cls_hires(l)
    end do

  end subroutine write_to_file_cl_mod 


end module cl_mod

