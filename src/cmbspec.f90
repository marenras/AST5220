program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  use cl_mod
  implicit none


  ! Varibles for time estimation 
  real :: start, finish 

  call cpu_time(start)

  ! Initialize time grids
  call initialize_time_mod
  call initialize_rec_mod
  call initialize_perturbation_eqns
  call integrate_perturbation_eqns  
  call compute_cls
  call write_to_file_cl_mod


  call cpu_time(finish)
  write(*,*) 'Finished. Time = ', finish-start


end program cmbspec
