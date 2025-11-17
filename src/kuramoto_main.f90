program kuramoto_simulation_main
    use mod_constants
    use mod_fft_utils
    use mod_kuramoto
    implicit none
    
    ! Simulation arrays
    real(8) :: omega(ly, lx), theta(ly, lx)
    complex(8) :: ct(ly, lx), fftk1(ly, lx), fftk2(ly, lx), fftk3(ly, lx)
    real(8) :: beta(2)
    
    ! FFT plans and RNG
    integer(8) :: plan_fwd, plan_bwd
    integer :: seed_size, iter
    integer, allocatable :: seed(:)
    
    ! Initialize parameters
    beta = [-pi/2.0d0, -pi/4.0d0]
    
    ! Random seed
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345
    call random_seed(put=seed)
    
    ! Setup domain and FFT
    call dfftw_plan_dft_2d(plan_fwd, ly, lx, ct, ct, -1, 64)
    call dfftw_plan_dft_2d(plan_bwd, ly, lx, ct, ct, 1, 64)
    call setup_kernels(fftk1, fftk2, fftk3, plan_fwd)
    
    ! Initialize simulation
    call initialize_fields(omega, theta)
    
    ! Main loop
    do iter = 0, 9999
        call save_field(theta, iter)
        if (mod(iter, 100) == 0) print *, 'Iteration:', iter
        call update_theta(theta, omega, fftk1, fftk2, fftk3, beta, &
                         plan_fwd, plan_bwd, ly, lx, order_n)
    end do
    
    ! Cleanup
    call dfftw_destroy_plan(plan_fwd)
    call dfftw_destroy_plan(plan_bwd)
    deallocate(seed)
    
end program kuramoto_simulation_main
