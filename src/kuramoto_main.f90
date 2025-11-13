program kuramoto_simulation_main
    use mod_constants
    use mod_fft_utils
    use mod_kuramoto
    implicit none
    
    ! Arrays for simulation
    real(8), dimension(ly, lx) :: omega, theta, r
    real(8), dimension(ly) :: ly_arr
    real(8), dimension(lx) :: lx_arr
    complex(8), dimension(ly, lx) :: ct  ! Needed for FFT plan creation
    complex(8), dimension(ly, lx) :: fftk1, fftk2, fftk3
    real(8), dimension(2) :: beta
    
    ! FFT plans
    integer(8) :: plan_forward, plan_backward
    
    ! Random number seed
    integer :: seed_size
    integer, allocatable :: seed(:)
    
    ! Loop counter
    integer :: iter
    
    ! Initialize beta parameters
    beta(1) = -pi / 2.0d0
    beta(2) = -pi / 4.0d0
    
    ! Initialize random seed
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345
    call random_seed(put=seed)
    
    ! Initialize spatial grid
    call initialize_grid(lx_arr, ly_arr, r)
    
    ! Create FFT plans
    call dfftw_plan_dft_2d(plan_forward, ly, lx, ct, ct, -1, 64)
    call dfftw_plan_dft_2d(plan_backward, ly, lx, ct, ct, 1, 64)
    
    ! Setup convolution kernels
    call setup_kernels(r, fftk1, fftk2, fftk3, plan_forward)
    
    ! Initialize fields
    call initialize_fields(omega, theta)
    
    ! Main simulation loop
    do iter = 0, 9999
        ! Save current state
        call save_field(theta, iter)
        
        ! Progress indicator
        if (mod(iter, 100) == 0) then
            print *, 'Iteration:', iter
        end if
        
        ! Update theta field
        call update_theta(theta, omega, fftk1, fftk2, fftk3, beta, &
                         plan_forward, plan_backward, ly, lx, order_n)
    end do
    
    ! Clean up FFT plans
    call dfftw_destroy_plan(plan_forward)
    call dfftw_destroy_plan(plan_backward)
    
    ! Deallocate seed array
    deallocate(seed)
    
end program kuramoto_simulation_main
