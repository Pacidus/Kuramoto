program kuramoto_simulation_main
    use mod_constants
    use mod_fft_utils
    use mod_kuramoto
    implicit none
    
    ! Simulation arrays
    real(8) :: omega(ly, lx), theta(ly, lx)
    complex(8) :: ct(ly, lx)
    complex(8), allocatable :: fft_kernels(:,:,:)
    
    ! Kernel configuration
    type(kernel_config), allocatable :: kernel_configs(:)
    integer :: n_kernels, status
    
    ! FFT plans and RNG
    integer(8) :: plan_fwd, plan_bwd
    integer :: seed_size, iter
    integer, allocatable :: seed(:)
    
    print *, '========================================='
    print *, 'Kuramoto Simulation with Dynamic Kernels'
    print *, '========================================='
    print *, ''
    
    ! Read kernel configurations from config file
    print *, 'Reading kernel configurations from config.txt...'
    call read_config_kernels('config.txt', kernel_configs, n_kernels, status)
    
    if (status /= 0) then
        print *, 'ERROR: Failed to read config file. Exiting.'
        stop
    end if
    
    print *, ''
    print *, 'Kernel Configuration Summary:'
    print *, '---------------------------------'
    do iter = 1, n_kernels
        print *, 'Kernel', iter, ':'
        print *, '  Order:   ', kernel_configs(iter)%order
        print *, '  File:    ', trim(kernel_configs(iter)%filename)
        print *, '  Beta:    ', kernel_configs(iter)%beta
        print *, '  Active:  ', kernel_configs(iter)%active
        print *, ''
    end do
    
    ! Random seed initialization
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 12345
    call random_seed(put=seed)
    
    ! Setup FFT plans
    print *, 'Setting up FFT plans...'
    call dfftw_plan_dft_2d(plan_fwd, ly, lx, ct, ct, -1, 64)
    call dfftw_plan_dft_2d(plan_bwd, ly, lx, ct, ct, 1, 64)
    
    ! Load and setup kernels
    print *, ''
    print *, 'Loading and processing kernels...'
    call setup_kernels_dynamic(fft_kernels, kernel_configs, n_kernels, &
                               plan_fwd, ly, lx, status)
    
    if (status /= 0) then
        print *, 'ERROR: Failed to setup kernels. Exiting.'
        call dfftw_destroy_plan(plan_fwd)
        call dfftw_destroy_plan(plan_bwd)
        deallocate(seed)
        if (allocated(kernel_configs)) deallocate(kernel_configs)
        if (allocated(fft_kernels)) deallocate(fft_kernels)
        stop
    end if
    
    ! Initialize simulation fields
    print *, ''
    print *, 'Initializing fields...'
    call initialize_fields(omega, theta)
    
    ! Main simulation loop
    print *, ''
    print *, '========================================='
    print *, 'Starting main simulation loop...'
    print *, '========================================='
    print *, ''
    
    do iter = 0, 9999
        ! Save field periodically
        if (mod(iter, 100) == 0) then
            call save_field(theta, iter)
            print *, 'Iteration:', iter, ' - Field saved'
        end if
        
        ! Update theta using dynamic kernels
        call update_theta_dynamic(theta, omega, fft_kernels, kernel_configs, &
                                 n_kernels, plan_fwd, plan_bwd, ly, lx)
    end do
    
    ! Final save
    call save_field(theta, 10000)
    
    print *, ''
    print *, '========================================='
    print *, 'Simulation completed successfully!'
    print *, '========================================='
    
    ! Cleanup
    call dfftw_destroy_plan(plan_fwd)
    call dfftw_destroy_plan(plan_bwd)
    deallocate(seed)
    deallocate(kernel_configs)
    deallocate(fft_kernels)
    
end program kuramoto_simulation_main
