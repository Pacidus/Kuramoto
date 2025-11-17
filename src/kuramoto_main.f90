program kuramoto_simulation_main
    use mod_constants
    use mod_fft_utils
    use mod_kuramoto
    implicit none

    !===========================================================================
    ! VARIABLE DECLARATIONS
    !===========================================================================

    ! Simulation fields
    real(8) :: omega(ly, lx)              ! Frequency field
    real(8) :: theta(ly, lx)              ! Phase field
    complex(8) :: ct(ly, lx)              ! Temporary complex array for FFT plans

    ! Kernel data structures
    complex(8), allocatable :: fft_kernels(:, :, :)      ! FFT-transformed kernels
    type(kernel_config), allocatable :: kernel_configs(:)  ! Kernel configurations
    integer :: n_kernels                  ! Number of kernels

    ! Configuration and control
    character(len=256) :: config_file     ! Path to configuration file
    integer :: status                     ! Error status flag

    ! FFT plans
    integer(8) :: plan_fwd                ! Forward FFT plan
    integer(8) :: plan_bwd                ! Backward FFT plan

    ! Random number generator
    integer :: seed_size                  ! Size of RNG seed array
    integer, allocatable :: seed(:)       ! RNG seed

    ! Loop control
    integer :: iter                       ! Iteration counter
    integer :: num_args                   ! Command line argument count

    !===========================================================================
    ! PROGRAM START
    !===========================================================================

    call print_header()

    !---------------------------------------------------------------------------
    ! Parse configuration file
    !---------------------------------------------------------------------------

    config_file = get_config_filename()

    !---------------------------------------------------------------------------
    ! Read and display kernel configurations
    !---------------------------------------------------------------------------

    print *, 'Reading kernel configurations...'
    call read_config_kernels(trim(config_file), kernel_configs, n_kernels, status)

    if (status /= 0) then
        print *, 'ERROR: Failed to read config file. Exiting.'
        stop 1
    end if

    call print_kernel_summary(kernel_configs, n_kernels)

    !---------------------------------------------------------------------------
    ! Initialize random number generator
    !---------------------------------------------------------------------------

    call initialize_rng(seed, seed_size)

    !---------------------------------------------------------------------------
    ! Setup FFT plans
    !---------------------------------------------------------------------------

    print *, 'Setting up FFT plans...'
    call dfftw_plan_dft_2d(plan_fwd, ly, lx, ct, ct, -1, 64)
    call dfftw_plan_dft_2d(plan_bwd, ly, lx, ct, ct, 1, 64)
    print *, 'FFT plans created successfully'

    !---------------------------------------------------------------------------
    ! Load and process kernels
    !---------------------------------------------------------------------------

    print *, ''
    print *, 'Loading and processing kernels...'
    call setup_kernels_dynamic(fft_kernels, kernel_configs, n_kernels, &
                               plan_fwd, ly, lx, status)

    if (status /= 0) then
        print *, 'ERROR: Failed to setup kernels. Exiting.'
        call cleanup_and_exit(plan_fwd, plan_bwd, seed, kernel_configs, &
                              fft_kernels, status_code=1)
    end if

    !---------------------------------------------------------------------------
    ! Initialize simulation fields
    !---------------------------------------------------------------------------

    print *, ''
    print *, 'Initializing fields...'
    call initialize_fields(omega, theta)
    print *, 'Fields initialized successfully'

    !---------------------------------------------------------------------------
    ! Main simulation loop
    !---------------------------------------------------------------------------

    call print_simulation_start()

    do iter = 0, 9999
        ! Save field periodically
        call save_field(theta, iter)

        if (mod(iter, 100) == 0) then
            print *, 'Iteration:', iter, ' - Field saved'
        end if

        ! Update theta using dynamic kernels
        call update_theta_dynamic(theta, omega, fft_kernels, kernel_configs, &
                                  n_kernels, plan_fwd, plan_bwd, ly, lx)
    end do

    ! Final save
    call save_field(theta, 10000)

    !---------------------------------------------------------------------------
    ! Cleanup and exit
    !---------------------------------------------------------------------------

    call print_completion()
    call cleanup_and_exit(plan_fwd, plan_bwd, seed, kernel_configs, &
                          fft_kernels, status_code=0)

contains

    !===========================================================================
    ! INTERNAL SUBROUTINES
    !===========================================================================

    !---------------------------------------------------------------------------
    subroutine print_header()
        !---------------------------------------------------------------------------
        ! Print program header
        !---------------------------------------------------------------------------
        print *, '========================================='
        print *, 'Kuramoto Simulation with Dynamic Kernels'
        print *, '========================================='
        print *, ''
    end subroutine print_header

    !---------------------------------------------------------------------------
    function get_config_filename() result(filename)
        !---------------------------------------------------------------------------
        ! Get configuration filename from command line or use default
        !---------------------------------------------------------------------------
        character(len=256) :: filename

        num_args = command_argument_count()

        if (num_args >= 1) then
            call get_command_argument(1, filename)
            print *, 'Using config file:', trim(filename)
        else
            filename = 'config.txt'
            print *, 'No config file specified, using default: config.txt'
            print *, 'Usage: ./kuramoto_main [config_file]'
        end if

        print *, ''
    end function get_config_filename

    !---------------------------------------------------------------------------
    subroutine print_kernel_summary(configs, n)
        !---------------------------------------------------------------------------
        ! Print summary of kernel configurations
        !---------------------------------------------------------------------------
        type(kernel_config), dimension(:), intent(in) :: configs
        integer, intent(in) :: n
        integer :: i

        print *, ''
        print *, 'Kernel Configuration Summary:'
        print *, '---------------------------------'

        do i = 1, n
            print *, 'Kernel', i, ':'
            print *, '  Order:   ', configs(i)%order
            print *, '  File:    ', trim(configs(i)%filename)
            print *, '  Beta:    ', configs(i)%beta
            print *, '  Active:  ', configs(i)%active
            print *, ''
        end do
    end subroutine print_kernel_summary

    !---------------------------------------------------------------------------
    subroutine initialize_rng(seed_array, size_out)
        !---------------------------------------------------------------------------
        ! Initialize random number generator with fixed seed
        !---------------------------------------------------------------------------
        integer, allocatable, intent(out) :: seed_array(:)
        integer, intent(out) :: size_out

        print *, ''
        print *, 'Initializing random number generator...'

        call random_seed(size=size_out)
        allocate (seed_array(size_out))
        seed_array = 12345
        call random_seed(put=seed_array)

        print *, 'RNG initialized with seed: 12345'
    end subroutine initialize_rng

    !---------------------------------------------------------------------------
    subroutine print_simulation_start()
        !---------------------------------------------------------------------------
        ! Print simulation start banner
        !---------------------------------------------------------------------------
        print *, ''
        print *, '========================================='
        print *, 'Starting main simulation loop...'
        print *, '========================================='
        print *, ''
    end subroutine print_simulation_start

    !---------------------------------------------------------------------------
    subroutine print_completion()
        !---------------------------------------------------------------------------
        ! Print completion message
        !---------------------------------------------------------------------------
        print *, ''
        print *, '========================================='
        print *, 'Simulation completed successfully!'
        print *, '========================================='
        print *, ''
    end subroutine print_completion

    !---------------------------------------------------------------------------
    subroutine cleanup_and_exit(plan_f, plan_b, seed_arr, configs, kernels, status_code)
        !---------------------------------------------------------------------------
        ! Clean up allocated resources and exit
        !---------------------------------------------------------------------------
        integer(8), intent(in) :: plan_f, plan_b
        integer, allocatable, intent(inout) :: seed_arr(:)
        type(kernel_config), allocatable, intent(inout) :: configs(:)
        complex(8), allocatable, intent(inout) :: kernels(:, :, :)
        integer, intent(in) :: status_code

        ! Destroy FFT plans
        call dfftw_destroy_plan(plan_f)
        call dfftw_destroy_plan(plan_b)

        ! Deallocate arrays
        if (allocated(seed_arr)) deallocate (seed_arr)
        if (allocated(configs)) deallocate (configs)
        if (allocated(kernels)) deallocate (kernels)

        ! Exit with appropriate status
        if (status_code /= 0) then
            stop 1
        end if
    end subroutine cleanup_and_exit

end program kuramoto_simulation_main
