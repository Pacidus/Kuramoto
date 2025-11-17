module mod_kuramoto
    use mod_constants
    implicit none
    
    ! Derived type to store kernel configuration
    type :: kernel_config
        integer :: order                    ! Coupling order (0, 1, 2, ...)
        character(len=256) :: filename      ! Path to kernel binary file
        real(8) :: beta                     ! Phase parameter
        logical :: active                   ! Whether this kernel should be used
    end type kernel_config
    
contains
    !==========================================================================
    subroutine read_config_kernels(config_file, kernel_configs, n_kernels, status)
    !--------------------------------------------------------------------------
    ! Read kernel configuration from config file
    !
    ! Arguments:
    !   config_file     : path to configuration file
    !   kernel_configs  : array of kernel configurations (allocatable)
    !   n_kernels       : number of kernels found
    !   status          : error status (0=success, >0=error)
    !--------------------------------------------------------------------------
        implicit none
        
        ! Arguments
        character(len=*), intent(in) :: config_file
        type(kernel_config), allocatable, intent(out) :: kernel_configs(:)
        integer, intent(out) :: n_kernels, status
        
        ! Local variables
        integer :: unit_num, io_status, line_count, i
        character(len=512) :: line
        character(len=256) :: filename_str
        integer :: order_val
        real(8) :: beta_val
        logical :: in_kernel_section
        type(kernel_config), allocatable :: temp_configs(:)
        
        ! Initialize
        status = 0
        n_kernels = 0
        in_kernel_section = .false.
        
        ! First pass: count kernel entries
        open(newunit=unit_num, file=trim(config_file), status='old', &
             action='read', iostat=io_status)
        
        if (io_status /= 0) then
            print *, 'Error: Unable to open config file: ', trim(config_file)
            status = 1
            return
        end if
        
        do
            read(unit_num, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            
            ! Trim leading spaces
            line = adjustl(line)
            
            ! Skip empty lines
            if (len_trim(line) == 0) cycle
            
            ! Check if we've entered kernel section (before skipping comments)
            if (index(line, '# Kernel configuration') > 0) then
                in_kernel_section = .true.
                cycle
            end if
            
            ! Skip comment lines (but we've already checked for section header)
            if (line(1:1) == '#') cycle
            
            ! If in kernel section and line contains comma, it's a kernel entry
            if (in_kernel_section .and. index(line, ',') > 0) then
                n_kernels = n_kernels + 1
            end if
        end do
        
        close(unit_num)
        
        if (n_kernels == 0) then
            print *, 'Warning: No kernel configurations found in config file'
            status = 2
            return
        end if
        
        ! Allocate array for kernel configs
        allocate(kernel_configs(n_kernels), stat=io_status)
        if (io_status /= 0) then
            print *, 'Error: Unable to allocate memory for kernel configs'
            status = 3
            return
        end if
        
        ! Second pass: read kernel data
        open(newunit=unit_num, file=trim(config_file), status='old', &
             action='read', iostat=io_status)
        
        in_kernel_section = .false.
        line_count = 0
        
        do
            read(unit_num, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            
            line = adjustl(line)
            
            ! Skip empty lines
            if (len_trim(line) == 0) cycle
            
            ! Check for kernel section header
            if (index(line, '# Kernel configuration') > 0) then
                in_kernel_section = .true.
                cycle
            end if
            
            ! Skip other comment lines
            if (line(1:1) == '#') cycle
            
            ! Parse kernel entry
            if (in_kernel_section .and. index(line, ',') > 0) then
                line_count = line_count + 1
                
                ! Parse: order, filename, beta
                call parse_kernel_line(line, order_val, filename_str, beta_val, io_status)
                
                if (io_status /= 0) then
                    print *, 'Warning: Failed to parse kernel line: ', trim(line)
                    cycle
                end if
                
                ! Store configuration
                kernel_configs(line_count)%order = order_val
                kernel_configs(line_count)%filename = trim(filename_str)
                kernel_configs(line_count)%beta = beta_val
                
                ! Check if kernel should be active (filename not "none")
                kernel_configs(line_count)%active = &
                    (trim(adjustl(filename_str)) /= 'none' .and. &
                     trim(adjustl(filename_str)) /= 'NONE')
            end if
        end do
        
        close(unit_num)
        
        print *, 'Successfully loaded', n_kernels, 'kernel configurations'
        
    end subroutine read_config_kernels
    
    !==========================================================================
    subroutine parse_kernel_line(line, order, filename, beta, status)
    !--------------------------------------------------------------------------
    ! Parse a kernel configuration line: order, filename, beta
    !--------------------------------------------------------------------------
        implicit none
        
        character(len=*), intent(in) :: line
        integer, intent(out) :: order, status
        character(len=*), intent(out) :: filename
        real(8), intent(out) :: beta
        
        integer :: comma1, comma2
        character(len=512) :: temp
        
        status = 0
        
        ! Find first comma
        comma1 = index(line, ',')
        if (comma1 == 0) then
            status = 1
            return
        end if
        
        ! Read order
        temp = line(1:comma1-1)
        read(temp, *, iostat=status) order
        if (status /= 0) return
        
        ! Find second comma
        comma2 = index(line(comma1+1:), ',') + comma1
        if (comma2 == comma1) then
            status = 2
            return
        end if
        
        ! Read filename (trim spaces)
        filename = adjustl(line(comma1+1:comma2-1))
        
        ! Read beta
        temp = line(comma2+1:)
        read(temp, *, iostat=status) beta
        
    end subroutine parse_kernel_line

    !==========================================================================
    subroutine read_binary_data(filepath, data, n_rows, n_cols, status)
    !--------------------------------------------------------------------------
    ! Read binary data file with header containing dimensions
    !
    ! Arguments:
    !   filepath : path to binary file
    !   data     : output array (allocatable)
    !   n_rows   : first dimension from header
    !   n_cols   : second dimension from header  
    !   status   : error status (0=success, >0=error)
    !--------------------------------------------------------------------------
        implicit none
        
        ! Arguments
        character(len=*), intent(in) :: filepath
        real(8), allocatable, intent(out) :: data(:,:)
        integer, intent(out) :: n_rows, n_cols
        integer, intent(out) :: status
        
        ! Local variables
        integer :: unit_num, io_status
        integer(kind=4) :: dim1, dim2  ! int32 for header dimensions
        
        ! Initialize status
        status = 0
        
        ! Open the binary file
        open(newunit=unit_num, file=trim(filepath), form='unformatted', &
             access='stream', status='old', action='read', iostat=io_status)
        
        if (io_status /= 0) then
            print *, 'Error: Unable to open file: ', trim(filepath)
            status = 1
            return
        end if
        
        ! Read header (2 int32 values)
        read(unit_num, iostat=io_status) dim1, dim2
        
        if (io_status /= 0) then
            print *, 'Error: Unable to read header from file'
            close(unit_num)
            status = 2
            return
        end if
        
        ! Convert to default integer kind
        n_rows = int(dim1)
        n_cols = int(dim2)
        
        ! Validate dimensions
        if (n_rows <= 0 .or. n_cols <= 0) then
            print *, 'Error: Invalid dimensions in header: ', n_rows, n_cols
            close(unit_num)
            status = 3
            return
        end if
        
        ! Allocate array
        allocate(data(n_rows, n_cols), stat=io_status)
        
        if (io_status /= 0) then
            print *, 'Error: Unable to allocate memory for data array'
            close(unit_num)
            status = 4
            return
        end if
        
        ! Read data (double precision values)
        read(unit_num, iostat=io_status) data
        
        if (io_status /= 0) then
            print *, 'Error: Unable to read data from file'
            deallocate(data)
            close(unit_num)
            status = 5
            return
        end if
        
        ! Close file
        close(unit_num)
        
        print *, 'Successfully loaded data with dimensions: ', n_rows, 'x', n_cols
        
    end subroutine read_binary_data

    !==========================================================================
    subroutine setup_kernels_dynamic(fft_kernels, kernel_configs, n_kernels, &
                                      plan_forward, ny, nx, status)
    !--------------------------------------------------------------------------
    ! Setup convolution kernels for the simulation using config file data
    !
    ! Arguments:
    !   fft_kernels     : output array of FFT-transformed kernels (allocatable)
    !   kernel_configs  : array of kernel configurations
    !   n_kernels       : number of kernels to load
    !   plan_forward    : FFTW forward plan
    !   ny, nx          : grid dimensions
    !   status          : error status (0=success, >0=error)
    !--------------------------------------------------------------------------
        use mod_fft_utils
        
        complex(8), dimension(:,:,:), allocatable, intent(out) :: fft_kernels
        type(kernel_config), dimension(:), intent(in) :: kernel_configs
        integer, intent(in) :: n_kernels, ny, nx
        integer(8), intent(in) :: plan_forward
        integer, intent(out) :: status
        
        real(8), allocatable, dimension(:, :) :: kernel_real
        real(8) :: kernel_sum, scaling_factor
        integer :: n_rows, n_cols, i, io_stat
        
        status = 0
        
        ! Allocate array for all FFT kernels
        allocate(fft_kernels(ny, nx, n_kernels), stat=io_stat)
        if (io_stat /= 0) then
            print *, 'Error: Unable to allocate memory for FFT kernels'
            status = 1
            return
        end if
        
        ! Load and process each kernel
        do i = 1, n_kernels
            if (.not. kernel_configs(i)%active) then
                print *, 'Skipping inactive kernel at order', kernel_configs(i)%order
                fft_kernels(:, :, i) = cmplx(0.0d0, 0.0d0, kind=8)
                cycle
            end if
            
            print *, 'Loading kernel', i, ':', trim(kernel_configs(i)%filename)
            print *, '  Order:', kernel_configs(i)%order, '  Beta:', kernel_configs(i)%beta
            
            ! Load kernel from file
            call read_binary_data(kernel_configs(i)%filename, kernel_real, &
                                 n_rows, n_cols, io_stat)
            
            if (io_stat /= 0) then
                print *, 'Error: Failed to load kernel:', trim(kernel_configs(i)%filename)
                status = 2
                return
            end if
            
            ! Check dimensions
            if (n_rows /= ny .or. n_cols /= nx) then
                print *, 'Error: Kernel dimensions mismatch. Expected:', ny, 'x', nx
                print *, '       Got:', n_rows, 'x', n_cols
                deallocate(kernel_real)
                status = 3
                return
            end if
            
            ! Apply ifftshift
            call ifftshift_real(kernel_real, ny, nx)
            
            ! Calculate scaling factor based on order
            ! Order 0: use k (first order coupling)
            ! Order n > 0: use n (higher order couplings)
            if (kernel_configs(i)%order == 0) then
                scaling_factor = k
            else
                scaling_factor = real(kernel_configs(i)%order, kind=8)
            end if
            
            ! Normalize and scale kernel
            kernel_sum = sum(abs(kernel_real))
            if (kernel_sum > 0.0d0) then
                kernel_real = kernel_real * scaling_factor * dt / kernel_sum
            else
                print *, 'Warning: Kernel sum is zero for', trim(kernel_configs(i)%filename)
            end if
            
            ! Convert to complex and FFT
            fft_kernels(:, :, i) = cmplx(kernel_real, 0.0d0, kind=8)
            call dfftw_execute_dft(plan_forward, fft_kernels(:,:,i), fft_kernels(:,:,i))
            fft_kernels(:, :, i) = conjg(fft_kernels(:, :, i))
            
            ! Clean up
            deallocate(kernel_real)
        end do
        
        print *, 'Successfully loaded', n_kernels, 'kernels'
        
    end subroutine setup_kernels_dynamic

    !==========================================================================
    subroutine setup_kernels(fft_kernel1, fft_kernel2, fft_kernel3, plan_forward)
    !--------------------------------------------------------------------------
    ! Setup convolution kernels for the simulation (LEGACY - kept for compatibility)
    !--------------------------------------------------------------------------
        use mod_fft_utils
        
        complex(8), dimension(ly, lx), intent(out) :: fft_kernel1, fft_kernel2, fft_kernel3
        integer(8), intent(in) :: plan_forward
        
        real(8), allocatable, dimension(:, :) :: kernel_real
        real(8) :: kernel_sum
        integer :: n_rows, n_cols, status

        ! Load and process first kernel
        call read_binary_data("kernel1.bin", kernel_real, n_rows, n_cols, status)
        call ifftshift_real(kernel_real, ly, lx)
        kernel_sum = sum(abs(kernel_real))
        kernel_real = kernel_real * k * dt / kernel_sum
        fft_kernel1 = cmplx(kernel_real, 0.0d0, kind=8)
        call dfftw_execute_dft(plan_forward, fft_kernel1, fft_kernel1)
        fft_kernel1 = conjg(fft_kernel1)
        
        ! Load and process second kernel
        call read_binary_data("kernel2.bin", kernel_real, n_rows, n_cols, status)
        call ifftshift_real(kernel_real, ly, lx)
        kernel_sum = sum(abs(kernel_real))
        kernel_real = kernel_real * dt / kernel_sum
        fft_kernel2 = cmplx(kernel_real, 0.0d0, kind=8)
        call dfftw_execute_dft(plan_forward, fft_kernel2, fft_kernel2)
        fft_kernel2 = conjg(fft_kernel2)
        
        ! Load and process third kernel
        call read_binary_data("kernel2.bin", kernel_real, n_rows, n_cols, status)
        call ifftshift_real(kernel_real, ly, lx)
        kernel_sum = sum(abs(kernel_real))
        kernel_real = kernel_real * 2.0d0 * dt / kernel_sum
        fft_kernel3 = cmplx(kernel_real, 0.0d0, kind=8)
        call dfftw_execute_dft(plan_forward, fft_kernel3, fft_kernel3)
        fft_kernel3 = conjg(fft_kernel3)
        
    end subroutine setup_kernels
    
    !==========================================================================
    subroutine initialize_fields(omega, theta)
    !--------------------------------------------------------------------------
    ! Initialize omega (frequencies) and theta (phases) fields
    !
    ! omega: Gaussian distribution with mean 0.1
    ! theta: Uniform distribution in [0, tau)
    !--------------------------------------------------------------------------
        real(8), dimension(ly, lx), intent(out) :: omega, theta
        real(8) :: rand1, rand2, normal_val, uniform_val
        integer :: i, j
        
        ! Initialize omega with Gaussian distribution
        do j = 1, lx
            do i = 1, ly
                call random_number(rand1)
                call random_number(rand2)
                ! Box-Muller transform for normal distribution
                normal_val = sqrt(-2.0d0 * log(rand1)) * cos(2.0d0 * pi * rand2)
                omega(i, j) = normal_val * dt * tau / 30.0d0 + 0.1d0
            end do
        end do
        
        ! Initialize theta with uniform distribution in [0, tau)
        do j = 1, lx
            do i = 1, ly
                call random_number(uniform_val)
                theta(i, j) = tau * uniform_val
            end do
        end do
        
    end subroutine initialize_fields
    
    !==========================================================================
    subroutine update_theta_dynamic(theta, omega, fft_kernels, kernel_configs, &
                                     n_kernels, plan_forward, plan_backward, ny, nx)
    !--------------------------------------------------------------------------
    ! Update theta field using dynamically loaded kernels
    !
    ! Arguments:
    !   theta, omega      : phase and frequency fields
    !   fft_kernels       : array of precomputed FFT kernels
    !   kernel_configs    : kernel configuration array
    !   n_kernels         : number of kernels
    !   plan_forward      : FFTW forward plan
    !   plan_backward     : FFTW backward plan  
    !   ny, nx            : field dimensions
    !--------------------------------------------------------------------------
        integer, intent(in) :: ny, nx, n_kernels
        real(8), dimension(ny, nx), intent(inout) :: theta
        real(8), dimension(ny, nx), intent(in) :: omega
        complex(8), dimension(ny, nx, n_kernels), intent(in) :: fft_kernels
        type(kernel_config), dimension(n_kernels), intent(in) :: kernel_configs
        integer(8), intent(in) :: plan_forward, plan_backward
        
        complex(8), dimension(ny, nx) :: complex_theta, complex_k, fft_work
        integer :: i, j, k_idx
        complex(8) :: phase_factor
        real(8) :: coupling_order

        ! Add intrinsic frequency
        theta = theta + omega
        
        ! Process each active kernel
        do k_idx = 1, n_kernels
            if (.not. kernel_configs(k_idx)%active) cycle
            
            coupling_order = real(kernel_configs(k_idx)%order, kind=8)
            
            ! Calculate exp(m*i*theta) where m is the coupling order + 1
            do j = 1, nx
                do i = 1, ny
                    complex_theta(i, j) = exp((coupling_order + 1.0d0) * ci * theta(i, j))
                end do
            end do
            
            ! Forward FFT
            fft_work = complex_theta
            call dfftw_execute_dft(plan_forward, fft_work, fft_work)
            
            ! Multiply by kernel in Fourier space
            fft_work = fft_work * fft_kernels(:, :, k_idx)
            
            ! Backward FFT
            call dfftw_execute_dft(plan_backward, fft_work, complex_k)
            complex_k = complex_k / real(nx * ny, kind=8)
            
            ! Apply phase factor (beta) if specified
            if (abs(kernel_configs(k_idx)%beta) > 1.0d-10) then
                phase_factor = exp(ci * kernel_configs(k_idx)%beta)
                complex_k = phase_factor * complex_k
            end if
            
            ! Update theta with coupling term
            ! This is Im(complex_k * conj(complex_theta))
            do j = 1, nx
                do i = 1, ny
                    theta(i, j) = theta(i, j) + &
                                 aimag(complex_k(i, j)) * real(complex_theta(i, j)) - &
                                 real(complex_k(i, j)) * aimag(complex_theta(i, j))
                end do
            end do
        end do
        
        ! Wrap theta to [0, tau)
        do j = 1, nx
            do i = 1, ny
                theta(i, j) = modulo(theta(i, j), tau)
            end do
        end do
        
    end subroutine update_theta_dynamic
    
    !==========================================================================
    subroutine update_theta(theta, omega, fft_kernel1, fft_kernel2, fft_kernel3, &
                           beta, plan_forward, plan_backward, ny, nx, order)
    !--------------------------------------------------------------------------
    ! Update theta field for one time step (LEGACY - kept for compatibility)
    !--------------------------------------------------------------------------
        integer, intent(in) :: ny, nx, order
        real(8), dimension(ny, nx), intent(inout) :: theta
        real(8), dimension(ny, nx), intent(in) :: omega
        complex(8), dimension(ny, nx), intent(in) :: fft_kernel1, fft_kernel2, fft_kernel3
        real(8), dimension(2), intent(in) :: beta
        integer(8), intent(in) :: plan_forward, plan_backward
        
        complex(8), dimension(ny, nx) :: complex_theta, complex_k, fft_work
        integer :: i, j
        complex(8) :: phase_factor

        ! Calculate exp(i*theta) for first order coupling
        do j = 1, nx
            do i = 1, ny
                complex_theta(i, j) = exp(ci * theta(i, j))
            end do
        end do
        
        ! First order coupling
        fft_work = complex_theta
        call dfftw_execute_dft(plan_forward, fft_work, fft_work)
        
        fft_work = fft_work * fft_kernel1
        
        call dfftw_execute_dft(plan_backward, fft_work, complex_k)
        complex_k = complex_k / real(nx * ny, kind=8)
        
        ! Update theta with first order contribution
        do j = 1, nx
            do i = 1, ny
                theta(i, j) = theta(i, j) + omega(i, j) + &
                             aimag(complex_k(i, j)) * real(complex_theta(i, j)) - &
                             real(complex_k(i, j)) * aimag(complex_theta(i, j))
            end do
        end do
        
        ! Higher order coupling (if order >= 1)
        if (order >= 1) then
            phase_factor = exp(ci * beta(1))
            
            ! Calculate exp(2*i*theta) for second order coupling
            do j = 1, nx
                do i = 1, ny
                    complex_theta(i, j) = exp(2.0d0 * ci * theta(i, j))
                end do
            end do
            
            fft_work = complex_theta
            call dfftw_execute_dft(plan_forward, fft_work, fft_work)
            
            fft_work = fft_work * fft_kernel2
            
            call dfftw_execute_dft(plan_backward, fft_work, complex_k)
            complex_k = phase_factor * complex_k / real(nx * ny, kind=8)
            
            ! Update theta with second order contribution
            do j = 1, nx
                do i = 1, ny
                    theta(i, j) = theta(i, j) + &
                                 aimag(complex_k(i, j)) * real(complex_theta(i, j)) - &
                                 real(complex_k(i, j)) * aimag(complex_theta(i, j))
                end do
            end do
        end if
        
        ! Wrap theta to [0, tau)
        do j = 1, nx
            do i = 1, ny
                theta(i, j) = modulo(theta(i, j), tau)
            end do
        end do
        
    end subroutine update_theta
    
    !==========================================================================
    subroutine save_field(theta, iteration)
    !--------------------------------------------------------------------------
    ! Save theta field to binary file
    !--------------------------------------------------------------------------
        real(8), dimension(ly, lx), intent(in) :: theta
        integer, intent(in) :: iteration
        character(len=100) :: filename
        integer :: unit_num
        
        write(filename, '(A,I4.4,A)') 'theta_', iteration, '.dat'
        open(newunit=unit_num, file=filename, form='unformatted', access='stream')
        write(unit_num) theta
        close(unit_num)
        
    end subroutine save_field
    
end module mod_kuramoto
