module mod_kuramoto
    use mod_constants
    implicit none
    
contains
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
    subroutine setup_kernels(fft_kernel1, fft_kernel2, fft_kernel3, plan_forward)
    !--------------------------------------------------------------------------
    ! Setup convolution kernels for the simulation
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
        kernel_real = kernel_real * dt / kernel_sum  ! Note: removed 1.0d0 multiplier
        fft_kernel2 = cmplx(kernel_real, 0.0d0, kind=8)
        call dfftw_execute_dft(plan_forward, fft_kernel2, fft_kernel2)
        fft_kernel2 = conjg(fft_kernel2)
        
        ! Load and process third kernel (same file as kernel2 but different scaling)
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
    subroutine update_theta(theta, omega, fft_kernel1, fft_kernel2, fft_kernel3, &
                           beta, plan_forward, plan_backward, ny, nx, order)
    !--------------------------------------------------------------------------
    ! Update theta field for one time step using Kuramoto model with coupling
    !
    ! Arguments:
    !   theta, omega      : phase and frequency fields
    !   fft_kernel[1-3]   : precomputed FFT kernels
    !   beta              : coupling parameters
    !   plan_forward      : FFTW forward plan
    !   plan_backward     : FFTW backward plan  
    !   ny, nx            : field dimensions
    !   order             : coupling order (0=first order only)
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
