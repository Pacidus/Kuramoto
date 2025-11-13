module mod_kuramoto
    use mod_constants
    implicit none
    
contains
    
    subroutine initialize_grid(lx_arr, ly_arr, r)
        ! Initialize spatial grid and distance matrix
        use mod_utils

        real(8), dimension(lx), intent(out) :: lx_arr
        real(8), dimension(ly), intent(out) :: ly_arr
        real(8), dimension(ly, lx), intent(out) :: r
        integer :: i, j
        
        ! Create x grid
        call linspace(-dl * real(lx,8), dl * real(lx,8), lx, lx_arr)
        
        ! Create y grid
        call linspace(-dl * real(ly,8), dl * real(ly,8), ly, ly_arr)
        
        ! Calculate distance matrix r
        do i = 1, lx
            do j = 1, ly
                r(j, i) = ly_arr(j)**2 + lx_arr(i)**2
            end do
        end do
        
        print *, 'Grid check:'
        print *, '  Lx range:', lx_arr(1), 'to', lx_arr(lx)
        print *, '  Ly range:', ly_arr(1), 'to', ly_arr(ly)
        print *, '  R min/max:', minval(r), maxval(r)
        
    end subroutine initialize_grid
    
    
    subroutine setup_kernels(r, fftk1, fftk2, fftk3, plan_forward)
        ! Setup convolution kernels for the simulation
        use mod_fft_utils
        
        real(8), dimension(ly, lx), intent(in) :: r
        complex(8), dimension(ly, lx), intent(out) :: fftk1, fftk2, fftk3
        integer(8), intent(in) :: plan_forward
        
        real(8), dimension(ly, lx) :: kernel_real
        real(8) :: kernel_sum
        integer :: i, j
        
        ! First kernel
        do i = 1, lx
            do j = 1, ly
                kernel_real(j, i) = exp(-r(j, i) / 0.03d0)
            end do
        end do
        call ifftshift_real(kernel_real, ly, lx)
        kernel_sum = sum(abs(kernel_real))
        kernel_real = kernel_real * k * dt / kernel_sum
        
        fftk1 = cmplx(kernel_real, 0.0d0, kind=8)
        call dfftw_execute_dft(plan_forward, fftk1, fftk1)
        fftk1 = conjg(fftk1)
        
        ! Second kernel
        do i = 1, lx
            do j = 1, ly
                kernel_real(j, i) = exp(-r(j, i) / 4.0d0) - 0.2d0 * exp(-r(j, i) / 8.0d0)
            end do
        end do
        call ifftshift_real(kernel_real, ly, lx)
        kernel_sum = sum(abs(kernel_real))
        kernel_real = kernel_real * 1.0d0 * dt / kernel_sum
        
        fftk2 = cmplx(kernel_real, 0.0d0, kind=8)
        call dfftw_execute_dft(plan_forward, fftk2, fftk2)
        fftk2 = conjg(fftk2)
        
        ! Third kernel
        do i = 1, lx
            do j = 1, ly
                kernel_real(j, i) = -exp(-r(j, i) / 10.0d0)
            end do
        end do
        call ifftshift_real(kernel_real, ly, lx)
        kernel_sum = sum(abs(kernel_real))
        kernel_real = kernel_real * 2.0d0 * dt / kernel_sum
        
        fftk3 = cmplx(kernel_real, 0.0d0, kind=8)
        call dfftw_execute_dft(plan_forward, fftk3, fftk3)
        fftk3 = conjg(fftk3)
        
    end subroutine setup_kernels
    
    
    subroutine initialize_fields(omega, theta)
        ! Initialize omega (frequencies) and theta (phases) fields
        real(8), dimension(ly, lx), intent(out) :: omega, theta
        real(8) :: u1, u2, normal_val, x
        integer :: i, j
        
        ! Initialize omega with Gaussian distribution
        do i = 1, lx
            do j = 1, ly
                call random_number(u1)
                call random_number(u2)
                normal_val = sqrt(-2.0d0*log(u1)) * cos(2.0d0*pi*u2)
                omega(j, i) = normal_val * dt * tau / 30.0d0 + 0.1d0
            end do
        end do
        
        ! Initialize theta with uniform distribution
        do i = 1, lx
            do j = 1, ly
                call random_number(x)
                theta(j, i) = tau * x
            end do
        end do
        
    end subroutine initialize_fields
    
    
    subroutine update_theta(theta, omega, fftk1, fftk2, fftk3, beta, &
                           plan_fwd, plan_bwd, ny, nx, n)
        ! Update theta field for one time step
        integer, intent(in) :: ny, nx, n
        real(8), dimension(ny, nx), intent(inout) :: theta
        real(8), dimension(ny, nx), intent(in) :: omega
        complex(8), dimension(ny, nx), intent(in) :: fftk1, fftk2, fftk3
        real(8), dimension(2), intent(in) :: beta
        integer(8), intent(in) :: plan_fwd, plan_bwd
        
        complex(8), dimension(ny, nx) :: ct, ck, fft_work
        integer :: i, j
        complex(8) :: phase_factor
        
        ! Calculate exp(i*theta)
        do j = 1, nx
            do i = 1, ny
                ct(i, j) = exp(ci * theta(i, j))
            end do
        end do
        
        ! First order coupling
        fft_work = ct
        call dfftw_execute_dft(plan_fwd, fft_work, fft_work)
        
        fft_work = fft_work * fftk1
        
        call dfftw_execute_dft(plan_bwd, fft_work, ck)
        ck = ck / (nx * ny)  
        
        ! Update theta with first order contribution
        do j = 1, nx
            do i = 1, ny
                theta(i, j) = theta(i, j) + omega(i, j) + &
                             aimag(ck(i, j)) * real(ct(i, j)) - &
                             real(ck(i, j)) * aimag(ct(i, j))
            end do
        end do
        
        ! Higher order coupling (if n >= 1)
        if (n >= 1) then
            phase_factor = exp(ci * beta(1))
            
            ! Calculate exp(2*i*theta)
            do j = 1, nx
                do i = 1, ny
                    ct(i, j) = exp(2.0d0 * ci * theta(i, j))
                end do
            end do
            
            fft_work = ct
            call dfftw_execute_dft(plan_fwd, fft_work, fft_work)
            
            fft_work = fft_work * fftk2
            
            call dfftw_execute_dft(plan_bwd, fft_work, ck)
            ck = phase_factor * ck / (nx * ny)
            
            ! Update theta with second order contribution
            do j = 1, nx
                do i = 1, ny
                    theta(i, j) = theta(i, j) + &
                                 aimag(ck(i, j)) * real(ct(i, j)) - &
                                 real(ck(i, j)) * aimag(ct(i, j))
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
    
    
    subroutine save_field(theta, iter)
        ! Save theta field to file
        real(8), dimension(ly, lx), intent(in) :: theta
        integer, intent(in) :: iter
        character(len=100) :: filename
        
        write(filename, '(A,I4.4,A)') 'theta_', iter, '.dat'
        open(unit=10, file=filename, form='unformatted', access='stream')
        write(10) theta
        close(10)
        
    end subroutine save_field
    
end module mod_kuramoto
