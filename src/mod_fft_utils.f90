module mod_fft_utils
    implicit none
    
contains
    
    subroutine ifftshift_real(arr, ny, nx)
        ! Performs inverse FFT shift on a real array
        integer, intent(in) :: ny, nx
        real(8), dimension(ny, nx), intent(inout) :: arr
        real(8), dimension(ny, nx) :: temp
        integer :: cy, cx, py, px
        
        cx = nx / 2
        cy = ny / 2
        px = nx - cx
        py = ny - cy
        
        temp = arr
        
        arr(1:py, 1:px) = temp(cy+1:ny, cx+1:nx)
        arr(1:py, px+1:nx) = temp(cy+1:ny, 1:cx)
        arr(py+1:ny, 1:px) = temp(1:cy, cx+1:nx)
        arr(py+1:ny, px+1:nx) = temp(1:cy, 1:cx)
        
    end subroutine ifftshift_real
    
end module mod_fft_utils
