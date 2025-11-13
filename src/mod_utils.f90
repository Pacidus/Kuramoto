module mod_utils
    implicit none
    
contains
    subroutine linspace(xmin, xmax, n, arr)
        real(8), intent(in) :: xmin, xmax
        integer, intent(in) :: n
        real(8), dimension(n), intent(out) :: arr
        integer :: i
        real(8) :: dl

        dl = (xmax - xmin)/real(n, 8)

        do i = 1, n
            arr(i) = xmin + (real(i, 8) * dl) 
        end do
    end subroutine linspace 
end module mod_utils
