module mod_constants
    implicit none
    
    ! Mathematical constants
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: tau = 2.0d0 * pi
    
    ! Grid parameters
    real(8), parameter :: dl = 2.0d-3
    real(8), parameter :: dt = 1.0d-1
    integer, parameter :: lx = 2048
    integer, parameter :: ly = 2048
    
    ! Simulation parameters
    integer, parameter :: order_n = 1
    real(8), parameter :: k = 3.0d0
    
    ! Complex unit
    complex(8), parameter :: ci = cmplx(0.0d0, 1.0d0, kind=8)
    
end module mod_constants
