module ConvertUnitsLib

use :: iso_c_binding ! for C/C++ interop
!real(c_double), bind(c) :: degF, degC
implicit none

public DegCtoF

contains

!
! Convert temperature degrees Celsius Fahrenheit
!
subroutine DegCtoF(degC, degF, n)&
    bind(c, name = "DegCtoF")

    integer, intent(in) :: n
    real(c_double), intent(in), dimension(n) :: degC
    real(c_double), intent(out), dimension(n) :: degF
    integer :: i

    do i = 1, n
        degF(i) = ( degC(i) * 1.8 ) + 32
    end do

end subroutine DegCtoF

! End of module
end module ConvertUnitsLib
