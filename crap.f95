module functionh
implicit none
contains
	real*16 function h(z) result(t)
		real*16 :: z
		
		t = (1.0 / (z ** z)) - z - 1.0
	end function h
	
	real*16 function dhdz(z) result(t)
		real*16 :: z
		
		t = -((z ** z) * (1.0 + log(z)) + 1.0)
	end function dhdz
	
	real*16 function minusdhdz(z) result(t)
		real*16 :: z
		
		t = (z ** z) * (1.0 + log(z)) + 1.0
	end function minusdhdz
end module functionh

program crap
use functionh
implicit none
! z will be in [e - 1, e]
	real*16 :: z = exp(1.0) - 1.0
	integer :: step
	
	write(*, *) 0, z, 1.0 / z
	do step = 1, 64
!		z = z - h(z) / dhdz(z)
		z = z + h(z) / minusdhdz(z)
		write(*, *) step, z, 1.0 / z
	end do
end program crap

