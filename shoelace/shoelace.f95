program main
use shoelace_implementation
implicit none
	integer :: perm(4, 24) = reshape( &
		(/ &
			1, 2, 3, 4, &
			1, 2, 4, 3, &
			1, 3, 2, 4, &
			1, 3, 4, 2, &
			1, 4, 2, 3, &
			1, 4, 3, 2, &
			2, 1, 3, 4, &
			2, 1, 4, 3, &
			2, 3, 1, 4, &
			2, 3, 4, 1, &
			2, 4, 1, 3, &
			2, 4, 3, 1, &
			3, 1, 2, 4, &
			3, 1, 4, 2, &
			3, 2, 1, 4, &
			3, 2, 4, 1, &
			3, 4, 1, 2, &
			3, 4, 2, 1, &
			4, 1, 2, 3, &
			4, 1, 3, 2, &
			4, 2, 1, 3, &
			4, 2, 3, 1, &
			4, 3, 1, 2, &
			4, 3, 2, 1 &
		/), (/ 4, 24 /) &
	)
	
	integer :: row
	real :: temp(2, 4)
		
	read (*, *) &
		temp(1, 1), temp(2, 1), &
		temp(1, 2), temp(2, 2), &
		temp(1, 3), temp(2, 3), &
		temp(1, 4), temp(2, 4)
	
	do row = 1, 24
		write (*, *) &
			intersect( &
				temp(1, perm(1, row)), temp(2, perm(1, row)), &
				temp(1, perm(2, row)), temp(2, perm(2, row)), &
				temp(1, perm(3, row)), temp(2, perm(3, row)), &
				temp(1, perm(4, row)), temp(2, perm(4, row))), &
			shoelace( &
				temp(1, perm(1, row)), temp(2, perm(1, row)), &
				temp(1, perm(2, row)), temp(2, perm(2, row)), &
				temp(1, perm(3, row)), temp(2, perm(3, row)), &
				temp(1, perm(4, row)), temp(2, perm(4, row)))
	end do	
end program main

