program main
use shoelace_implementation
implicit none
	integer :: perm(24, 4) = reshape( &
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
		/), (/ 24, 4 /), order = (/ 2, 1 /) &
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
				temp(1, perm(row, 1)), temp(2, perm(row, 1)), &
				temp(1, perm(row, 2)), temp(2, perm(row, 2)), &
				temp(1, perm(row, 3)), temp(2, perm(row, 3)), &
				temp(1, perm(row, 4)), temp(2, perm(row, 4))), &
			shoelace( &
				temp(1, perm(row, 1)), temp(2, perm(row, 1)), &
				temp(1, perm(row, 2)), temp(2, perm(row, 2)), &
				temp(1, perm(row, 3)), temp(2, perm(row, 3)), &
				temp(1, perm(row, 4)), temp(2, perm(row, 4)))
	end do	
end program main

