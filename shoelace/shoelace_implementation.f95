module shoelace_implementation
implicit none
contains
	real function shoelace(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy)
	implicit none
		real, intent(in) :: Ax, Ay, Bx, By, Cx, Cy, Dx, Dy
		
		shoelace = 0.5 * ( &
			(Ax * By + Bx * Cy + Cx * Dy + Dx * Ay) - &
			(Bx * Ay + Cx * By + Dx * Cy + Ax * Dy) &
		)
	end function shoelace
	
	logical function cross(X1, Y1, X2, Y2, X3, Y3)
	implicit none
		real, intent(in) :: X1, Y1, X2, Y2, X3, Y3
		
		cross = (Y3 - Y1) * (X2 - X1) > (Y2 - Y1) * (X3 - X1)
	end function cross
	
	logical function intersect(Ax, Ay, Bx, By, Cx, Cy, Dx, Dy)
	implicit none
		real, intent(in) :: Ax, Ay, Bx, By, Cx, Cy, Dx, Dy
		
		intersect = &
			(cross(Ax, Ay, Cx, Cy, Dx, Dy) .neqv. cross(Bx, By, Cx, Cy, Dx, Dy)) .and. &
			(cross(Ax, Ay, Bx, By, Cx, Cy) .neqv. cross(Ax, Ay, Bx, By, Dx, Dy))
	end function intersect
end module shoelace_implementation

