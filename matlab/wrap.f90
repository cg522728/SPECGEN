module anode_wrap
  use iso_c_binding, only: c_double, c_int
  use anode, only:  anode_int

  implicit none
contains
	subroutine c_anode_int(z, e, v, i, a, iout) bind(c)
		integer(c_int), intent(in) :: z
		real(c_double), intent(in) :: e
		real(c_double), intent(in) :: v
		real(c_double), intent(in) :: i
		real(c_double), intent(in) :: a
		real(c_double), intent(out) :: iout
		iout = anode_int(z, e, v, i, a)
	end subroutine c_anode_int
end module anode_wrap
