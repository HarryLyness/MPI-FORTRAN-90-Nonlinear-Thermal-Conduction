subroutine negu(u,uflag)
  use header 
  implicit none
  include 'mpif.h'
  ! Subroutine to see to check if the updated solution becomes negative 
  type(vector), intent(in) :: u
  integer, intent(out) :: uflag 
  integer :: ierr
  uflag = 0
  ! check to see if smallest u component is negative
  if (minval(u%xx(u%ibeg:u%iend))< 0)then 
     uflag = 1
  end if 
  ! checks all other processes for negative solutions
  call mpi_allreduce(uflag,uflag,1,mpi_integer,mpi_max,mpi_comm_world,ierr)
  ! if negative solutions are found, then the returned uflag should be greater than 0
end subroutine negu
