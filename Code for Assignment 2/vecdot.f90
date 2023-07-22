!
!  Function to calculate the scalar product of two vectors, 
!  i.e to calculate d = a*b (parallel version).
!
!------------------------------------------------------------------------
function Vec_Dot(a,b) result(d)
  use header 
  implicit none 
  include "mpif.h"
  type(Vector) :: a,b
  real(kind = 8) :: vec_local,ddot
  real(kind = 8) :: d
  integer :: ierr,n_loc
  n_loc= a%iend - a%ibeg + 1
  !vec_local is the solution to the distributed dot product (unique to processor)
  vec_local = 0.0_8
  ! calculate the local vector product for distributed vectors 'a' and 'b'. 
  vec_local = ddot(n_loc,a%xx(a%ibeg),1,b%xx(b%ibeg),1)
  ! calculate the total global vector dot product result from all processes using 'mpi_sum' and return the result to all
  ! processes 
  call mpi_allreduce(vec_local,d,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr) 


end function Vec_Dot
  
