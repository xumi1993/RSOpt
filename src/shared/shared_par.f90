module constants
  include "constants.h"

  real(kind=dp), parameter :: PI = 3.141592653589793
  
end module constants

module shared_par
  use constants
  implicit none
  
  integer :: loglevel, LID
  character(len=MAX_STRING_LEN) :: log_fname
end module