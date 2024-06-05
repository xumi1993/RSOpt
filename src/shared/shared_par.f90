module constants
  include "constants.h"

  real(kind=dp), parameter :: PI = 3.141592653589793
  character(len=MAX_NAME_LEN) :: model_basename = 'model_iter.h5'
  character(len=MAX_NAME_LEN) :: final_basename = 'final_model.h5'
  integer, parameter :: m_store = 5, iter_start = 0
  
end module constants

module shared_par
  use constants
  implicit none
  
  integer :: loglevel, LID
  character(len=MAX_STRING_LEN) :: log_fname, model_fname, final_fname
end module