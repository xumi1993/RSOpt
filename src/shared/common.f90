module commonlib
  include 'constants.h'

contains

  subroutine exit_main(msg)
    character(len=*), intent(in) :: msg
    print *, trim(msg)
    stop
  end subroutine exit_main
end module commonlib