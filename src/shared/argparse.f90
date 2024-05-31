module argparse
  use shared_par
  implicit none
contains

  subroutine argparse_tomo(fname)
    character(len=MAX_STRING_LEN),dimension(:), allocatable :: args
    character(len=MAX_STRING_LEN) :: arg, value
    character(len=MAX_STRING_LEN),intent(out) :: fname
    integer :: i,nargs,m, nopt=2

    nargs = command_argument_count()
    allocate(args(nargs))
    do i = 1, nargs
      call get_command_argument(i, args(i)) 
    enddo
    if (nargs==0 .or. any(args == '-h')) then
      write(*,*)'Usage: surfatt_tomo -i para_file [-h]'
      write(*,*)''
      write(*,*)'Adjoint-state travel time tomography for surface wave'
      write(*,*)''
      write(*,*)'required arguments:'
      write(*,*)' -i para_file  Path to parameter file in yaml format'
      write(*,*)''
      write(*,*)'optional arguments:'
      write(*,*)' -h            Print help message'
    endif
    m = 0
    do i = 1, nargs
      arg = args(i)
      if (arg(1:2) == '-h') then
        stop
      elseif (arg(1:2) == '-i') then
        m = m+1
        fname = args(i+1)
      endif
    enddo
    if (m<1) then
      stop 'not enough arguments'
    elseif(m>nopt) then
      stop 'too many arguments'
    endif
  end subroutine argparse_tomo

end module