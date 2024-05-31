program rsopt
  use para, rsp => rs_para_global
  use argparse

  implicit none

  character(len=MAX_STRING_LEN) :: input_fname

  call argparse_tomo(input_fname)

  call rsp%read(input_fname)

end program rsopt