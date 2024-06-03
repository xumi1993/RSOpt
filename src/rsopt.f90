program rsopt
  use para, rsp => rs_para_global
  use data, sd => surf_data_global
  use inv
  use argparse
  use setup_att_log

  implicit none

  type(RSInv) :: rsi
  character(len=MAX_STRING_LEN) :: input_fname

  call argparse_tomo(input_fname)

  call rsp%read(input_fname)

  call setuplog(rsp%output%log_level)

  call sd%read()

  call rsi%init()

  call rsi%do_inversion()

end program rsopt