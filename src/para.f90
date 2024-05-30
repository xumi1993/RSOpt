module para
  include 'constants.h'
  use commonlib
  use stdlib_logger, only: debug_level, information_level
implicit none

  type, public :: para_data
    character(len=MAX_STRING_LEN) :: data_surf_file, data_rf_file
    character(len=MAX_STRING_LEN), dimension(2) :: data_name = ['surf', 'rf']
    real(kind=dp), dimension(2) :: weight = [1.0_dp, 1.0_dp]
  end type

  type, public :: para_domain
    real(kind=dp) :: zmax
    integer :: dz
  end type

  type, public :: para_output
    character(len=MAX_STRING_LEN) :: output_dir
    logical :: output_syn_data
    integer :: log_level
  end type

  type, public :: para_inversion
    integer :: n_iter, optim_method
    real(kind=dp) :: tol, step_length, max_shrink, sigma
  end type

  type, public :: RSPara
    type(para_data) :: data
    type(para_domain) :: domain
    type(para_output) :: output
    type(para_inversion) :: inversion
    contains
    procedure :: read => read_para_file
  end type

contains

  subroutine read_para_file(self, fname)
    use yaml, only: parse, error_length
    use yaml_types, only: type_node, type_dictionary, type_error, real_kind, &
                        type_list, type_list_item, type_scalar
    class(RSPara), intent(inout) :: self
    character(len=MAX_STRING_LEN), intent(in) :: fname
    character(len=error_length) :: error
    class(type_node), pointer :: root
    class(type_dictionary), pointer :: data_sec, output, domain,inversion
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    type (type_error), pointer :: io_err
    integer :: stat
    character(len=MAX_STRING_LEN) :: errmsg

    root => parse(fname, error = error)
      if (error/='') call exit_mpi(myrank, error)

      select type (root)
      class is (type_dictionary)
        data_sec => root%get_dictionary('data', required=.true., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%data%data_surf_file = data_sec%get_string('data_surf_file', default='', error=io_err)
        this%data%data_rf_file = data_sec%get_string('data_rf_file', default='', error=io_err)
        list => data_sec%get_list('weight', error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        call read_real_list(list, this%data%weight)

        output => root%get('output', required=.true., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%output%output_dir = output%get_string('output_dir', default='', error=io_err)
        this%output%output_syn_data = output%get_logical('output_syn_data', default=.true., error=io_err)

        domain => root%get('domain', required=.true., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%domain%zmax = domain%get_real('zmax', default=0., error=io_err)
        this%domain%dz = domain%get_integer('dz', default=0, error=io_err)

        inversion => root%get('inversion', required=.true., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%inversion%n_iter = inversion%get_integer('n_iter', default=0, error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%inversion%optim_method = inversion%get_integer('optim_method', default=0, error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%inversion%tol = inversion%get_real('tol', default=0., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%inversion%step_length = inversion%get_real('step_length', default=0., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%inversion%max_shrink = inversion%get_real('max_shrink', default=0., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)
        this%inversion%sigma = inversion%get_real('sigma', default=5., error=io_err)
        if (associated(io_err)) call exit_main(io_err%message)

      end select
      call root%finalize()
      deallocate(root)

  end subroutine

  subroutine read_real_list(list, list_out)
    use yaml_types, only: type_scalar, type_list, type_list_item
    class (type_list), pointer :: list
    class (type_list_item), pointer :: item
    real(kind=dp), dimension(:), intent(out) :: list_out
    integer :: i
    
    item => list%first
    i = 1
    do while(associated(item))
      select type (element => item%node)
      class is (type_scalar)
        list_out(i) = element%to_real(0.)
        item => item%next
        i = i + 1
      end select
    enddo
  end subroutine read_real_list

end module