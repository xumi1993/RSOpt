module data
  use para, rsp => rs_para_global
  use commonlib
  use utils

  implicit none

  type, private :: vel_data
    real(kind=dp), dimension(:), allocatable :: period, velocity
    integer :: np
  end type vel_data

  type, public :: SurfData
    type(vel_data), dimension(2) :: vdata
    logical, dimension(2) :: vel_type = [.false., .false.]
    integer, dimension(2) :: igr = [0, 1]
    integer :: iwave = 2
    contains
    procedure :: read => read_surf_data
  end type SurfData

  type(SurfData), public :: surf_data_global

  contains
  
  subroutine read_surf_data(this)
    class(SurfData), intent(inout) :: this
    real(kind=dp) :: per, vel
    character(len=MAX_NAME_LEN) :: veltype
    integer :: ierr, i

    open(IIN, file=trim(rsp%data%data_surf_file), status='old',iostat=ierr)
    if (ierr /= 0) call exit_main('cannot open'//trim(rsp%data%data_surf_file))

    do 
      read(IIN, *, iostat=ierr) per, vel, veltype
      if (veltype == 'ph') then
        call append(this%vdata(1)%period, per)
        call append(this%vdata(1)%velocity, vel)
      elseif (veltype == 'gr') then
        call append(this%vdata(2)%period, per)
        call append(this%vdata(2)%velocity, vel)
      endif
      if (ierr /=0) exit
    enddo
    do i = 1, 2
      if (size(this%vdata(i)%period) /= 0) this%vel_type(i) = .true.
    enddo

  end subroutine read_surf_data

end module data