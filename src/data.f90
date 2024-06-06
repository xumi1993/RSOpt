module data
  use para, rsp => rs_para_global
  use commonlib
  use utils

  implicit none

  type, private :: mode_data
    real(kind=dp), dimension(:), allocatable :: period, velocity, syn
    integer :: np
  end type

  type, private :: vel_data
    integer, dimension(:), allocatable :: mode
    type(mode_data), dimension(5) :: mdata
    integer :: nmode
  end type vel_data

  type, public :: SurfData
    type(vel_data), dimension(2) :: vdata
    logical, dimension(2) :: vel_type = [.false., .false.]
    character(len=2), dimension(2) :: vel_name = ['ph', 'gr']
    integer, dimension(2) :: igr = [0, 1]
    integer :: iwave = 2
    contains
    procedure :: read => read_surf_data, write => write_surf_data
  end type SurfData

  type(SurfData), public :: surf_data_global

  contains
  
  subroutine read_surf_data(this)
    class(SurfData), intent(inout) :: this
    real(kind=dp) :: per, vel
    integer :: mode, mode_read
    character(len=MAX_NAME_LEN) :: veltype
    integer :: ierr, i, idx

    open(IIN, file=trim(rsp%data%data_surf_file), status='old',iostat=ierr)
    if (ierr /= 0) call exit_main('cannot open'//trim(rsp%data%data_surf_file))
    
    do 
      read(IIN, *, iostat=ierr) per, vel, veltype, mode_read
      mode = mode_read + 1
      if (mode > 4) call exit_main('The RFOpt only support highest mode of 5')
      if (veltype == 'ph') then
        if ( .not. allocated(this%vdata(1)%mode)) then
          call append(this%vdata(1)%mode, mode)  
        else if ( .not. any(this%vdata(1)%mode == mode)) then
          call append(this%vdata(1)%mode, mode)
        endif
        idx = find_loc(this%vdata(1)%mode, mode)
        call append(this%vdata(1)%mdata(idx)%period, per)
        call append(this%vdata(1)%mdata(idx)%velocity, vel)
      elseif (veltype == 'gr') then
        if ( .not. allocated(this%vdata(2)%mode) ) then
          call append(this%vdata(2)%mode, mode)  
        else if ( .not. any(this%vdata(2)%mode == mode)) then
          call append(this%vdata(2)%mode, mode)
        endif
        idx = find_loc(this%vdata(2)%mode, mode)
        call append(this%vdata(2)%mdata(idx)%period, per)
        call append(this%vdata(2)%mdata(idx)%velocity, vel)
      endif
      if (ierr /=0) exit
    enddo
    do i = 1, 2
      if ( allocated(this%vdata(i)%mode) ) then
        this%vel_type(i) = .true.
        this%vdata(i)%nmode = size(this%vdata(i)%mode)
        do idx = 1, this%vdata(i)%nmode
          this%vdata(i)%mdata(idx)%np = size(this%vdata(i)%mdata(idx)%period)
        enddo
      endif
    enddo

  end subroutine read_surf_data

  subroutine write_surf_data(this, filename)
    class(SurfData), intent(in) :: this
    integer :: i, j, k
    integer :: ierr
    character(len=*), intent(in) :: filename
    integer :: nmode, np
    real(kind=dp) :: per, vel
    integer :: mode
    character(len=2) :: veltype

    open(IOUT, file=filename, status='unknown',iostat=ierr)
    if (ierr /= 0) call exit_main('cannot open'//trim(filename))

    do i = 1, 2
      if (this%vel_type(i)) then
        nmode = this%vdata(i)%nmode
        do j = 1, nmode
          mode = this%vdata(i)%mode(j)
          np = this%vdata(i)%mdata(j)%np
          do k = 1, np
            per = this%vdata(i)%mdata(j)%period(k)
            vel = this%vdata(i)%mdata(j)%syn(k)
            veltype = this%vel_name(i)
            write(IOUT, *) per, vel, veltype, mode-1
          enddo
        enddo
      endif
    enddo
    close(IOUT)
  end subroutine write_surf_data

end module data