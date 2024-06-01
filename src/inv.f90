module inv

  use para, rsp => rs_para_global
  use data, sd => surf_data_global
  use hdf5_interface
  use surfker

  implicit none

  type, public :: RSInv
    real(kind=dp), dimension(:), allocatable :: init_vs, zgrids, misfits, vsinv
    real(kind=dp) :: maxupdate
    integer :: nz
    contains
    procedure :: init => initialize_inversion
    procedure, private :: steepest_descent
  end type RSInv

  integer :: itype, iter
  real(kind=dp), dimension(:), allocatable :: gradient_vs
  type(hdf5_file) :: h5f

  contains

  subroutine initialize_inversion(this)
    class(RSInv), intent(inout) :: this
    
    this%zgrids = arange(0._dp, rsp%domain%zmax, rsp%domain%dz)
    this%nz = size(this%zgrids)
    this%init_vs = linspace(rsp%domain%init_vel_range(1), &
                            rsp%domain%init_vel_range(2), &
                            this%nz)
    this%misfits = zeros(rsp%inversion%n_iter)
    this%maxupdate = rsp%inversion%step_length

    model_fname = trim(rsp%output%output_dir)//'/'//trim(model_basename)
    call h5f%open(model_fname, status='new', action='write')
    call h5f%add('/depth', this%zgrids)
    call h5f%add('/vs_000', this%init_vs)
  end subroutine

  subroutine do_inversion(this)
    class(RSInv), intent(inout) :: this

    real(kind=dp), dimension(:,:), allocatable :: sen_vs, sen_vp, sen_rho
    real(kind=dp), dimension(this%nz) :: update, sen,update_total
    real(kind=dp), dimension(:), allocatable :: tmp
    real(kind=dp) :: chi
    integer :: ip

    this%vsinv = this%init_vs
    do iter = 1, rsp%inversion%n_iter
      update = 0._dp
      gradient_vs = zeros(this%nz)
      do itype = 1, 2
        if (.not. sd%vel_type(itype)) cycle
        sen_vs = zeros(sd%vdata(itype)%np, this%nz)
        sen_vp = zeros(sd%vdata(itype)%np, this%nz)
        sen_rho = zeros(sd%vdata(itype)%np, this%nz)
        tmp = zeros(sd%vdata(itype)%np)
        call fwdsurf1d(this%vsinv,sd%iwave,sd%igr(itype),&
                       sd%vdata(itype)%period,this%zgrids,tmp)
        chi = 0.5*sum((sd%vdata(itype)%velocity-tmp)**2)
        this%misfits(iter) = this%misfits(iter) + chi

        call depthkernel1d(this%vsinv,this%nz,sd%iwave,&
                            sd%igr(itype),sd%vdata(itype)%np,&
                            sd%vdata(itype)%period,this%zgrids,&
                            sen_vs, sen_vp, sen_rho)
        update = 0.
        do ip = 1, sd%vdata(itype)%np
          sen = sen_vs(ip,:) + sen_vp(ip,:)*dalpha_dbeta(this%vsinv) + &
                sen_rho(ip,:)*drho_dalpha(empirical_vp(this%vsinv))*dalpha_dbeta(this%vsinv)
          update = update - sen * (sd%vdata(itype)%velocity(ip)-tmp(ip))
        enddo
        update = update / sd%vdata(itype)%np
        update = smooth_1(update, this%zgrids, rsp%inversion%sigma)
        gradient_vs = gradient_vs + update * 0.5_dp
      enddo
      if (rsp%inversion%optim_method == 0) then
        call this%steepest_descent()
      endif
    enddo

  end subroutine
  
  subroutine steepest_descent(this)
    class(RSInv), intent(inout) :: this
    character(len=MAX_NAME_LEN) :: key

    write(key, '("/gradient_",i3)') iter
    call h5f%add(trim(key), gradient_vs)
    if (iter > 1 .and. this%misfits(iter) > this%misfits(iter-1)) then
      this%maxupdate = this%maxupdate * rsp%inversion%max_shrink
    endif
    gradient_vs = this%maxupdate * gradient_vs/maxval(abs(gradient_vs))
    this%vsinv = this%vsinv * (1-gradient_vs)
  end subroutine
end module inv