module inv

  use para, rsp => rs_para_global
  use data, sd => surf_data_global
  use hdf5_interface
  use surfker
  use setup_att_log

  implicit none

  type, public :: RSInv
    real(kind=dp), dimension(:), allocatable :: init_vs, zgrids, misfits, vsinv
    real(kind=dp) :: maxupdate
    integer :: nz
    contains
    procedure :: init => initialize_inversion, do_inversion
    procedure, private :: steepest_descent, line_search
  end type RSInv

  integer :: itype, iter
  real(kind=dp), dimension(:), allocatable :: gradient_vs
  type(hdf5_file) :: h5f
  character(len=MAX_NAME_LEN) :: module_name='INV'
  character(len=MAX_STRING_LEN) :: msg

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

    final_fname = trim(rsp%output%output_dir)//'/'//trim(final_basename)
    model_fname = trim(rsp%output%output_dir)//'/'//trim(model_basename)
    call h5f%open(model_fname, status='new', action='write')
    call h5f%add('/depth', this%zgrids)
    call h5f%add('/vs_000', this%init_vs)
    call h5f%close()
  end subroutine

  subroutine do_inversion(this)
    class(RSInv), intent(inout) :: this

    real(kind=dp), dimension(:,:), allocatable :: sen_vs, sen_vp, sen_rho
    real(kind=dp), dimension(this%nz) :: update, sen,update_total
    real(kind=dp), dimension(:), allocatable :: syn
    real(kind=dp) :: chi
    integer :: ip, imode

    this%vsinv = this%init_vs
    do iter = 1, rsp%inversion%n_iter
      write(msg, '("Starting ",i3,"th iteration")') iter-1
      call write_log(msg, 1, module_name)
      gradient_vs = zeros(this%nz)
      do itype = 1, 2
        if (.not. sd%vel_type(itype)) cycle

        ! init update_total for different velcity type
        update_total = 0._dp
        do imode = 1, sd%vdata(itype)%nmode

          ! init update for different mode
          update = 0._dp

          ! init matrix of sensitivity kernels
          sen_vs = zeros(sd%vdata(itype)%mdata(imode)%np, this%nz)
          sen_vp = zeros(sd%vdata(itype)%mdata(imode)%np, this%nz)
          sen_rho = zeros(sd%vdata(itype)%mdata(imode)%np, this%nz)

          ! init vector of synthetic data
          syn = zeros(sd%vdata(itype)%mdata(imode)%np)

          ! Calculate synthetic data
          call fwdsurf1d(this%vsinv,sd%iwave,sd%igr(itype),sd%vdata(itype)%mode(imode),&
                        sd%vdata(itype)%mdata(imode)%period,this%zgrids,syn)

          ! calculate misfit 
          chi = 0.5*sum((sd%vdata(itype)%mdata(imode)%velocity-syn)**2)
          this%misfits(iter) = this%misfits(iter) + chi

          ! compute sensitivity kernels
          call depthkernel1d(this%vsinv,this%nz,sd%iwave,sd%igr(itype),&
                            sd%vdata(itype)%mode(imode),sd%vdata(itype)%mdata(imode)%np,&
                            sd%vdata(itype)%mdata(imode)%period,this%zgrids,&
                            sen_vs, sen_vp, sen_rho)

          ! calculate misfit kernels
          do ip = 1, sd%vdata(itype)%mdata(imode)%np
            sen = sen_vs(ip,:) + sen_vp(ip,:)*dalpha_dbeta(this%vsinv) + &
                  sen_rho(ip,:)*drho_dalpha(empirical_vp(this%vsinv))*dalpha_dbeta(this%vsinv)
            update = update + sen * (sd%vdata(itype)%mdata(imode)%velocity(ip)-syn(ip))
          enddo
          update = update / sd%vdata(itype)%mdata(imode)%np

          ! sum misfit kernels of this model to update_total
          update_total = update_total + update
        enddo

        write(msg, '("Total misfit of ",a,": ",F0.4," (",F0.4,"%)")') &
            sd%vel_name(itype), this%misfits(iter), 100*this%misfits(iter)/this%misfits(1)
        call write_log(msg, 1, module_name)

        ! do smooth of kernel to obtain the gradient.
        update_total = smooth_1(update_total, this%zgrids, rsp%inversion%sigma)
        gradient_vs = gradient_vs + update_total * 0.5_dp
      enddo
      
      ! optimization
      if (rsp%inversion%optim_method == 0) then
        call this%steepest_descent()
      else 
        call this%line_search()
      endif
    enddo
    call h5write(final_fname, '/vs', this%vsinv)
  end subroutine
  
  subroutine steepest_descent(this)
    class(RSInv), intent(inout) :: this
    character(len=MAX_NAME_LEN) :: key

    write(key, '("/gradient_",i3.3)') iter-1
    call h5write(model_fname, trim(key), gradient_vs)
    if (iter > 1 .and. this%misfits(iter) > this%misfits(iter-1)) then
      this%maxupdate = this%maxupdate * rsp%inversion%max_shrink
    endif
    gradient_vs = this%maxupdate * gradient_vs/maxval(abs(gradient_vs))
    this%vsinv = this%vsinv * (1+gradient_vs)
    write(key, '("/vs_",i3.3)') iter
    call h5write(model_fname, trim(key), this%vsinv)
  end subroutine

  subroutine line_search(this)
    class(RSInv), intent(inout) :: this
    character(len=MAX_NAME_LEN) :: key
    real(kind=dp), dimension(:), allocatable :: direction, tmp_vs, syn
    integer :: sub_iter, imode
    real(kind=dp) :: chi

    write(key, '("/gradient_",i3.3)') iter-1
    call h5write(model_fname, trim(key), gradient_vs)

    if (iter == 1) then
      direction = gradient_vs
    else
      call get_lbfgs_direction(direction)
    endif

    this%maxupdate = rsp%inversion%step_length
    do sub_iter = 1, rsp%inversion%max_sub_niter
      write(msg, '(a,i3.3,a)') 'Sub-iteration ',sub_iter, ' for line search.'
      call write_log(msg, 0, module_name)
      ! build tmp velocity model
      direction = this%maxupdate * direction/maxval(abs(direction))
      tmp_vs = this%vsinv * (1+direction)

      ! do forward simulation
      chi = 0._dp
      do itype = 1, 2
        if (.not. sd%vel_type(itype)) cycle
        do imode = 1, sd%vdata(itype)%nmode
          call fwdsurf1d(tmp_vs,sd%iwave,sd%igr(itype),sd%vdata(itype)%mode(imode),&
                         sd%vdata(itype)%mdata(imode)%period,this%zgrids,syn)
          chi = chi + 0.5*sum((syn-sd%vdata(itype)%mdata(imode)%velocity)**2)
        enddo
      enddo
      if (chi < this%misfits(iter)) then
        write(msg, '(a,F0.4," is ok, break line search")') 'Step length of ',this%maxupdate
        call write_log(msg, 1, module_name)
        exit
      else
        write(msg, '(a,F0.4,a,F0.4)') 'Misfit ',chi, ' larger than ', this%misfits(iter)
        call write_log(msg, 0, module_name)
        call write_log('step length is too large', 0, module_name)
        this%maxupdate = this%maxupdate * rsp%inversion%max_shrink
      endif
    enddo

    this%vsinv = tmp_vs
    write(key, '("/vs_",i3.3)') iter
    call h5write(model_fname, trim(key), this%vsinv)

  end subroutine

  subroutine get_lbfgs_direction(direction)

    real(kind=dp), dimension(:), allocatable, intent(out) :: direction
    real(kind=dp), dimension(:), allocatable :: gradient0,gradient1,model0,model1,&
                                                q_vector,r_vector
    real(kind=dp), dimension(:,:), allocatable :: gradient_diff, model_diff
    real(kind=dp), dimension(:), allocatable :: p, a
    real(kind=dp) :: p_k_up_sum, p_k_down_sum, p_k, b
    integer :: istore, i, nz, iter_store, nstore, this_iter
    integer, dimension(:), allocatable :: idx_iter

    call write_log('Get L-BFGS direction...', 1, module_name)
    this_iter = iter - 1
    iter_store = this_iter-m_store
    if (iter_store <= iter_start) iter_store = iter_start
    nstore = this_iter-iter_store

    call get_gradient(this_iter, q_vector)
    nz = size(q_vector)
    gradient_diff = zeros(nstore, nz)
    model_diff = zeros(nstore, nz)
    idx_iter = zeros(nstore)
    p = zeros(nstore)
    a = zeros(nstore)
    i = 0
    do istore = this_iter-1,iter_store,-1
      i = i+1
      idx_iter(i) = istore
      call get_gradient(istore, gradient0)
      call get_gradient(istore+1, gradient1)
      call get_model(istore, model0)
      call get_model(istore+1, model1)

      model_diff(i,:) = model1 - model0
      gradient_diff(i,:) = gradient1 - gradient0

      p(i) = 1/sum(model_diff(i,:)*gradient_diff(i,:))
      a(i) = sum(model_diff(i,:)*q_vector)*p(i)
      q_vector = q_vector - a(i)*gradient_diff(i,:)
    enddo
    p_k_up_sum = sum(model_diff(1,:)*gradient_diff(1,:))
    p_k_down_sum = sum(gradient_diff(1,:)*gradient_diff(1,:))
    p_k = p_k_up_sum/p_k_down_sum
    r_vector = p_k*q_vector

    do istore = iter_store,this_iter-1
      i = find_loc(idx_iter, istore)
      b = sum(gradient_diff(i,:)*r_vector)*p(i)
      r_vector = r_vector + model_diff(i,:)*(a(i)-b)
    enddo
    direction = -1.0_dp * r_vector

  end subroutine get_lbfgs_direction

  subroutine get_gradient(it, gradient)
    integer, intent(in) :: it
    real(kind=dp), dimension(:), allocatable, intent(out) :: gradient
    character(len=MAX_NAME_LEN) :: name

    write(name, '("/gradient_",I3.3)') it
    call h5read(model_fname, name, gradient)
  end subroutine get_gradient

  subroutine get_model(it, model)
    integer, intent(in) :: it
    real(kind=dp), dimension(:), allocatable, intent(out) :: model
    character(len=MAX_NAME_LEN) :: name

    write(name, '("/vs_",I3.3)') it
    call h5read(model_fname, name, model)
  end subroutine

end module inv