module surfker

  use shared_par
  use RayleighWaveKernel
  use utils
  implicit none
  integer, parameter :: iflsph=1
contains

  subroutine depthkernel1d(vel,nz,iwave,igr,mode,kmaxRc,tRc,depz,&
                          sen_vsRc,sen_vpRc,sen_rhoRc)
    integer ::  nz, kmaxRc, iwave,igr, mode
    real(kind=dp), dimension(nz), intent(in) :: vel, depz
    real(kind=dp), dimension(kmaxRc), intent(in) :: tRc
    real(kind=dp), dimension(kmaxRc,nz), intent(inout) :: sen_vsRc, sen_vpRc, sen_rhoRc
    real(kind=dp), dimension(kmaxRc) :: cp,cg
    real, dimension(nz) :: vpz,vsz,rhoz
    real, dimension(:), allocatable :: rvp,rvs,rrho,rthk
    real(kind=dp), dimension(:), allocatable :: dcdar,dcdbr,dcdhr,dcdrr,ur,uz,tr,tz
    integer :: mmax, i

    vsz = real(vel)
    mmax = nz+1
    allocate(dcdar(mmax), dcdbr(mmax), dcdhr(mmax), dcdrr(mmax))
    dcdar = 0._dp
    dcdbr = 0._dp
    dcdhr = 0._dp
    dcdrr = 0._dp
    allocate(ur(mmax), uz(mmax), tr(mmax), tz(mmax))
    allocate(rvp(mmax), rvs(mmax), rrho(mmax), rthk(mmax))
    call get_vprho(vsz, nz, vpz, rhoz)
    call refinegrid(vpz, vsz, rhoz, real(depz), nz, rvp, rvs, rrho, rthk)
    if (iwave == 2 .and. igr == 0) then
      call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,igr,kmaxRc,&
                      tRc,cp)
      do i = 1, kmaxRc
        call sregn96(rthk,rvp,rvs,rrho,mmax,tRc(i),cp(i),cg(i),&
                  ur,uz,tr,tz,dcdar,dcdbr,dcdhr,dcdrr,iflsph)
                  ! thk,vp,vs,rhom,nlayer,&
                  !   t,cp,cg,dispu,dispw,stressu,stressw,&
                  !   dc2da,dc2db,dc2dh,dc2dr,iflsph)
        sen_vsRc(i,1:nz) = dcdbr(1:nz)
        sen_vpRc(i,1:nz) = dcdar(1:nz)
        sen_rhoRc(i,1:nz) = dcdrr(1:nz)
      enddo
    elseif (iwave == 2 .and. igr == 1) then
      block
        real(kind=dp), dimension(kmaxRc) :: t1, t2, cp, c1, c2, cg
        real(kind=dp), dimension(mmax) :: dudar, dudbr, dudhr, dudrr
        real(kind=dp) :: dt = 0.01
        t1 = tRc * (1+0.5*dt)
        t2 = tRc * (1-0.5*dt)
        call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,0,kmaxRc,&
                          tRc,cp)
        call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,0,kmaxRc,&
                          t1,c1)
        call surfdisp96(rthk,rvp,rvs,rrho,mmax,iflsph,iwave,mode,0,kmaxRc,&
                          t2,c2)
        do i = 1, kmaxRc
          call sregnpu(rthk,rvp,rvs,rrho,mmax,tRc(i),cp(i),cg(i),&
                      ur,uz,tr,tz,t1(i),c1(i),t2(i),c2(i),&
                      dcdar,dcdbr,dcdhr,dcdrr,&
                      dudar,dudbr,dudhr,dudrr,iflsph)
                      !  (thk,vp,vs,rhom,nlayer,&
                      !   t,cp,cg,dispu,dispw,stressu,stressw,&
                      !   t1,cp1,t2,cp2,&
                      !   dc2da,dc2db,dc2dh,dc2dr,&
                      !   du2da,du2db,du2dh,du2dr,&
                      !   iflsph)
          sen_vsRc(i,1:nz) = dudbr(1:nz)
          sen_vpRc(i,1:nz) = dudar(1:nz)
          sen_rhoRc(i,1:nz) = dudrr(1:nz)  
        enddo
      end block
    else
      stop 'kernel1D: Only rayleigh wave is supported for now.'
    endif

  end subroutine depthkernel1d

  subroutine fwdsurf1d(vel,iwave,igr,mode,tRc,depz,svel)
    integer, intent(in) :: iwave,igr, mode
    real(kind=dp), dimension(:), intent(in) :: vel, depz, tRc
    real(kind=dp), dimension(:), allocatable, intent(out) :: svel
    real, dimension(:), allocatable :: vpz,vsz,rhoz,rthk
    integer :: mmax, nz, kmaxRc, kk

    nz = size(depz)
    kmaxRc = size(tRc)
    mmax=nz
    vsz = zeros(nz)
    vpz = zeros(nz)
    rhoz = zeros(nz)
    rthk = zeros(nz)
    svel = zeros(kmaxRc)
    vsz(1:nz)=real(vel(1:nz))
    ! some other emperical relationship maybe better, 
    ! This is from Tomas M.Brocher 2005 BSSA
    call get_vprho(vsz, nz, vpz, rhoz)
    do kk=1,mmax-1
      rthk(kk) = depz(kk+1)-depz(kk)
    enddo
    !!half space
    rthk(mmax) = 0.
    ! call refineGrid2LayerMdl(minthk,mmax,depz,vpz,vsz,rhoz,rmax,rdep,&
    !   rvp,rvs,rrho,rthk)
    call surfdisp96(rthk,vpz,vsz,rhoz,mmax,iflsph,iwave,mode,igr,kmaxRc,&
                    tRc,svel)
    ! svel(1:kmaxRc)=cgRc(1:kmaxRc)

  end subroutine fwdsurf1d

  subroutine refinegrid(vp, vs, rho, dep, nz, rvp, rvs, rrho, rthk)
    integer :: nz
    real, dimension(nz), intent(in) :: vp, vs, rho, dep
    real, dimension(nz+1), intent(out) :: rvp, rvs, rrho, rthk
    integer kk, mmax

    mmax = nz+1
    do kk=1,nz
      rvs(kk) = vs(kk)
      rvp(kk) = vp(kk)
      if (kk == nz) then
        rthk(kk) = rthk(kk-1)
      else
        rthk(kk) = dep(kk+1)-dep(kk)
      endif
      rrho(kk) = rho(kk)
    enddo
    rthk(mmax) = 0.
    rvp(mmax) = rvp(nz)
    rvs(mmax) = rvs(nz)
    rrho(mmax) = rrho(nz)

  end subroutine refinegrid

  subroutine get_vprho(vsz, nz, vpz, rhoz)
    integer :: nz
    real, dimension(nz), intent(in) :: vsz
    real, dimension(nz), intent(out) :: vpz, rhoz
    
    vpz=0.9409 + 2.0947*vsz - 0.8206*vsz**2+ &
        0.2683*vsz**3 - 0.0251*vsz**4

    rhoz=1.6612*vpz - 0.4721*vpz**2 + &
        0.0671*vpz**3 - 0.0043*vpz**4 + & 
        0.000106*vpz**5

  end subroutine 
end module surfker
