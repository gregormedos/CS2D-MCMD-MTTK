module parameters
  !ensemble parameters
  integer :: N                  !number of particles (N)
  real*8 :: dens                !density
  real*8 :: temp                !temperature
  !simulation parameters
  real*8 :: lbox                !length of box
  real*8 :: a1,a2               !LJ potential parameters
  real*8 :: a3,a4,a5            !CS potential parameters
  integer :: ran                !random seed
  integer :: nseries            !number of series of sampling (MAX 20)
  integer :: nekv               !equilibration length
  !MC parameters
  integer :: ncycl              !number of cycles = N trial moves
  integer :: nsampl             !number of trial moves per sample (same as N)
  real*8 :: delta               !MAX random displacement by sigma
  !MD parameters
  integer :: nsteps             !number of time steps
  real*8 :: tstep,dt2,d2t2      !time step, half step, half squared step
  real*8 :: dt4,dt8             !quarter step, eighth step
  real*8 :: tt                  !temperature time constant
  real*8 :: NfkT                !number of degrees of freedom * temperature
  real*8, dimension(10) :: Q    !heat bath coupling masses (10)
  integer :: nframe             !number of time steps per frame for movie
  !number Pi
  real*8 :: Pi
  PARAMETER (Pi=4.0d0*datan(1.0d0))
endmodule parameters

module particles
  real*8, dimension(:,:), allocatable :: pos               !matrix of position vectors
  real*8, dimension(:,:), allocatable :: vel               !matrix of velocity vectors
  real*8, dimension(:,:), allocatable :: acc               !matrix of acceleration vectors
  real*8 :: ekin                                           !kinetic energy
  real*8 :: tem                                            !equipartion theorem temperature
  real*8 :: etot                                           !total internal energy
  real*8, dimension(10) :: vksi                            !heat bath velocities
  real*8, dimension(10) :: G                               !heat bath forces
endmodule particles

module potential
  real*8, dimension(:,:), allocatable :: pot               !matrix of potential energies
  real*8, dimension(:,:), allocatable :: tla               !matrix of virial pressures
  real*8, dimension(:,:), allocatable :: dis               !matrix of distances
  real*8, dimension(:,:,:), allocatable :: rrr             !matrix of distance vectors
  real*8 :: rcut                                           !distance of potential energy cutoff
endmodule potential

module mc
  real*8 :: u,uu,uu2                !potential energy, excess internal energy, excess internal energy squared
  real*8 :: p,pp                    !virial pressure
  integer :: m                      !number of samples
  integer :: accx                   !number of accepted trial moves for displacement
  real*8 :: kx                      !ratio of accepted trial moves for displacement
endmodule mc

module correlation
  real*8, dimension(10000) :: gr                !correlation function g(r)
  real*8 :: grinterv                            !g(r) radial interval
  integer :: ngr                                !cycles/steps per g(r)
  integer :: mmm                                !number of samples
  real*8 :: top                                 !translational order parameter
endmodule correlation

module diffusion
  real*8, dimension(:), allocatable :: msd               !mean squared displacement msd(n) or msd(t)
  integer :: nmsd                                        !every n steps for msd(n) or msd(t)
  integer :: nmsdpoints                                  !number of msd(n) or msd(t) points
  real*8 :: msdlin(10)                                   !linear regression of msd(n) or msd(t) with intercept term
  real*8 :: msdlinzero(8)                                !linear regression of msd(n) or msd(t) without intercept term
  real*8, dimension(:,:), allocatable :: intpos, difpos  !initial positions and current positions inside or outside box
endmodule diffusion


program core
  use parameters
  use particles
  use potential
  use mc
  use correlation
  use diffusion
  implicit none
  real*8 :: time1,time2                         !CPU_time
  integer, dimension(8) :: values1,values2      !value(1)=year, value(2)=month, value(3)=day, value(4)=time difference with UTC in minutes, value(5)=hour, value(6)=minute, value(7)=second, value(8)=milisecond
  character :: tekst
  character(8) :: fmt
  integer :: i,j
  real*8 :: cv                                  !heat capacity
  real*8 :: therm(25),dtherm(25),thermo(25,20)  !thermodynamic quantities
  real*8 :: dgr(10000),ggr(10000,20)            !correlation function
  real*8, allocatable :: mmsd(:,:), dmsd(:)     !mean squared displacement

  call cpu_time(time1)
  call date_and_time(VALUES=values1)
  open(1002,file='!md.log')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'start date and time YYYY, MM, DD, UTC, HH, MIN, SS'
  write(1002,'(8i5)')values1

  open(1001,file='md_input')
  read(1001,*)tekst
  read(1001,*)tekst,N
  read(1001,*)tekst,dens
  read(1001,*)tekst,temp
  read(1001,*)tekst
  read(1001,*)tekst,a1,a2
  read(1001,*)tekst,a3,a4,a5
  read(1001,*)tekst,ran
  read(1001,*)tekst,nseries
  read(1001,*)tekst,nekv
  read(1001,*)tekst
  read(1001,*)tekst,ncycl
  read(1001,*)tekst,nsampl
  read(1001,*)tekst,delta
  read(1001,*)tekst
  read(1001,*)tekst,nsteps
  read(1001,*)tekst,tstep
  read(1001,*)tekst,tt
  read(1001,*)tekst,nframe
  read(1001,*)tekst
  read(1001,*)tekst,grinterv
  read(1001,*)tekst,ngr
  read(1001,*)tekst,nmsd
  close(1001)

  allocate( pos(N,2) )
  allocate( vel(N,2) )
  allocate( acc(N,2) )
  allocate( pot(N,N) )
  allocate( tla(N,N) )
  allocate( dis(N,N) )
  allocate( rrr(N,N,2) )
  lbox=dsqrt(dble(N)/dens)
  rcut=5.0d0*a2              !CUTOFF
  a4=1.0d0/(a2*a4)
  a4=a4**2
  dt2=tstep/2.0d0
  d2t2=tstep**2/2.0d0
  dt4=tstep/4.0d0
  dt8=tstep/8.0d0
  NfkT=2.0d0*dble(N)*temp
  Q(1)=NfkT*tt**2
  do i=2,10
    Q(i)=temp*tt**2
  enddo
  vksi=0.0d0
  G(1)=0.0d0
  do i=1,9
    G(i+1)=(Q(i)*vksi(i)**2-temp)/Q(i+1)
  enddo

  write(1002,*)'N=',N
  write(1002,*)'V=',lbox**2
  write(1002,*)'T=',temp
  write(1002,*)'density=',dens
  write(1002,*)'length of box=',lbox
  write(1002,*)'random seed=',ran
  write(1002,*)'MAX random displacement by sigma=',delta
  write(1002,*)'time step=',tstep
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'Monte Carlo (MC) and Molecular Dynamics (MD) simulation of 2D Core-Softened (CS) discs'

  !----------------------------------------------------------
  !MC
  !----------------------------------------------------------
  nmsdpoints=ncycl/nmsd
  allocate( msd(nmsdpoints) )
  allocate( mmsd(nmsdpoints,20) )
  allocate( dmsd(nmsdpoints) )
  allocate( intpos(N,2) )
  allocate( difpos(N,2) )

  call random_pos()
  call snapshot('sta')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'particles randomly put in box'

  call minimize()
  call snapshot('min')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'bad contacts removed'
  write(1002,*)'acceptance ratio kx=',kx

  call mcequil()
  call snapshot('ekc')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MC equilibration successful'
  write(1002,*)'uu/N=',uu/dble(N)
  write(1002,*)'pp=',pp
  write(1002,*)'acceptance ratio kx=',kx

  open(223,file='!results')
  write(223,*)'MC'
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MC simulation initialized'
  do j=1,nseries
    call mcseries(j)
    write(fmt,'(a1,i2.2)')'c',j
    call snapshot(trim(fmt))
    cv=(uu2-uu**2)/temp**2
    thermo(1,j)=uu/dble(N)
    thermo(2,j)=pp
    thermo(3,j)=cv/dble(N)
    thermo(4,j)=top
    thermo(5,j)=msdlin(1)/4.0d0
    thermo(6,j)=msdlinzero(1)/4.0d0
    ggr(:,j)=gr(:)
    mmsd(:,j)=msd(:)
    write(223,'(8e16.7)')temp,dens,uu/dble(N),pp,cv/dble(N),top,msdlin(1)/4.0d0,msdlinzero(1)/4.0d0
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*) j,'-th MC run successful'
    write(1002,*)'uu/N=',uu/dble(N)
    write(1002,*)'pp=',pp
    write(1002,*)'cv/N=',cv/dble(N)
    write(1002,*)'top=',top
    write(1002,*)'D=',msdlin(1)/4.0d0,msdlinzero(1)/4.0d0
    write(1002,*)'acceptance ratio kx=',kx
  enddo

  open(202,file='mc_data')
  therm=0.0d0
  dtherm=0.0d0
  do i=1,6
    do j=1,nseries
      therm(i)=therm(i)+thermo(i,j)
    enddo
    therm(i)=therm(i)/dble(nseries)
    do j=1,nseries
      dtherm(i)=dtherm(i)+(thermo(i,j)-therm(i))**2
    enddo
    dtherm(i)=dsqrt(dtherm(i)/dble(nseries))
    write(202,'(2e16.7)')therm(i),dtherm(i)
  enddo
  close(202)

  open(202,file='!mc_data')
  write(202,'(8e16.7)')temp,dens,therm(1),therm(2),therm(3),therm(4),therm(5),therm(6)
  close(202)

  open(202,file='mc_corr')
  gr=0.0d0
  dgr=0.0d0
  do i=1,10000
    do j=1,nseries
      gr(i)=gr(i)+ggr(i,j)
    enddo
    gr(i)=gr(i)/dble(nseries)
    do j=1,nseries
      dgr(i)=dgr(i)+(ggr(i,j)-gr(i))**2
    enddo
    dgr(i)=dsqrt(dgr(i)/dble(nseries))
    write(202,'(3f16.7)')(dble(i)-0.5d0)*grinterv,gr(i),dgr(i)
  enddo
  close(202)

  open(202,file='mc_msd')
  msd=0d0
  dmsd=0d0
  do i=1,nmsdpoints
    do j=1,nseries
      msd(i)=msd(i)+mmsd(i,j)
    enddo
    msd(i)=msd(i)/dble(nseries)
    do j=1,nseries
      dmsd(i)=dmsd(i)+(mmsd(i,j)-msd(i))**2
    enddo
    dmsd(i)=dsqrt(dmsd(i)/dble(nseries))
    write(202,'(3f16.7)')dble(i*nmsd)*1d0,msd(i),dmsd(i)
  enddo
  close(202)

  call lin_regression_msd(1d0,1,nmsdpoints)
  open(202,file='!mc_msdlin_stat')
  write(202,'(10e16.7)')msdlin(1),msdlin(2),msdlin(3),msdlin(4),msdlin(5),msdlin(6),msdlin(7),msdlin(8),msdlin(9),msdlin(10)
  close(202)
  open(202,file='mc_msdlin')
  write(202,'(2f16.7)')0d0,msdlin(2)
  write(202,'(2f16.7)')dble(nmsdpoints*nmsd)*1d0,msdlin(1)*dble(nmsdpoints*nmsd)*1d0+msdlin(2)
  close(202)
  open(202,file='!mc_msdlinzero_stat')
  write(202,'(8e16.7)')msdlinzero(1),msdlinzero(2),msdlinzero(3),msdlinzero(4),msdlinzero(5),msdlinzero(6),msdlinzero(7),&
  msdlinzero(8)
  close(202)
  open(202,file='mc_msdlinzero')
  write(202,'(2f16.7)')0d0,0d0
  write(202,'(2f16.7)')dble(nmsdpoints*nmsd)*1d0,msdlinzero(1)*dble(nmsdpoints*nmsd)*1d0
  close(202)

  deallocate( msd )
  deallocate( mmsd )
  deallocate( dmsd )
  deallocate( intpos )
  deallocate( difpos )
  !------------------------------------------------------------------
  !MD
  !------------------------------------------------------------------
  nmsdpoints=nsteps/nmsd
  allocate( msd(nmsdpoints) )
  allocate( mmsd(nmsdpoints,20) )
  allocate( dmsd(nmsdpoints) )
  allocate( intpos(N,2) )
  allocate( difpos(N,2) )

  call init_pos()
  call snapshot('int')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'positions initialized'

  call init_vel()
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'velocities initialized'

  call mdequil()
  call snapshot('ekd')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MD equilibration successful'
  write(1002,*)'uu/N=',uu/dble(N)
  write(1002,*)'pp=',pp

  write(223,*)'MD'
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MD simulation initialized'
  do j=1,nseries
    call mdseries(j)
    write(fmt,'(a1,i2.2)')'d',j
    call snapshot(trim(fmt))
    cv=(uu2-uu**2)/temp**2
    thermo(1,j)=uu/dble(N)
    thermo(2,j)=pp
    thermo(3,j)=cv/dble(N)
    thermo(4,j)=top
    thermo(5,j)=msdlin(1)/4.0d0
    thermo(6,j)=msdlinzero(1)/4.0d0
    ggr(:,j)=gr(:)
    mmsd(:,j)=msd(:)
    write(223,'(8e16.7)')temp,dens,uu/dble(N),pp,cv/dble(N),top,msdlin(1)/4.0d0,msdlinzero(1)/4.0d0
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*) j,'-th MD run successful'
    write(1002,*)'uu/N=',uu/dble(N)
    write(1002,*)'pp=',pp
    write(1002,*)'cv/N=',cv/dble(N)
    write(1002,*)'top=',top
    write(1002,*)'D=',msdlin(1)/4.0d0,msdlinzero(1)/4.0d0
  enddo
  close(223)

  open(202,file='md_data')
  therm=0.0d0
  dtherm=0.0d0
  do i=1,6
    do j=1,nseries
      therm(i)=therm(i)+thermo(i,j)
    enddo
    therm(i)=therm(i)/dble(nseries)
    do j=1,nseries
      dtherm(i)=dtherm(i)+(thermo(i,j)-therm(i))**2
    enddo
    dtherm(i)=dsqrt(dtherm(i)/dble(nseries))
    write(202,'(2e16.7)')therm(i),dtherm(i)
  enddo
  close(202)

  open(202,file='!md_data')
  write(202,'(8e16.7)')temp,dens,therm(1),therm(2),therm(3),therm(4),therm(5),therm(6)
  close(202)

  open(202,file='md_corr')
  gr=0.0d0
  dgr=0.0d0
  do i=1,10000
    do j=1,nseries
      gr(i)=gr(i)+ggr(i,j)
    enddo
    gr(i)=gr(i)/dble(nseries)
    do j=1,nseries
      dgr(i)=dgr(i)+(ggr(i,j)-gr(i))**2
    enddo
    dgr(i)=dsqrt(dgr(i)/dble(nseries))
    write(202,'(3f16.7)')(dble(i)-0.5d0)*grinterv,gr(i),dgr(i)
  enddo
  close(202)

  open(202,file='md_msd')
  msd=0d0
  dmsd=0d0
  do i=1,nmsdpoints
    do j=1,nseries
      msd(i)=msd(i)+mmsd(i,j)
    enddo
    msd(i)=msd(i)/dble(nseries)
    do j=1,nseries
      dmsd(i)=dmsd(i)+(mmsd(i,j)-msd(i))**2
    enddo
    dmsd(i)=dsqrt(dmsd(i)/dble(nseries))
    write(202,'(3f16.7)')dble(i*nmsd)*tstep,msd(i),dmsd(i)
  enddo
  close(202)

  call lin_regression_msd(tstep,1,nmsdpoints)
  open(202,file='!md_msdlin_stat')
  write(202,'(10e16.7)')msdlin(1),msdlin(2),msdlin(3),msdlin(4),msdlin(5),msdlin(6),msdlin(7),msdlin(8),msdlin(9),msdlin(10)
  close(202)
  open(202,file='md_msdlin')
  write(202,'(2f16.7)')0d0,msdlin(2)
  write(202,'(2f16.7)')dble(nmsdpoints*nmsd)*tstep,msdlin(1)*dble(nmsdpoints*nmsd)*tstep+msdlin(2)
  close(202)
  open(202,file='!md_msdlinzero_stat')
  write(202,'(8e16.7)')msdlinzero(1),msdlinzero(2),msdlinzero(3),msdlinzero(4),msdlinzero(5),msdlinzero(6),msdlinzero(7),&
  msdlinzero(8)
  close(202)
  open(202,file='md_msdlinzero')
  write(202,'(2f16.7)')0d0,0d0
  write(202,'(2f16.7)')dble(nmsdpoints*nmsd)*tstep,msdlinzero(1)*dble(nmsdpoints*nmsd)*tstep
  close(202)

  call cpu_time(time2)
  call date_and_time(VALUES=values2)
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'CPU simulation time=',time2-time1
  write(1002,*)'start and finish date and time YYYY, MM, DD, UTC, HH, MIN, SS'
  write(1002,'(8i5)')values1
  write(1002,'(8i5)')values2
  time2=(values2(8)-values1(8))/1000.0d0+values2(7)-values1(7)
  time2=time2+(values2(6)-values1(6))*60.0d0
  time2=time2+(values2(5)-values1(5))*60.0d0*60.0d0
  time2=time2+(values2(3)-values1(3))*60.0d0*60.0d0*24.0d0
  write(1002,*)'real simulation time=',time2
  close(1002)

  deallocate( pos )
  deallocate( vel )
  deallocate( acc )
  deallocate( pot )
  deallocate( tla )
  deallocate( dis )
  deallocate( rrr )
  deallocate( msd )
  deallocate( mmsd )
  deallocate( dmsd )
  deallocate( intpos )
  deallocate( difpos )
endprogram core

!----------------------------------------------------------------------------------------------------
!Generates random initial positions for particles inside the box without overlapping
!----------------------------------------------------------------------------------------------------
subroutine random_pos()
  use parameters
  use particles
  implicit none
  real*8 :: ran3,image !functions
  real*8 :: a
  real*8 :: x,y,r
  integer :: i,j,k

  !first particle
  a=ran3(ran)
  pos(1,1)=(a-0.5d0)*lbox
  a=ran3(ran)
  pos(1,2)=(a-0.5d0)*lbox
  
  !prevent overlapping
  do i=2,N
    j=0
    do while (j.lt.0.5)
      a=ran3(ran)
      pos(i,1)=(a-0.5d0)*lbox
      a=ran3(ran)
      pos(i,2)=(a-0.5d0)*lbox
      j=1
      do k=1,i-1
        x=image(pos(i,1),pos(k,1),lbox)
        y=image(pos(i,2),pos(k,2),lbox)
        r=dsqrt(x**2+y**2)
        if (r.lt.0.5d0) j=0
      enddo
    enddo
  enddo
endsubroutine random_pos

!----------------------------------------------------------------------------------------------------
!Ran3 random number generator generates real numbers [0,1)
!----------------------------------------------------------------------------------------------------
function ran3(idum)
  implicit none
  integer :: idum
  integer :: MBIG,MSEED,MZ
  real*8 :: ran3,FAC
  PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.0d0/dble(MBIG))
  integer :: i,iff,ii,inext,inextp,k
  integer :: mj,mk,ma(55)
  SAVE iff,inext,inextp,ma
  DATA iff /0/
  
  if (idum.lt.0.or.iff.eq.0) then
    iff=1
    mj=MSEED-iabs(idum)
    mj=mod(mj,MBIG)
    ma(55)=mj
    mk=1
    do i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if (mk.lt.MZ) mk=mk+MBIG
      mj=ma(ii)
    enddo
    do k=1,4
      do i=1,55
        ma(i)=ma(i)-ma(1+mod(i+30,55))
        if (ma(i).lt.MZ) ma(i)=ma(i)+MBIG
      enddo
    enddo
    inext=0
    inextp=31
    idum=1
  endif
  inext=inext+1
  if (inext.eq.56) inext=1
  inextp=inextp+1
  if (inextp.eq.56) inextp=1
  mj=ma(inext)-ma(inextp)
  if (mj.lt.MZ) mj=mj+MBIG
  ma(inext)=mj
  ran3=mj*FAC
  return
endfunction ran3

!----------------------------------------------------------------------------------------------------
!Minimum image convention for periodic boundary conditions
!----------------------------------------------------------------------------------------------------
function image(xa,xb,ll)
  implicit none
  real*8 :: xa,xb,ll
  real*8 :: xc
  real*8 :: image

  xc=xb-xa
  xc=xc-ll*dnint(xc/ll)      !rounds double to nearest double
  image=xc
  return
endfunction image

!----------------------------------------------------------------------------------------------------
!Removes bad contacts
!----------------------------------------------------------------------------------------------------
subroutine minimize()
  use parameters
  use mc
  implicit none
  integer :: i

  accx=0

  open(202,file='mc_min')
  call interactions()
  do i=1,10000*N
    call mcmove1(0)
    if (mod(i,nsampl).eq.0) then
      write(202,'(3f16.7)')dble(i),u,p
    endif
  enddo
  close(202)

  kx=dble(accx)/dble(10000*N)
endsubroutine minimize

!----------------------------------------------------------------------------------------------------
!Calculates interactions between all particles
!----------------------------------------------------------------------------------------------------
subroutine interactions()
  use parameters
  use particles
  use potential
  use mc
  implicit none
  real*8 :: image,cspot,tl !functions
  real*8 :: x,y,r
  integer :: i,j

  !potential module
  pot=0.0d0
  tla=0.0d0
  dis=0.0d0
  rrr=0.0d0

  !mc module
  u=0.0d0
  p=0.0d0

  !interactions
  do i=1,N-1
    do j=i+1,N
      x=image(pos(i,1),pos(j,1),lbox)
      y=image(pos(i,2),pos(j,2),lbox)
      r=dsqrt(x**2+y**2)
      if (r.lt.rcut) then
        !potential module
        pot(i,j)=cspot(r,a1,a2,a3,a4,a5,tl)
        tla(i,j)=tl
        dis(i,j)=r
        rrr(i,j,1)=x
        rrr(i,j,2)=y
        pot(j,i)=pot(i,j)
        tla(j,i)=tl
        dis(j,i)=r
        rrr(j,i,1)=-x
        rrr(j,i,2)=-y
        !mc module
        u=u+pot(i,j)
        p=p+tla(i,j)
      else
        !potential module
        dis(i,j)=r
        rrr(i,j,1)=x
        rrr(i,j,2)=y
        dis(j,i)=r
        rrr(j,i,1)=-x
        rrr(j,i,2)=-y
      endif
    enddo
  enddo
endsubroutine interactions

!----------------------------------------------------------------------------------------------------
!Calculates the CS potential between two particles
!----------------------------------------------------------------------------------------------------
function cspot(r,a1,a2,a3,a4,a5,p)
  implicit none
  !using reduced units
  real*8 :: r
  real*8 :: a1,a2
  !a1=epsilon = absolute value of the minimum value of the LJ potential (depth of the potential well)
  !a2=sigma = distance at which the potential becomes positive
  real*8 :: a3,a4,a5
  !a3=u0 = height of the gaussian ramp maximum
  !a4=c0 = width of the gaussian ramp
  !a4=(1/(a2*a4))**2=(1/(sigma*c0))**2 = width factor in the gaussian function
  !a5=r0 = distance at which the gaussian ramp reaches its maximum
  real*8 :: x,gauss
  real*8 :: cspot,p !CS potential, virial pressure

  x=a2/r
  x=x**6
  gauss=a3*dexp(-a4*(r-a5)**2)
  cspot=4.0d0*a1*x*(x-1.0d0)+gauss
  p=48.0d0*a1*x*(x-0.5d0)+2.0d0*a4*(r-a5)*gauss*r          ! -dU/dr * r = F * r
  return
endfunction cspot

!----------------------------------------------------------------------------------------------------
!Trial move for 1 particle
!----------------------------------------------------------------------------------------------------
subroutine mcmove1(ii)
  use parameters
  use particles
  use potential
  use mc
  use diffusion
  implicit none
  integer :: ii !switch
  real*8 :: ran3,image,cspot,tl !functions
  real*8 :: a
  real*8 :: pos1(2),pot1(N),tla1(N),dis1(N)
  real*8 :: dx,dy,x,y,r
  real*8 :: du,dp
  integer :: j,k

  !random particle
  a=ran3(ran)
  k=int(a*N)+1        !floors double to integer
  
  !new position
  a=ran3(ran)
  dx=(a-0.5d0)*delta*a2
  pos1(1)=pos(k,1)+dx
  pos1(1)=image(0.0d0,pos1(1),lbox)
  a=ran3(ran)
  dy=(a-0.5d0)*delta*a2
  pos1(2)=pos(k,2)+dy
  pos1(2)=image(0.0d0,pos1(2),lbox)

  !new interactions
  do j=1,N
    if (j.eq.k) then
      pot1(j)=0.0d0
      tla1(j)=0.0d0
      dis1(j)=0.0d0
    else
      x=image(pos1(1),pos(j,1),lbox)
      y=image(pos1(2),pos(j,2),lbox)
      r=dsqrt(x**2+y**2)
      if (r.lt.rcut) then
        pot1(j)=cspot(r,a1,a2,a3,a4,a5,tl)
        tla1(j)=tl
        dis1(j)=r
      else
        pot1(j)=0.0d0
        tla1(j)=0.0d0
        dis1(j)=r
      endif
    endif
  enddo

  !change in potential energy
  du=0.0d0
  dp=0.0d0
  do j=1,N
    du=du+pot1(j)-pot(k,j)
    dp=dp+tla1(j)-tla(k,j)
  enddo

  !accept or reject
  a=ran3(ran)
  if(a.lt.dexp(-du/temp)) then
    pos(k,:)=pos1(:)
    if (ii.eq.1) then
      difpos(k,1)=difpos(k,1)+dx
      difpos(k,2)=difpos(k,2)+dy
    endif
    do j=1,N
      !potential module
      pot(k,j)=pot1(j)
      tla(k,j)=tla1(j)
      dis(k,j)=dis1(j)
      pot(j,k)=pot1(j)
      tla(j,k)=tla1(j)
      dis(j,k)=dis1(j)
    enddo
    !mc module
    u=u+du
    p=p+dp
    accx=accx+1
  endif
endsubroutine mcmove1

!----------------------------------------------------------------------------------------------------
!MC equilibration
!----------------------------------------------------------------------------------------------------
subroutine mcequil()
  use parameters
  use mc
  implicit none
  integer :: i

  uu=0.0d0
  pp=0.0d0
  m=0
  accx=0

  open(202,file='mc_ekv')
  call interactions()
  do i=1,nekv*N
    call mcmove1(0)
    if (mod(i,nsampl).eq.0) then
      write(202,'(3f16.7)')dble(i),u,p
      uu=uu+u
      pp=pp+p
      m=m+1
    endif
  enddo
  close(202)

  uu=uu/dble(m)
  pp=pp/dble(m)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
  kx=dble(accx)/dble(nekv*N)
endsubroutine mcequil

!----------------------------------------------------------------------------------------------------
!MC sampling
!----------------------------------------------------------------------------------------------------
subroutine mcseries(idum)
  use parameters
  use particles
  use mc
  use correlation
  use diffusion
  implicit none
  integer :: idum
  character(8) :: fmt
  real*8 :: c                                   !numerical constant
  integer :: i

  uu=0.0d0
  uu2=0.0d0
  pp=0.0d0
  m=0
  accx=0

  gr=0.0d0
  mmm=0
  top=0.0d0

  intpos=pos
  difpos=pos

  write(fmt,'(i3.3)')idum
  open(202,file='mc_sam'//trim(fmt))
  call interactions()
  do i=1,ncycl*N
    call mcmove1(1)
    if (mod(i,nsampl).eq.0) then
      write(202,'(3f16.7)')dble(i),u,p
      uu=uu+u
      uu2=uu2+u**2
      pp=pp+p
      m=m+1
      if (mod(i,ngr*N).eq.0) then
        call corr()
        mmm=mmm+1
      endif
      if (mod(m,nmsd).eq.0) then
        call mean_squared_displacement(m/nmsd)
      endif
    endif
  enddo
  close(202)

  uu=uu/dble(m)
  uu2=uu2/dble(m)
  pp=pp/dble(m)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
  kx=dble(accx)/dble(ncycl*N)
  
  do i=1,10000
    c=Pi*dble(2*i-1)*grinterv**2
    gr(i)=gr(i)/dble(N-1)/dble(mmm)/dens/c
  enddo
  write(fmt,'(i3.3)')idum
  open(10,file='mc_corr'//trim(fmt))
  do i=1,10000
    write(10,'(2f16.7)')(dble(i)-0.5d0)*grinterv,gr(i)
  enddo
  close(10)

  call trans_order_parameter()

  write(fmt,'(i3.3)')idum
  open(10,file='mc_msd'//trim(fmt))
  do i=1,nmsdpoints
    write(10,'(2f16.7)')dble(i*nmsd)*1d0,msd(i)
  enddo
  close(10)

  call lin_regression_msd(1.0d0,1,nmsdpoints)
  write(fmt,'(i3.3)')idum
  open(10,file='mc_msdlin_stat'//trim(fmt))
  write(10,'(10e16.7)')msdlin(1),msdlin(2),msdlin(3),msdlin(4),msdlin(5),msdlin(6),msdlin(7),msdlin(8),msdlin(9),msdlin(10)
  close(10)
  open(10,file='mc_msdlin'//trim(fmt))
  write(10,'(2f16.7)')0d0,msdlin(2)
  write(10,'(2f16.7)')dble(nmsdpoints*nmsd)*1d0,msdlin(1)*dble(nmsdpoints*nmsd)*1d0+msdlin(2)
  close(10)
  open(10,file='mc_msdlinzero_stat'//trim(fmt))
  write(10,'(8e16.7)')msdlinzero(1),msdlinzero(2),msdlinzero(3),msdlinzero(4),msdlinzero(5),msdlinzero(6),msdlinzero(7),&
  msdlinzero(8)
  close(10)
  open(10,file='mc_msdlinzero'//trim(fmt))
  write(10,'(2f16.7)')0d0,0d0
  write(10,'(2f16.7)')dble(nmsdpoints*nmsd)*1d0,msdlinzero(1)*dble(nmsdpoints*nmsd)*1d0
  close(10)
endsubroutine mcseries

!----------------------------------------------------------------------------------------------------
!Pair correlation function g(r) -- Radial distribution function
!----------------------------------------------------------------------------------------------------
subroutine corr()
  use parameters
  use potential
  use correlation
  implicit none
  integer :: i,j
  integer :: k
  real*8 :: r

  do i=1,N-1
    do j=i+1,N
      r=dis(i,j)
      k=int(r/grinterv)        !floors double to integer
      if (k.lt.10000) then
        gr(k+1)=gr(k+1)+2.0d0
      endif
    enddo
  enddo
endsubroutine corr

!--------------------------------------------------------
!translational order parameter
!--------------------------------------------------------
subroutine trans_order_parameter()
  use parameters
  use correlation
  implicit none
  integer :: i
  integer :: npoints
  real*8 :: r

  npoints=0
  r=0d0
  top=0d0

  do while(r.le.(lbox/2.0d0))
    npoints=npoints+1
    r=(dble(npoints)-0.5d0)*grinterv
  enddo

  !Trapeze method of integration
  do i=1,npoints-1
    top=top+(abs(gr(i)-1.0d0)+abs(gr(i+1)-1.0d0))*grinterv/2.0d0
  enddo

  top=sqrt(dens)*top
endsubroutine trans_order_parameter

!---------------------------------------------------------------
!mean squared displacment msd(n) or msd(t)
!---------------------------------------------------------------
subroutine mean_squared_displacement(point)
  use parameters
  use diffusion
  implicit none
  integer :: point
  integer :: i
  real*8 :: sd

  sd=0d0
  do i=1,N
    sd=sd+(difpos(i,1)-intpos(i,1))**2+(difpos(i,2)-intpos(i,2))**2
  enddo

  !average
  sd=sd/dble(N)

  msd(point)=sd
endsubroutine mean_squared_displacement

!--------------------------------------------------------
!linear regression of msd
!--------------------------------------------------------
subroutine lin_regression_msd(step,initial,final)
  use diffusion
  implicit none
  real*8 :: step
  integer :: i,initial,final,npoints,df
  real*8 :: x,y
  real*8 :: vx,vy,vx2,vxy,vy2
  real*8 :: b,a,sb,sa,r2,sy,f,ssr,sse

  vx=0d0
  vy=0d0
  vx2=0d0
  vxy=0d0
  vy2=0d0

  do i=initial,final
    x=dble(i*nmsd)*step !CYCLE STEP or TIME STEP
    y=msd(i)
    vx=vx+x
    vy=vy+y
    vx2=vx2+x**2
    vxy=vxy+x*y
    vy2=vy2+y**2
  enddo
  
  !averages
  npoints=final-initial+1
  vx=vx/dble(npoints)
  vy=vy/dble(npoints)
  vx2=vx2/dble(npoints)
  vxy=vxy/dble(npoints)
  vy2=vy2/dble(npoints)

  !linear regression with zero intercept term
  b=(vxy-vx*vy)/(vx2-vx**2)
  a=vy-b*vx
  ssr=dble(npoints)*b*(vxy-vx*vy)
  sse=dble(npoints)*(vy2-vy**2)-ssr
  df=npoints-2
  f=ssr/sse*dble(df)
  sy=dsqrt(sse/df)
  r2=ssr/(ssr+sse)
  sb=sy/dsqrt(npoints*(vx2-vx**2))
  sa=sy*dsqrt(1/npoints+vx**2/(npoints*(vx2-vx**2)))

  msdlin(1)=b
  msdlin(2)=a
  msdlin(3)=sb
  msdlin(4)=sa
  msdlin(5)=r2
  msdlin(6)=sy
  msdlin(7)=f
  msdlin(8)=df
  msdlin(9)=ssr
  msdlin(10)=sse

  !linear regression without intercept term
  b=vxy/vx2
  ssr=dble(npoints)*b*vxy
  sse=dble(npoints)*vy2-ssr
  df=npoints-1
  f=ssr/sse*dble(df)
  sy=dsqrt(sse/df)
  r2=ssr/(ssr+sse)
  sb=sy/dsqrt(npoints*vx2)

  msdlinzero(1)=b
  msdlinzero(2)=sb
  msdlinzero(3)=r2
  msdlinzero(4)=sy
  msdlinzero(5)=f
  msdlinzero(6)=df
  msdlinzero(7)=ssr
  msdlinzero(8)=sse
endsubroutine lin_regression_msd

!---------------------------------------------------------------------------------------------
!Initializes positions for a homogeneous fluid
!---------------------------------------------------------------------------------------------
subroutine init_pos()
  use parameters
  use particles
  implicit none
  integer :: i
  real*8 :: image !function
  real*8 :: sumpos(2)

  !setting the sum of positions to zero
  sumpos=1.0d0
  do while (sum(sumpos*sumpos).gt.1.0d-5)
    sumpos(1)=sum(pos(:,1))/dble(N)
    sumpos(2)=sum(pos(:,2))/dble(N)
    do i=1,N
      pos(i,1)=image(sumpos(1),pos(i,1),lbox)
      pos(i,2)=image(sumpos(2),pos(i,2),lbox)
    enddo
  enddo
endsubroutine init_pos

!-------------------------------------------------------------------------------------------------------------------------------
!Initializes velocities for an isotropic fluid
!Box-Muller transform for Maxwell-Boltzmann distribution
!-------------------------------------------------------------------------------------------------------------------------------
subroutine init_vel()
  use parameters
  use particles
  implicit none
  integer :: i
  real*8 :: ran3 !function
  real*8 :: sigma,x1,x2,s,z1,z2
  real*8 :: sumvel(2)
  real*8 :: rescale

  !Box-Muller transform
  sigma=dsqrt(temp)
  do i=1,N
    s=2.0d0
    do while ((s.eq.0.0d0).or.(s.ge.1.0d0))
      x1=2.0d0*ran3(ran)-1.0d0
      x2=2.0d0*ran3(ran)-1.0d0
      s=x1**2+x2**2
    enddo
    z1=x1*dsqrt(-2.0d0*dlog(s)/s)
    z2=x2*dsqrt(-2.0d0*dlog(s)/s)
    vel(i,1)=sigma*z1
    vel(i,2)=sigma*z2
  enddo

  !setting the sum of velocities to zero
  sumvel=1.0d0
  do while (sum(sumvel*sumvel).gt.1.0d-5)
    sumvel(1)=sum(vel(:,1))/dble(N)
    sumvel(2)=sum(vel(:,2))/dble(N)
    vel(:,1)=vel(:,1)-sumvel(1)
    vel(:,2)=vel(:,2)-sumvel(2)
  enddo

  !calculate the temperature
  call temperature()

  !rescale the particle velocities to the desired temperature
  rescale=dsqrt(temp/tem)
  vel=vel*rescale
endsubroutine init_vel

!----------------------------------------------------------------------------------------------------
!Calculates current temperature using the equipartition theorem
!----------------------------------------------------------------------------------------------------
subroutine temperature()
  use parameters
  use particles
  implicit none
  integer :: i

  ekin=0.0d0  
  do i=1,N
    ekin=ekin+(vel(i,1)**2+vel(i,2)**2)/2.0d0
  enddo
  !(x,y) -> 2 degrees of freedom
  !1/2 T for each degree of freedom
  !2N degrees of freedom (N particles with (x,y))
  tem=ekin/dble(N)
endsubroutine temperature

!----------------------------------------------------------------------------------------------------
!MD equilibration
!----------------------------------------------------------------------------------------------------
subroutine mdequil()
  use parameters
  use particles
  use mc
  implicit none
  integer :: i
  real*8 :: t

  t=0.0d0
  uu=0.0d0
  pp=0.0d0
  m=0

  open(202,file='md_ekv')
  call forces()
  call temperature()
  do i=1,nekv
    call berendsen()
    t=t+tstep
    etot=ekin+u
    write(202,'(6f16.7)')t,u,ekin,etot,p,tem
  enddo
  call forces()
  call temperature()
  do i=1,nekv
    call nose_hoover_chain(0)
    t=t+tstep
    etot=ekin+u
    write(202,'(6f16.7)')t,u,ekin,etot,p,tem
    uu=uu+u
    pp=pp+p
    m=m+1
  enddo
  close(202)

  uu=uu/dble(m)
  pp=pp/dble(m)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
endsubroutine mdequil

!----------------------------------------------------------------------
!Calculates forces from interactions
!----------------------------------------------------------------------
subroutine forces()
  use parameters
  use particles
  use potential
  implicit none
  integer :: i,j
  real*8 :: f(N,N,2)

  call interactions()
  f=0.0d0
  do i=1,N-1
    do j=i+1,N
      f(i,j,1)=tla(i,j)*rrr(i,j,1)/dis(i,j)**2  
      f(i,j,2)=tla(i,j)*rrr(i,j,2)/dis(i,j)**2
      f(j,i,1)=-f(i,j,1)
      f(j,i,2)=-f(i,j,2)
    enddo
  enddo
  do j=1,N
    acc(j,1)=sum(f(:,j,1))
    acc(j,2)=sum(f(:,j,2))
  enddo
endsubroutine forces

!-------------------------------------------------------------------------
!Velocity Verlet algorithm for calculating velocities
!Berendsen rescaling thermostat with temperature time constant tt
!Newton's equations of motion
!-------------------------------------------------------------------------
subroutine berendsen()
  use parameters
  use particles
  implicit none
  integer :: i
  real*8 :: image !function
  real*8 :: rescale

  do i=1,N
    !half-update the particle velocities by tstep/2
    vel(i,1)=vel(i,1)+acc(i,1)*dt2
    vel(i,2)=vel(i,2)+acc(i,2)*dt2
    !update the particle positions by tstep
    pos(i,1)=pos(i,1)+vel(i,1)*tstep
    pos(i,2)=pos(i,2)+vel(i,2)*tstep
    pos(i,1)=image(0.0d0,pos(i,1),lbox)
    pos(i,2)=image(0.0d0,pos(i,2),lbox)
  enddo

  !update the forces
  call forces()

  do i=1,N
    !update the particle velocities by tstep/2
    vel(i,1)=vel(i,1)+acc(i,1)*dt2
    vel(i,2)=vel(i,2)+acc(i,2)*dt2
  enddo
  
  !update the temperature
  call temperature()

  !rescale the particle velocities to the desired temperature
  rescale=dsqrt(1.0d0+tstep/tt*(temp/tem-1.0d0))
  vel=vel*rescale
endsubroutine berendsen

!-------------------------------------------------------------------------------------------------------
!Velocity Verlet algorithm for calculating velocities
!Nose-Hoover chain thermostat with 10 heat bath coupling massess Q
!Extended Lagrangian mechanics -- Martyna-Tuckerman-Tobias-Klein equations of motion
!-------------------------------------------------------------------------------------------------------
subroutine nose_hoover_chain(ii)
  use parameters
  use particles
  use diffusion
  implicit none
  integer :: ii
  integer :: i
  real*8 :: image !function
  real*8 :: dx,dy

  !half-update the thermostat
  call nhc_thermostat()
  !half-update the particle velocities
  do i=1,N
    vel(i,1)=vel(i,1)+acc(i,1)*dt2
    vel(i,2)=vel(i,2)+acc(i,2)*dt2
  enddo
  !update the particle positions
  do i=1,N
    dx=vel(i,1)*tstep
    dy=vel(i,2)*tstep
    pos(i,1)=pos(i,1)+dx
    pos(i,2)=pos(i,2)+dy
    if (ii.eq.1) then
      difpos(i,1)=difpos(i,1)+dx
      difpos(i,2)=difpos(i,2)+dy
    endif
    pos(i,1)=image(0.0d0,pos(i,1),lbox)
    pos(i,2)=image(0.0d0,pos(i,2),lbox)
  enddo
  !update the forces
  call forces()
  !update the particle velocities
  do i=1,N
    vel(i,1)=vel(i,1)+acc(i,1)*dt2
    vel(i,2)=vel(i,2)+acc(i,2)*dt2
  enddo
  !update the thermostat
  call temperature()
  call nhc_thermostat()
endsubroutine nose_hoover_chain

!-------------------------------------------------------------------------------------------------------
!Nose-Hoover chain thermostat with 10 heat bath coupling massess Q
!-------------------------------------------------------------------------------------------------------
subroutine nhc_thermostat()
  use parameters
  use particles
  implicit none
  integer :: i
  real*8 :: X
  real*8 :: AA

  G(1)=(2.0d0*ekin-NfkT)/Q(1)
  vksi(10)=vksi(10)+G(10)*dt4
  do i=1,9
    X=dexp(-vksi(11-i)*dt8)
    vksi(10-i)=vksi(10-i)*X
    vksi(10-i)=vksi(10-i)+G(10-i)*dt4
    vksi(10-i)=vksi(10-i)*X
  enddo

  AA=dexp(-vksi(1)*dt2)
  !scale the particle velocities
  do i=1,N
    vel(i,1)=vel(i,1)*AA
    vel(i,2)=vel(i,2)*AA
  enddo
  !update the temperature
  ekin=ekin*AA**2
  tem=tem*AA**2
  
  G(1)=(2.0d0*ekin-NfkT)/Q(1)
  do i=1,9
    X=dexp(-vksi(i+1)*dt8)
    vksi(i)=vksi(i)*X
    vksi(i)=vksi(i)+G(i)*dt4
    vksi(i)=vksi(i)*X
    G(i+1)=(Q(i)*vksi(i)**2-temp)/Q(i+1)
  enddo
  vksi(10)=vksi(10)+G(10)*dt4
endsubroutine nhc_thermostat

!----------------------------------------------------------------------------------------------------
!MD sampling
!----------------------------------------------------------------------------------------------------
subroutine mdseries(idum)
  use parameters
  use particles
  use mc
  use correlation
  use diffusion
  implicit none
  integer :: idum
  character(8) :: fmt
  real*8 :: c                                   !numerical constant
  integer :: i
  real*8 :: t

  t=0.0d0
  uu=0.0d0
  uu2=0.0d0
  pp=0.0d0
  m=0

  gr=0.0d0
  mmm=0
  top=0.0d0

  intpos=pos
  difpos=pos

  write(fmt,'(i3.3)')idum
  open(202,file='md_sam'//trim(fmt))
  call forces()
  call temperature()
  do i=1,nsteps
    call nose_hoover_chain(1)
    t=t+tstep
    etot=ekin+u
    write(202,'(6f16.7)')t,u,ekin,etot,p,tem
    uu=uu+u
    uu2=uu2+u**2
    pp=pp+p
    m=m+1
    if (mod(i,ngr).eq.0) then
      call corr()
      mmm=mmm+1
    endif
    if ( (idum.eq.1).and.(mod(i,nframe).eq.0).and.(i/nframe.lt.501) ) then
      write(fmt,'(i2.2,i6.6)')idum,i/nframe
      call frame(fmt)
    endif
    if (mod(i,nmsd).eq.0) then
      call mean_squared_displacement(i/nmsd)
    endif
  enddo
  close(202)

  uu=uu/dble(m)
  uu2=uu2/dble(m)
  pp=pp/dble(m)
  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp

  do i=1,10000
    c=Pi*dble(2*i-1)*grinterv**2
    gr(i)=gr(i)/dble(N-1)/dble(mmm)/dens/c
  enddo
  write(fmt,'(i3.3)')idum
  open(10,file='md_corr'//trim(fmt))
  do i=1,10000
    write(10,'(2f16.7)')(dble(i)-0.5d0)*grinterv,gr(i)
  enddo
  close(10)

  call trans_order_parameter()

  write(fmt,'(i3.3)')idum
  open(10,file='md_msd'//trim(fmt))
  do i=1,nmsdpoints
    write(10,'(2f16.7)')dble(i*nmsd)*tstep,msd(i)
  enddo
  close(10)

  call lin_regression_msd(tstep,1,nmsdpoints)
  write(fmt,'(i3.3)')idum
  open(10,file='md_msdlin_stat'//trim(fmt))
  write(10,'(10e16.7)')msdlin(1),msdlin(2),msdlin(3),msdlin(4),msdlin(5),msdlin(6),msdlin(7),msdlin(8),msdlin(9),msdlin(10)
  close(10)
  open(10,file='md_msdlin'//trim(fmt))
  write(10,'(2f16.7)')0d0,msdlin(2)
  write(10,'(2f16.7)')dble(nmsdpoints*nmsd)*tstep,msdlin(1)*dble(nmsdpoints*nmsd)*tstep+msdlin(2)
  close(10)
  open(10,file='md_msdlinzero_stat'//trim(fmt))
  write(10,'(8e16.7)')msdlinzero(1),msdlinzero(2),msdlinzero(3),msdlinzero(4),msdlinzero(5),msdlinzero(6),msdlinzero(7),&
  msdlinzero(8)
  close(10)
  open(10,file='md_msdlinzero'//trim(fmt))
  write(10,'(2f16.7)')0d0,0d0
  write(10,'(2f16.7)')dble(nmsdpoints*nmsd)*tstep,msdlinzero(1)*dble(nmsdpoints*nmsd)*tstep
  close(10)
endsubroutine mdseries

!-------------------------------------------------------
!Snapshot
!-------------------------------------------------------
subroutine snapshot(fmt)
  use parameters
  use particles
  implicit none
  character(3) :: fmt
  integer :: i

  open(111,file='snapshot_'//fmt)
  do i=1,N
    write(111,'(3f16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
  enddo
  close(111)
endsubroutine snapshot

!--------------------------------------------------------
!Movie
!--------------------------------------------------------
subroutine frame(fmt)
  use parameters
  use particles
  implicit none
  character(8) :: fmt
  integer :: i

  open(111,file='frame_'//fmt)
  do i=1,N
    write(111,'(3f16.7)')pos(i,1)/lbox,pos(i,2)/lbox,a2/2.0d0/lbox
  enddo
  close(111)
endsubroutine frame