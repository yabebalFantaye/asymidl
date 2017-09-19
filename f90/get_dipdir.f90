program get_dipdir

  use healpix_types
  use alm_tools
  use pix_tools


  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: ierr,ntasks,unit,me,kk,zz,request
  integer(i4b) ::cnt,s0,sn,src,dest,tag,stat,cnt1,tag1,cnt2,tag2

  integer(i4b) :: nsim,nspots,nbins,nell, npix,nside,nlbands,nbin2bin,nside_spots
  integer(i4b) :: i, j, b1, b2,sim,band_width,dl,ind(1),nbj
  real(dp):: Aps,zbounds(1:2),dmean,th,ph
  real(dp) :: diff_mean, var_sum,mres_aps

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::nct_pp, nsim_pp,ell,bins
  REAL(DP), DIMENSION(:), ALLOCATABLE ::  pp
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sres,bres,dbres, res, dres,mnres,dnres,mres,w8
  !REAL(SP), DIMENSION(:,:), ALLOCATABLE ::theta,phi,theta0,phi0
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE ::thph,thph0

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: mcres
  COMPLEX(SPC), DIMENSION(:,:,:), ALLOCATABLE :: alm
  REAL(SP), DIMENSION(:), ALLOCATABLE :: map, mm


  if (iargc() > 0) then
     print*, '=---------------------------------------------'
     print*, ' parameters should be passed in fortran unformatted file'
     print*, ' the code expects the file temp/params.unf to exist. File contains'
     print*, ' params in LONG integer (4 byte) :  dl, nlbands, nbins, nell,nside,nsim,nspots'
     print*, '  douple precision pointsource floor: Aps'
     print*, '---------------------------------------------='
     stop
  endif


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me

  !necessary integers are alphabetically ordered
  open(unit,file='temp/params.unf',status='old',form='unformatted', iostat=ierr)
  if (ierr/=0 .and. me==0) then 
     print*, '=---------------------------------------------'
     print*, ' parameters should be passed in fortran unformatted file'
     print*, ' the code expects the file temp/params.unf to exist. File contains'
     print*, ' params in LONG integer (4 byte) :  dl, nlbands, nbins, nell,nside,nsim,nspots'
     print*, ' douple precision pointsource floor: Aps'
     print*, '---------------------------------------------='
     stop
  endif
  read(unit) nbj, dl, nlbands, nbins, nell,nside,nsim,nspots,nbin2bin
  close(unit)

  if (me==0) print*, 'nbj, dl, nlbands, nbins, nell,nside,nsim,nspots',nbj,2*dl, nlbands, nbins, nell,nside,nsim,nspots

  open(unit,file='temp/Aps.unf',status='old',form='unformatted', iostat=ierr)
  read(unit) Aps
  close(unit)
  if (me==0) print*,'Aps = ',real(Aps)



  ALLOCATE(nsim_pp(0:ntasks-1))
  ALLOCATE(nct_pp(0:ntasks-1))
  nsim_pp=(nsim+1)/ntasks
  if(mod(nsim+1,ntasks).ne.0) nsim_pp(0:mod(nsim+1,ntasks)-1)=nsim_pp(0:mod(nsim+1,ntasks)-1)+1
  if (me==0) print*,'sum(nsim_pp), nsim_pp = ', sum(nsim_pp), nsim_pp
  nct_pp = 2*nlbands*nsim_pp

  npix = nside2npix(nside)
  nside_spots = npix2nside(nspots)
  ALLOCATE(w8(1:2*nside_spots,1))
  ALLOCATE(map(0:npix-1))
  allocate(alm(1,0:1,0:1))
  w8=1d0
  zbounds = (/ -1.0_dp, 1.0_dp /)


 allocate(dbres(0:nlbands-1,0:nspots-1), bres(0:nlbands-1,0:nspots-1))
 allocate(res(0:nbins-1,0:nspots-1), mres(0:nbins-1,0:nspots-1))
  allocate(dres(0:nbins-1,0:nspots-1), mnres(0:nbins-1,0:nspots-1), dnres(0:nbins-1,0:nspots-1))
  allocate(mcres(0:nbins-1, 0:nspots-1, 0:nsim-1))

  allocate(bins(0:nbins-1),ell(0:nell-1), pp(nell))


  sn = nsim_pp(me)
  allocate(thph0(0:nlbands-1,0:sn-1,1:2)) !theta and phi per processor
  allocate(thph(0:nlbands-1,0:nsim,1:2)) !the first row data and nsim mc simulations
  allocate(mm(0:nspots-1))

  if (me==0) print*, 'get_dipdir: loading data from tem/*...'

  open(unit,file='temp/mcres.unf',status='old',form='unformatted')
  read(unit) mcres
  close(unit)

  open(unit,file='temp/mres.unf',status='old',form='unformatted')
  read(unit) mres
  close(unit)

  open(unit,file='temp/mnres.unf',status='old',form='unformatted')
  read(unit) mnres
  close(unit)

  open(unit,file='temp/res.unf',status='old',form='unformatted')
  read(unit) res
  close(unit)

  open(unit,file='temp/dres.unf',status='old',form='unformatted')
  read(unit) dres
  close(unit)

  open(unit,file='temp/dnres.unf',status='old',form='unformatted')
  read(unit) dnres
  close(unit)



  open(unit,file='temp/pixwin.unf',status='old',form='unformatted')
  read(unit) pp
  close(unit)

  open(unit,file='temp/bins.unf',status='old',form='unformatted')
  read(unit) bins
  close(unit)

  open(unit,file='temp/ell.unf',status='old',form='unformatted')
  read(unit) ell
  close(unit)

!   if (me==0) then
!      print*, '------------------'
!      print*, 'res: ', res(0:10,0:1)
!      print*, '------------------'
!      print*, 'mres: ', mres(0:10,0:1)
!      print*, '------------------'
!      print*, 'dres: ', dres(0:10,0:1)
!      print*, '------------------'
!      print*, 'mnres: ', mnres(0:10,0:1)
!      print*, '------------------'
!      print*, 'dnres: ',dnres(0:10,0:1)
!      print*, '------------------'
!      stop
!   endif

  if (me==0) print*, ' loading done! '


  do kk=0,nsim_pp(me)-1

     if (me.ne.0) then
        sim=kk+sum(nsim_pp(0:me-1))
     else
        sim=kk
     endif

     if (mod(real(sim), real(nsim/10.)) == 0.0) print*, 'sim  = ',sim

     do j=0,nspots-1 
        
        do i=0,nlbands-1 

           b1=i*nbin2bin
           if (i==0) b1 = nbj

           b2=(i+1)*nbin2bin-1
           if (b2<b1) stop 'nbj must be less than nbin2bin-1'

           if (i > 2)  b2=b2+1
           if (i > 2) b1=b1+1
           if (i >= 6)  b2=b2+1
           if (i > 6)  b1=b1+1
           if (i >= 10)  b2=b2+1
           if (i > 10)  b1=b1+1
           if (i >= 14)  b2=b2+1
           if (i > 14)  b1=b1+1

           band_width = b2-b1+1


           diff_mean = 0d0
           var_sum=0d0
           do zz = b1,b2
              if (sim == 0) then  !this is the data

                 mres_aps = mres(zz,j) + Aps/pp(bins(zz))**2
                 diff_mean = diff_mean + (res(zz,j) - mres_aps)/dble(band_width)
                 var_sum = var_sum + dres(zz,j)**2 + (2d0*mres_aps*mnres(zz,j))/(2d0*(ell(bins(zz))+dl)+1d0) + dnres(zz,j)**2/2d0

              else !these are from Gaussian CMB only simulations

                 diff_mean = diff_mean + (mcres(zz,j,sim-1)-mres(zz,j))/dble(band_width) 
                 var_sum = var_sum + dres(zz,j)**2

              endif
           enddo

           !print*, 'bres = ',real(diff_mean/sqrt(var_sum))
           bres(i,j) = diff_mean/sqrt(var_sum/dble(band_width))
           if (sim==0) dbres(i,j) = bres(i,j)

          !if (i==0 .and. j==0 .and. sim==0) print*,aps, diff_mean, var_sum, bres(j,i)
        enddo
     enddo

     do i=0,nlbands-1 
        
        !get dipole  
        mm = real(bres(i,:))
        
        call map2alm(nside_spots, 1, 1, mm, alm, zbounds, w8)
        call alm2map(nside,1,1,alm,map)


        ind = maxloc(map)
        call pix2ang_ring(nside,ind(1)-1,th,ph)

        !if (sim==0) then 
           !print*, mm(0:100)
           !print*, 'ind, th, ph', ind, th, ph
        !endif

        thph0(i, kk,1)=real(th)
        thph0(i, kk,2)=real(ph)


     enddo

  enddo

  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  if (me==0) print*, 'MPI send revcieve theta & phi ..'

  ! call MPI_SENDRECV (
  ! sendbuf,sendcount,sendtype,dest,sendtag, 
  ! ...... recvbuf,recvcount,recvtype,source,recvtag, 
  ! ...... comm,status,ierr)

  if (ntasks>1) then

     dest=0  
     if (me .ne. 0) then
        tag = me
        cnt=nct_pp(me)
        call  MPI_SSEND(thph0,cnt,MPI_REAL,dest,tag,MPI_COMM_WORLD,ierr) 
     else
        do i=0,ntasks-1

           if (i==0) then
              s0 = 0
              sn = nsim_pp(0)-1
              thph(:,s0:sn,1:2)=thph0
           else

              src = i
              tag = i
              cnt=nct_pp(i)
              s0 = sum(nsim_pp(0:i-1))
              sn = sum(nsim_pp(0:i))-1
              !print*, 'i, s0,sn,cnt, tag, nlbands: ',i, s0,sn,cnt, tag, nlbands

              call MPI_RECV(thph(:,s0:sn,1:2),cnt,MPI_REAL,src,tag,MPI_COMM_WORLD,stat,ierr)
           endif
        enddo
     endif
    
  else
     thph = thph0
  endif

  if (me==0) print*, 'MPI send revcieve theta & phi done!'


  if (me==0) then
     if (me==0) print*, ' writing temp/[bres,theta,phi] .. '

     open(unit,file='temp/bres.unf',status='unknown',form='unformatted')
     write(unit) dbres
     close(unit)


     open(unit,file='temp/theta.unf',status='unknown',form='unformatted')
     write(unit) dble(thph(:,:,1))
     close(unit)

     open(unit,file='temp/phi.unf',status='unknown',form='unformatted')
     write(unit) dble(thph(:,:,2))
     close(unit)
  endif

  CALL MPI_FINALIZE(ierr)
  if (me==0) print*, ' get_dipdir done! '
end program get_dipdir

!=----------------------

function dmean(x) result (mean)

  use healpix_types
  integer(i4b) n
  real(dp)::x(:)
  real(dp)::mean
  n = size(x)
  mean = sum(x)/dble(n)
  
end function dmean

! function rmean(n,x) result mean
! integer(i4b) n
! real(sp)::x(n),mean

! mean = sum(x)/real(n)

! end function rmean


! function imean(n,x) result mean
! integer(i4b) n
! integer(i4b) ::x(n),mean

! mean = sum(x)/n

! end function imean
