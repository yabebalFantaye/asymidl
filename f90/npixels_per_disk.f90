program npixels_perdisk

  use healpix_types
  use alm_tools
  use pix_tools
  use extension


  IMPLICIT NONE

  ! INCLUDE 'mpif.h'
  INTEGER(I4B) :: ierr,ntasks,unit,me
  integer(i4b) ::cnt,s0,sn,src,dest,tag,stat,cnt1,tag1,cnt2,tag2

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::nct_pp, nsim_pp,ell,bins
  REAL(DP), DIMENSION(:), ALLOCATABLE ::  pp


  integer(i4b) :: nspots,npix,nside,nside_spots, nrad
  integer(i4b) :: i, j,nlist
  real(dp):: vec(1:3),th, ph
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rval
  INTEGER(I4B), DIMENSION(:,:), allocatable :: npixmat
  INTEGER(I4B), DIMENSION(:), allocatable :: listpix


  character(len=100) :: strn
  character(len=500) :: fname



  ! CALL MPI_INIT(ierr)
  ! CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
  ! CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  ! ALLOCATE(nsim_pp(0:ntasks-1))
  ! ALLOCATE(nct_pp(0:ntasks-1))
  ! nsim_pp=(nsim+1)/ntasks
  ! if(mod(nsim+1,ntasks).ne.0) nsim_pp(0:mod(nsim+1,ntasks)-1)=nsim_pp(0:mod(nsim+1,ntasks)-1)+1
  ! if (me==0) print*,'sum(nsim_pp), nsim_pp = ', sum(nsim_pp), nsim_pp
  ! nct_pp = 2*nlbands*nsim_pp

  unit=12+me

  nside=2048
  nrad=16
  nside_spots=16
  nspots=3072


  allocate(rval(1:nrad))
  open(unit,file='/mn/owl1/d3/yabebalf/planck/variance_radvec.unf',status='old',form='unformatted')
  read(unit) rval
  close(unit)


  allocate(npixmat(0:nspots-1,1:nrad))
  allocate(listpix(0:12*nside**2-1))

  do i=1,nrad 
     print*, 'irad, radDeg,radR = ',i,rval(i),rval(i)*pi/180d0
     do j=0,nspots-1   
        call pix2vec_ring(nside_spots,j,vec)      
        call query_disc( nside, vec, rval(i)*pi/180d0, listpix, nlist)
        npixmat(j,i)= nlist
     enddo
  enddo


  fname='/mn/owl1/d3/yabebalf/planck/common_files/n_pixels_nside2048_nspot3072_nrad16.unf'

  open(unit,file=trim(fname),status='unknown',form='unformatted')
  write(unit) npixmat
  close(unit)



  ! call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  ! CALL MPI_FINALIZE(ierr)


  if (me==0) print*, ' npixels_per_disk done! '
end program npixels_perdisk

!=----------------------

function dmean(x) result (mean)

  use healpix_types
  integer(i4b) n
  real(dp)::x(:)
  real(dp)::mean
  n = size(x)
  mean = sum(x)/dble(n)
  
end function dmean

