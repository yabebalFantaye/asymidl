program healpix_map_dot_product

  use healpix_types
  use alm_tools
  use pix_tools
  use extension


  IMPLICIT NONE

  !INCLUDE 'mpif.h'

  INTEGER(I4B) :: ierr,ntasks,unit,me,kk,zz,request,unnorm
  integer(i4b) ::ipix,nside
  
  REAL(DP), DIMENSION(:), ALLOCATABLE :: map
  INTEGER(I4B) :: ordering,degree
  REAL(DP), DIMENSION(0:2) :: vec,vecpix
  character(len=100) :: strn

  me=0
  unit=99
  nside=2048
  call ang2vec_ring,nside,41.75*pi/180d0,263.85*pi/180d0,vec !cmb dipole direction


  if (me==0) print*, 'nside, vec: ',nside,vec

  allocate(map(12*nside**2))

  do ipix=0,npix-1
     call pix2vec_ring(nside,ipix,vecpix)
     map(ipix) = sum(vec*vecpix)
  enddo


  open(unit,file='map_dot_product.unf',status='unknown',form='unformatted')
  write(unit) dble(map)
  close(unit)

  if (me==0) print*, ' healpix_map_dot_product done! '

end program healpix_map_dot_product
