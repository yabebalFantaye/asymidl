program get_dipdir

  use healpix_types
  use alm_tools
  use pix_tools


  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: ierr,ntasks,unit,me,kk,zz,request, unit2
  integer(i4b) ::cnt,s0,sn,src,dest,tag,stat,cnt1,tag1,cnt2,tag2

  integer(i4b) :: nspots,npix,nside, nside_spots,isnest,nlist,ispot
  integer(i4b) :: i, j
  real(dp) :: rad,vec(3)


  INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::nct_pp, nspots_pp,listpix
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE ::thph,thph0
  
  character(len=30) :: strn,str1,str2,str3
  character(len=300) :: fname



  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  unit=12+me
  unit2=11

  open(unit,file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest/nside_map.unf',status='old',form='unformatted', iostat=ierr)
  read(unit) nside
  close(unit)

  open(unit,file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest/nside_spots.unf',status='old',form='unformatted', iostat=ierr)
  read(unit) nside_spots
  close(unit)

  open(unit,file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest/rad_listpix.unf',status='old',form='unformatted', iostat=ierr)
  read(unit) rad
  close(unit)

  isnest=1
  ! open(unit,file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest/isnest.unf',status='old',form='unformatted', iostat=ierr)
  ! read(unit) isnest
  ! close(unit)


  if (me==0) print*, 'nside, nside_spot, rad, isnest: ',nside, nside_spots, rad, isnest

  write(str1,*)nside
  write(str2,*)nside_spots
  write(str3,*)int(rad)


  rad=rad*pi/180d0 !now in radian
  npix = nside2npix(nside)
  nspots = nside2npix(nside_spots)

  

  ALLOCATE(nspots_pp(0:ntasks-1))
  nspots_pp=(nspots)/ntasks
  if(mod(nspots,ntasks).ne.0) nspots_pp(0:mod(nspots,ntasks)-1)=nspots_pp(0:mod(nspots,ntasks)-1)+1
  if (me==0) print*,'sum(nspots_pp), nspots_pp = ', sum(nspots_pp), nspots_pp


  allocate(listpix(0:npix-1))

  !compute at each cpu listpix per ispot pixel
  do kk=0,nspots_pp(me)-1

     if (me.ne.0) then
        ispot=kk+sum(nspots_pp(0:me-1))
     else
        ispot=kk
     endif
        
     call pix2vec_nest(nside_spots,ispot,vec)

     call query_disc(nside,vec,rad,listpix,nlist)

     !print*, 'writing me, ispot, nlist: ',me, ispot , nlist ,listpix(0:2)

     write(strn,*)ispot
     open(unit,file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest/listpix_'//trim(adjustl(strn))//'.unf',status='unknown',form='unformatted')
     write(unit) nlist
     write(unit) listpix(0:nlist-1)
     close(unit)

     listpix=0

  enddo

  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  deallocate(listpix)


  if (me==0) then

     fname='/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest_nin'//trim(adjustl(str1))//'_nout'//trim(adjustl(str2))//'_rad'//trim(adjustl(str3))//'.unf'
     open(unit2,file=fname,status='unknown',form='unformatted')

     do ispot=0,nspots-1

        write(strn,*)ispot
        open(unit,file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest/listpix_'//trim(adjustl(strn))//'.unf',status='unknown',form='unformatted')

        read(unit) nlist

        !print*, 'reading ispot, nlist: ',ispot, nlist

        allocate(listpix(nlist))

        read(unit) listpix
        close(unit)

        write(unit2) nlist
        write(unit2) listpix(0:nlist-1)

        deallocate(listpix)
        call system('rm '//'/mn/owl1/d3/yabebalf/planck/common_files/listpix_nest/listpix_'//trim(adjustl(strn))//'.unf')
        
     enddo

     close(unit2)

  endif





  CALL MPI_FINALIZE(ierr)
  if (me==0) print*, ' get_dipdir done! '

end program get_dipdir

