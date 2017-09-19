program get_dipdir

  use healpix_types
  use alm_tools
  use pix_tools
  use modif_pixtools

  use extension
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER(I4B) :: ierr,ntasks,unit,me,kk,zz,request, unit2
  integer(i4b) ::cnt,s0,sn,src,dest,tag,stat,cnt1,tag1,cnt2,tag2

  integer(i4b) :: nspots,npix,nside, nside_spots,isnest,nlist,ispot
  integer(i4b) :: i, j, nrad
  real(dp) :: rad,vec(3)
  real(dp),dimension(:), allocatable :: radvec
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE ::thph,thph0

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::nct_pp, nspots_pp
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE ::listpix,nlistvec
  
  
  character(len=30) :: strn,str1,str2,str3
  character(len=300) :: fname,rootdir,pfile,pwddir

  !check directory or file
  LOGICAL        ::fexist
  CHARACTER (255)::C_DIRSPEC

  !timing 
  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1


  CALL MPI_INIT(ierr)

  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, me, ierr)

  ! INQUIRE (DIRECTORY=".", DIRSPEC=pwddir, EXIST=fexist)
  ! if (me==0) print*,'current directory: '//trim(adjustl(pwddir))

proc_number=me

  unit=12+me
  unit2=11

  if (iargc() > 0) then
     call getarg(1,pfile)
     open(unit,file=pfile,status='old',form='unformatted', iostat=ierr)
     read(unit) nside
     read(unit) nside_spots
     read(unit) nrad
     allocate(radvec(nrad))
     read(unit) radvec
     close(unit)
  else
     nside=2048
     nside_spots=16
     nrad=16
     allocate(radvec(nrad))
     !radvec=(/50,70,90/) !
     !radvec=(/10,12,14,16,18,20,22,30,50,70,90/) !
     radvec=(/1,2,4,6,8,10,12,14,16,18,20,22,30,50,70,90/) !
     !radvec=(/70,90/) !1,2,4,6,8,10,12,14,16,18,20,22,30,50,70
  endif


  isnest=0


  !if (me==0) print*, 'nside, nside_spot, rad, isnest: ',nside, nside_spots, nrad, isnest
  !if (me==0) print*,'radvec: ',radvec

  write(str1,*)nside
  write(str2,*)nside_spots
  write(str3,*)int(nrad)

  rootdir = 'listpix_ring_nin'//trim(adjustl(str1))//'_nout'//trim(adjustl(str2))//'_nrad'//trim(adjustl(str3))//'/'
  if (me==0) call system('mkdir -p '//trim(adjustl(rootdir)))

  if (me==0) print*, 'all ispot nlist will be saved in: '//trim(rootdir)
  

  radvec=radvec*DEG2RAD !now in radian
  npix = nside2npix(nside)
  nspots = nside2npix(nside_spots)

  

  ALLOCATE(nspots_pp(0:ntasks-1))
  nspots_pp=(nspots)/ntasks
  if(mod(nspots,ntasks).ne.0) nspots_pp(0:mod(nspots,ntasks)-1)=nspots_pp(0:mod(nspots,ntasks)-1)+1
  !if (me==0) print*,'sum(nspots_pp), nspots_pp = ', sum(nspots_pp), nspots_pp


  allocate(listpix(0:npix-1),nlistvec(nrad))


  !compute at each cpu listpix per ispot pixel
  do kk=0,nspots_pp(me)-1

     if (me.ne.0) then
        ispot=kk+sum(nspots_pp(0:me-1))
     else
        ispot=kk
     endif

     write(strn,*)ispot

     inquire(file=trim(adjustl(rootdir))//'/listpix_'//trim(adjustl(strn))//'.unf',exist=fexist)
     !if (me==ntasks-1) print*, 'file, exist: ',trim(adjustl(rootdir))//'/listpix_'//trim(adjustl(strn))//'.unf',fexist


     !estimate time
     call wall_clock_time(time0)
     call cpu_time(ptime0)
     

     !if (.not. fexist) then
        !if (isnest==1) call pix2vec_nest(nside_spots,ispot,vec)

        call pix2vec_ring(nside_spots,ispot,vec)
        call query_disc_radvec(nside,vec,radvec,listpix,nlist,nlistvec)
                ! 
        open(unit,file=trim(adjustl(rootdir))//'/listpix_'//trim(adjustl(strn))//'.unf',status='unknown',form='unformatted')
        write(unit) nrad, nlist
        write(unit) nlistvec
        write(unit) listpix(0:nlist-1)
        close(unit)

     ! else
     !    open(unit,file=trim(adjustl(rootdir))//'/listpix_'//trim(adjustl(strn))//'.unf',status='old',form='unformatted')
     !    read(unit) nrad, nlist
     !    read(unit) nlistvec
     !    read(unit) listpix(0:nlist-1)
     !    close(unit)
     !endif


     call wall_clock_time(time1)
     call cpu_time(ptime1)
     clock_time = time1 - time0
     ptime      = ptime1 - ptime0
     print*,'query disk done, ispot, time=',ispot,clock_time
     

     listpix=0

  enddo

  !call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  ! deallocate(listpix)
  ! deallocate(nlistvec)
  ! deallocate(radvec)

  ! if (me==0) then

  !    fname='/mn/owl1/d3/yabebalf/planck/common_files/listpix_ring_nin'//trim(adjustl(str1))//'_nout'//trim(adjustl(str2))//'_rad'//trim(adjustl(str3))//'.unf'
  !    open(unit2,file=fname,status='unknown',form='unformatted')

  !    do ispot=0,nspots-1

  !       write(strn,*)ispot
  !       open(unit,file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_ring/listpix_'//trim(adjustl(strn))//'.unf',status='unknown',form='unformatted')

  !       read(unit) nlist

  !       !print*, 'reading ispot, nlist: ',ispot, nlist

  !       allocate(listpix(nlist))

  !       read(unit) listpix
  !       close(unit)

  !       write(unit2) nlist
  !       write(unit2) listpix(0:nlist-1)

  !       deallocate(listpix)
  !       call system('rm '//'/mn/owl1/d3/yabebalf/planck/common_files/listpix_ring/listpix_'//trim(adjustl(strn))//'.unf')
        
  !    enddo

  !    close(unit2)

  ! endif





  CALL MPI_FINALIZE(ierr)
  if (me==0) print*, ' get_listpix done! '

end program get_dipdir

