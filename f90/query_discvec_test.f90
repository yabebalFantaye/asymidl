  program test


  use healpix_types
  use alm_tools
  use pix_tools
  use extension

  use fitstools
  use head_fits

  use modif_pixtools
  USE long_intrinsic, only: long_size
  use misc_utils,only: assert, assert_alloc, fatal_error, wall_clock_time, string, strupcase, brag_openmp

   IMPLICIT  NONE
   integer(i4b), parameter :: nrad=16
   real(dp), dimension(1:nrad) :: radvec

   real(sp), dimension(:,:), allocatable :: mask

   real(dp) :: rad, vector0(3), theta, phi
   integer(i4b) :: ispot, unit
   integer(i4b) :: nlist,nside,nside_spots,nspots, npix
   integer(i4b) :: j, ip, np, k, ntotlist, irad
   integer(i4b), dimension(0:nrad-1) :: nlistvec
   integer(i4b), dimension(nrad,2) :: npixvalid
   integer(i4b), dimension(:), allocatable :: listpix,allrad_listpix

  !timing 
  real     (SP) :: clock_time, time0, time1
  real     (SP) :: ptime, ptime0, ptime1

  CHARACTER(LEN=500) :: numstr,filenum,filename

   !radvec = (/1,2,4,6,10/)*DEG2RAD
  radvec=(/1,2,4,6,8,10,12,14,16,18,20,22,30,50,70,90/)*DEG2RAD
   print*,'----------------'
   print*, 'radvec: ',radvec
   print*,'----------------'

   unit=12

   nside = 2048
   nside_spots=16
   npix=nside2npix(nside)
   nlist=npix
   
   allocate(mask(0:npix-1,1))
   filename = '/mn/owl1/d3/yabebalf/planck/DX11/dx11_v2/dx11_v2_common_int_mask_005a_2048.fits'
   call input_map(filename,mask,npix,1)   


   allocate(listpix(0:npix-1),allrad_listpix(0:npix-1))

  call wall_clock_time(time0)
  call cpu_time(ptime0)

  ispot=1240
  !call pix2vec_ring(nside_spots,ispot,vector0)
  !call query_disc_radvec ( nside, vector0, radvec, allrad_listpix, ntotlist,nlistvec)

  filename='/mn/owl1/d3/yabebalf/planck/common_files/listpix_ring_nin2048_nout16_nrad16/listpix_1240.unf'
  open(unit,file=trim(filename),status='old',form='unformatted')
  read(unit) k, ntotlist
  read(unit) nlistvec
  read(unit) allrad_listpix(0:ntotlist-1)
  close(unit)

print*, k

   print*,'==============='
   print*, ''
   print*, 'modified: nlist, sum(nlistvec): ',ntotlist, sum(nlistvec)
   print*, ''
   print*,'==============='



   k=0
   ip=1
   ntotlist=0
   
   do irad = 1,nrad
      
      ip=ntotlist !previous nlist
      ntotlist=sum(nlistvec(0:irad-1))
      npixvalid(irad,1) = ntotlist

      !now listpix contains only those pixels that are to be evaluated      
      do j=ip,ntotlist-1
         if ( mask(allrad_listpix(j),1) .gt. 0.8 ) then 
            listpix(k)=allrad_listpix(j)
            k=k+1
         endif
      enddo
      nlist=k !number of valid pixels
      
      npixvalid(irad,2)=nlist

      if (irad==1 .or. irad==5) then
         print*,'==============='
         print*, ''
         print*, 'original: rad, nlist,nvalid, mean(listpix/nvalid): ',radvec(1),ntotlist,dble(nlist), sum(dble(listpix(0:k-1))/dble(k))
         print*, ''
         print*,'==============='
      endif

   enddo


  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0

  print*,'==============='
  print*, ''
  write(*,*) " modified: Clock and CPU time [s]: ", clock_time, ptime
  print*, ''
  print*,'==============='


   print*,'==============='
   print*, ''
   print*, 'modified: nlist, sum(nlistvec): ',nlist, sum(nlistvec)
   print*, 'nlistvec',nlistvec
   print*, 'allrad_ntotlist: ',npixvalid(:,1)
   print*, 'allrad_nvalidlist: ',npixvalid(:,2)
   print*, ''
   print*,'==============='


   !===================================================
   !-----------------------------------------------------
   !===================================================

  call wall_clock_time(time0)
  call cpu_time(ptime0)


  !disk 1
  call query_disc(nside, vector0, radvec(1), allrad_listpix, ntotlist)
  !now listpix contains only those pixels that are to be evaluated      
  k=0
  ip=0
  do j=ip,ntotlist-1
     if ( mask(allrad_listpix(j),1) .gt. 0.8 ) then 
        listpix(k)=allrad_listpix(j)
        k=k+1
     endif
  enddo
  print*,'==============='
  print*, ''
  print*, 'original: rad, nlist,nvalid, mean(listpix/nvalid): ',radvec(1),ntotlist,k, sum(dble(listpix(0:k-1))/dble(k))
  print*, ''
  print*,'==============='



  ! disk 5
  call query_disc(nside, vector0, radvec(5), allrad_listpix, ntotlist)
  k=0
  ip=0
  do j=ip,ntotlist-1
     if ( mask(allrad_listpix(j),1) .gt. 0.8 ) then 
        listpix(k)=allrad_listpix(j)
        k=k+1
     endif
  enddo
  print*,'==============='
  print*, ''
  print*, 'original: rad, nlist,nvalid, mean(listpix/nvalid): ',radvec(5),ntotlist,k, sum(dble(listpix(0:k-1))/dble(k))
  print*, ''
  print*,'==============='


  call wall_clock_time(time1)
  call cpu_time(ptime1)
  clock_time = time1 - time0
  ptime      = ptime1 - ptime0
  write(*,*) " original: Clock and CPU time [s]: ", clock_time, ptime
   !----------------------------------

  

end   program test
