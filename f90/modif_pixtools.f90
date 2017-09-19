module modif_pixtools

  use healpix_types
  use alm_tools
  use pix_tools
  USE long_intrinsic
  use misc_utils
  use extension

  IMPLICIT NONE
  INTEGER(KIND=i4b), private, PARAMETER :: ns_max4=8192     ! 2^13
  INTEGER(KIND=i4b), private, PARAMETER :: ns_max8=268435456! 2^28
  integer(KIND=i4b) ::proc_number
contains



! #ifdef DOI8B
!   subroutine query_disc_radvec_8( nside, vector0, radvec, listpix, nlist, nlistvec, nest, inclusive)
!     integer(i4b), parameter :: MKD = I8B
! #else
  subroutine query_disc_radvec ( nside, vector0, radvec, listpix, nlist, nlistvec, nest, inclusive)
    integer(i4b), parameter :: MKD = I4B
! #endif
    !Computes the pixels in a disk array. 
    !Fast because we loop over pixels only once for the big disk.
    !The smaller disks are determined along the way.
    !assumes the disks are ordered in ascending order
    !Yabebal Fantaye 27 Nov 2014
    !=======================================================================
    integer(kind=I4B), intent(in)                 :: nside
    integer(kind=I4B), intent(out)                 :: nlist
    real(kind=DP),     intent(in), dimension(1:)  :: vector0
    real(kind=DP),     intent(in), dimension(1:)  :: radvec
    integer(kind=MKD), intent(out), dimension(0:) :: listpix
    integer(kind=MKD), intent(out) ,dimension(1:) :: nlistvec
    integer(kind=I4B), intent(in), optional       :: nest
    integer(kind=I4B), intent(in), optional       :: inclusive

    INTEGER(KIND=i4B) :: nrad,imaxrad,irad,ngr
    INTEGER(kind=I4B) :: irmin, irmax, iz, ip, nir, nr, nrh
    REAL(kind=DP) :: norm_vect0
    REAL(kind=DP) :: z0, radius_eff, fudge
    REAL(kind=DP) :: phi0, dphi
    REAL(kind=DP) :: rlat0, rlat1, rlat2, zmin, zmax, z
    integer(kind=MKD) :: npix, list_size, nlost, ilist
    INTEGER(kind=I4B) :: status
    character(len=*), parameter :: code = "QUERY_DISC"
    logical(kind=LGT) :: do_inclusive
    integer(kind=I4B) :: my_nest
!    real(dp), allocatable, dimension(:) :: ztab, dphitab
    real(dp), dimension(1:4*nside-1) :: ztab, dphitab
    real(dp), dimension(:), allocatable :: zlist, dphilist
    integer(i4b), dimension(:,:),allocatable :: ringphi
    integer(i4b), dimension(:,:,:),allocatable :: allringphi
    integer(i4b), dimension(:),allocatable :: allngr
    integer(i4b) :: nsideh, nsboost
    real(dp) :: radiush, radius

    !=======================================================================

    nrad = size(radvec)
    radius=radvec(nrad)   !maxval(radvec)
    imaxrad=nrad  !maxloc(radvec)


    !allocate(ringphi(1:3, 1:4*nside-1))
    allocate(allringphi(1:3, 1:4*nside-1,1:nrad))
    allocate(allngr(nrad))

    list_size = long_size(listpix)
    !     ---------- check inputs ----------------
    npix = nside2npix(nside)

    if (radius < 0.0_dp .or. radius > PI) then
       write(unit=*,fmt="(a)") code//"> the angular radius is in RADIAN "
       write(unit=*,fmt="(a)") code//"> and should lie in [0,Pi] "
       call fatal_error("> program abort ")
    endif

    do_inclusive = .false.
    if (present(inclusive)) then
       if (inclusive == 1) do_inclusive = .true.
    endif

    my_nest = 0
    if (present(nest)) then
       if (nest == 0 .or. nest == 1) then
          my_nest = nest
       else
          print*,code//"> NEST should be 0 or 1"
          call fatal_error("> program abort ")
       endif
    endif

    ! radius_eff = radius
    ! if (do_inclusive) radius_eff = fudge_query_radius(nside, radius)

    !     ---------- circle center -------------
    norm_vect0 =  SQRT(DOT_PRODUCT(vector0,vector0))
    z0 = vector0(3) / norm_vect0

    do irad=1,nrad
       radius_eff = radvec(irad)
       if (do_inclusive) radius_eff = fudge_query_radius(nside, radvec(irad))
       
       !     --- coordinate z of highest and lowest points in the disc ---
       rlat0  = ASIN(z0)    ! latitude in RAD of the center
       rlat1  = rlat0 + radius_eff
       rlat2  = rlat0 - radius_eff
       if (rlat1 >=  halfpi) then
          zmax =  1.0_dp
       else
          zmax = SIN(rlat1)
       endif
       irmin = ring_num(nside, zmax)
       irmin = MAX(1, irmin - 1) ! start from a higher point, to be safe
       
       if (rlat2 <= -halfpi) then
          zmin = -1.0_dp
       else
          zmin = SIN(rlat2)
       endif
       irmax = ring_num(nside, zmin)
       irmax = MIN(4*nside-1, irmax + 1) ! go down to a lower point
       
       nr = irmax-irmin+1 ! in [1, 4*Nside-1]
       do iz = irmin, irmax
          ztab(iz-irmin+1) = ring2z(nside, iz)
       enddo


       call discphirange_at_z(vector0, radius_eff, ztab, nr, dphitab, phi0)
       call pixels_on_edge(nside, irmin, irmax, phi0, dphitab, allringphi(:,:,irad), allngr(irad),irad)

       ! allringphi(:,:,irad)=ringphi
       ! allngr(irad)=ngr

       if (do_inclusive) then
          ! sample edge pixels at larger Nside
          nsboost = 16
          nsideh = min(NS_MAX8, nside * int(nsboost,i8b))
          radiush = fudge_query_radius(nsideh, radius, quadratic=.true.)
          
          irmin = ring_num(nsideh, zmax)
          irmax = ring_num(nsideh, zmin)
          nrh = irmax - irmin + 1
          allocate(zlist(1:nrh), dphilist(1:nrh))
          do iz = irmin, irmax
             zlist(iz-irmin+1) = ring2z(nsideh, iz)
          enddo
          
          call discphirange_at_z(vector0, radiush, zlist, nrh, dphilist, phi0)
          call check_edge_pixels(nside, nsboost, irmin, irmax, phi0, dphilist, allringphi(:,:,irad), allngr(irad))
          ! allringphi(:,:,irad)=ringphi
          ! allngr(irad)=ngr
          deallocate(zlist, dphilist)
       endif

    enddo  !rad loop

    call discedge2fulldisc_vec(nside, allringphi, allngr, listpix, nlist, nlistvec,imaxrad)

    !deallocate(ringphi)
    deallocate(allringphi)
    deallocate(allngr)

    if (my_nest == 1) then
       do ip=0_MKD, nlist-1
          call ring2nest(nside, listpix(ip), listpix(ip))
       enddo
    endif



    
    return

! #ifdef DOI8B
!   end subroutine query_disc_radvec_8
! #else
  end subroutine query_disc_radvec
!#endif



  !=======================================================================
! #ifdef DOI8B
!   subroutine discedge2fulldisc_vec_8( nside, ringphi,ngr,list,nlist,nlistvec,ibigrad)
!     integer(i4b), parameter :: MKD = I8B
! #else
  subroutine discedge2fulldisc_vec(nside, ringphi, ngr, list,nlist,nlistvec,ibigrad)
    integer(i4b), parameter :: MKD = I4B
!#endif
    !It takes ringphi which is determined for radvec disk boundaries
    !ngr are the minimum and maximum ring numbers for the different disk boundaries
    !Modified by yabebal 27 Nov 2014
    !=======================================================================
    integer(i4b), intent(in) :: nside
    integer(i4b), dimension(1:,1:,1:), intent(in) :: ringphi
    integer(i4b), dimension(1:),        intent(in) :: ngr
    integer(MKD), dimension(1:),       intent(out) :: list
    integer(MKD), dimension(1:),       intent(out) :: nlistvec
    !
    integer(i4b),                                optional :: ibigrad
    !
    integer(i4b), dimension(:),  allocatable:: iradvec
    integer(i4b), dimension(:,:),allocatable:: templist,npradvec
    !
    integer(MKD), intent(out) :: nlist
    integer(MKD)              :: inlist,fnlist
    !
    integer(i4b) :: irad, nrad, imaxrad,iii,npix
    integer(i4b) :: my_alow,my_ahi, pr_alow,pr_ahi
    integer(i4b) :: kk,jj,ii, jdisk, pr_low, pr_hi,nprad,nplow,nphi
    !
    integer(i4b) :: j, jr_min, jr_max, nj, nr, kshift, np, my_low, my_hi,sjr_min, sjr_max,snj
    integer(i4b) :: ir, i, ip,ismall
    integer(i8b) :: npc, listsize
    logical(lgt) ::smalldisk
    !=======================================================================
    
    nrad=size(ngr)
    npix = size(list)

    allocate(templist(1:npix,nrad))
    allocate(npradvec(nrad,2))

    imaxrad=nrad !assuming the radvec is given in ascending order
    if (present(ibigrad)) imaxrad=ibigrad


    listsize = long_size(list)
    nlist = 0_MKD
    nlistvec = 0_MKD
    if (ngr(imaxrad) == 0) then ! no valid rings
       list(1) = -1
       return
    endif


    jr_min = ringphi(1, 1,imaxrad)
    jr_max = ringphi(1, ngr(imaxrad),imaxrad)
    nj = jr_max - jr_min + 1
    do j=0, nj-1
       ir = jr_min + j             ! current ring, in [1, nl4-1]

       !nr is maximum number of pixels on a ring
       !npc number of pixels on the current ring and above
       call pixels_per_ring(nside, ir, nr, kshift, npc)

       !my_low = ringphi(2, j+1,imaxrad)
       my_alow = ringphi(2, j+1,imaxrad)
       my_ahi = ringphi(3, j+1,imaxrad)

       if (my_ahi >= my_alow) then

          np = my_ahi - my_alow+1
          np = min(np, nr)  
          
          if (nlist + np > listsize) then
             print*,'Pixel query: too many pixels found for output list provided.'
             print*,'truncated at ',nlist
             return
          endif

         
          !find those disks which contains the current ring
          iii=0
          do irad=1,nrad-1
             sjr_min = ringphi(1, 1,irad)
             sjr_max = ringphi(1, ngr(irad),irad)
             if ( (ir>=sjr_min) .and. (ir<=sjr_max) ) then !height check                
                iii=iii+1
                jdisk = ir-sjr_min
                npradvec(iii,1)=irad
                npradvec(iii,2)=jdisk
                ! if (ir>=1 .and. ir<=57 .and. proc_number==8) print*,'small: ',irad,pr_alow,pr_ahi,nprad   !
             endif
          enddo
          ! !all pixels are of course in the big disk - at iii


          do i=0, np-1
             !get the smaller disk that contains the current pixel                              
             smalldisk=.false.             
             do kk=1,iii
                irad=npradvec(kk,1)
                jdisk=npradvec(kk,2)

                !raw hi and low
                pr_alow=ringphi(2,jdisk+1, irad)
                pr_ahi=ringphi(3, jdisk+1, irad)
                
                !Check if the current pixel belongs to the smaller disks
                if ( (my_alow+i >= pr_alow)  .and. (my_alow+i <= pr_ahi) ) then !side check
                   !small disk i 
                   ismall=my_alow + i - pr_alow
                   ip = modulo(pr_alow + ismall, nr)

                   nlistvec(irad)=nlistvec(irad) + 1
                   templist(nlistvec(irad),irad) = npc - nr + ip ! fill final list
                   smalldisk=.true.
                   exit  !at the first small disk that contains this pixel
                endif

             enddo !kk loop

             !big disk
             if (.not. smalldisk) then
                ip = modulo(my_alow + i, nr)
                nlistvec(imaxrad)=nlistvec(imaxrad) + 1
                templist(nlistvec(imaxrad),imaxrad) = npc - nr + ip ! fill final list
            endif
             
          enddo !ip loop


          nlist = nlist + np

       endif

    enddo !ir ring loop


    if (nlist == 0) then
       list(1) = -1
    else

       ip=1       
       np=nlistvec(1)
       list(ip:np)=templist(1:nlistvec(1),1)
       !if (proc_number==8) print*, 'ip, np',ip,np

       do irad=2,nrad          
          ip=sum(nlistvec(1:irad-1))+1
          np=sum(nlistvec(1:irad))          
          list(ip:np)=templist(1:nlistvec(irad),irad)          
          !if (proc_number==8) print*, 'irad,ip, np',irad,ip,np
       enddo

    endif
    
    deallocate(templist)
    deallocate(npradvec)

    Return

! #ifdef DOI8B
!   end subroutine discedge2fulldisc_vec_8

! #else
  end subroutine discedge2fulldisc_vec
!#endif



  !=======================================================================
  subroutine pixels_on_edge(nside, irmin, irmax, phi0, dphi, ringphi, ngr,irad)
  !=======================================================================
  !=======================================================================
    integer(i4b),                intent(in) :: nside, irmin, irmax
    real(dp),                    intent(in) :: phi0
    real(dp),     dimension(1:), intent(in) :: dphi
    integer(i4b), dimension(1:,1:), intent(out) :: ringphi
    integer(i4b),                   intent(out) :: ngr
    integer(i4b), parameter :: badvalue = -9999
    integer(i4b) :: nrings, ir, k, thisring, npr
    real(dp)    , parameter :: zero = 0.0_dp
    integer(i4b) :: iphi_low, iphi_hi, kshift,iii
    integer(i4b),optional :: irad
    real(dp) :: shift

    iii=1
    !if (present(irad)) iii=irad
    !if (proc_number==8) print*, 'thisring,irad, low, hi',thisring, irad, iphi_low, iphi_hi
    nrings = irmax - irmin + 1

    ngr = 0
    do thisring = irmin, irmax
       ir = thisring - irmin + 1
       call pixels_per_ring(nside, thisring, npr, kshift)

       shift = kshift * 0.5_dp
       iphi_low = ceiling (npr * (phi0 - dphi(ir)) / TWOPI - shift)
       iphi_hi  = floor   (npr * (phi0 + dphi(ir)) / TWOPI - shift)

       if (iphi_hi >= iphi_low) then ! pixel center in range
          ngr = ngr + 1
          ringphi(1, ngr) = thisring
          !ringphi(2, ngr) = modulo(iphi_low, npr)
          !ringphi(3, ngr) = modulo(iphi_hi,  npr)
          ringphi(2, ngr) = iphi_low
          ringphi(3, ngr) = iphi_hi
       endif
    enddo

  end subroutine pixels_on_edge
  !=======================================================================
  subroutine pixels_per_ring(nside, ring, npr, kshift, npnorth)
  !=======================================================================
    ! for a given Nside and ring index in [1,4*Nside-1], 
    ! returns the number of pixels in ring, their shift (0 or 1) in azimuth
    ! and the number of pixels in current ring and above (=North)
    !
    ! NB: 'rings' 0 and 4*Nside respectively are North and South Poles
    !=======================================================================
    integer(i4b), intent(in) :: nside, ring
    integer(i4b), intent(out) :: npr, kshift
    integer(i8b), intent(out), optional :: npnorth
    integer(i8b) :: ncap, npix, ir
    
    ! number of pixels in current ring
    npr = min(nside, ring, 4*nside-ring) * 4
    ! shift
    kshift = mod(ring + 1, 2) ! 1 for even, 0 for odd
    if (nside == 1) kshift = 1 - kshift ! except for Nside=1
    if (npr < 4*nside) kshift = 1 ! 1 on polar cap
    ! Number of pixels in current ring and above
    if (present(npnorth)) then
       if (ring <= nside) then ! in North cap
          npnorth = ring*(ring+1_i8b)*2_i8b
       elseif (ring <= 3*nside) then ! in Equatorial region
          ncap = nside*(nside+1_i8b)*2_i8b
          ir = ring-nside
          npnorth = ncap + 4_i8b*nside*ir
       else ! in South cap
          npix = nside2npix(nside)
          ir = 4_i8b*nside-ring - 1 ! count ring from south
          npnorth = npix - ir*(ir+1_i8b)*2_i8b
       endif
    endif
    return
  end subroutine pixels_per_ring
  !=======================================================================
  subroutine check_edge_pixels(nside, nsboost, irmin, irmax, phi0, dphi, ringphi, ngr)
  !=======================================================================
    integer(i4b), intent(in) :: nside, nsboost, irmin, irmax
    real(dp),                       intent(in) :: phi0
    real(dp),     dimension(1:),    intent(in) :: dphi
    integer(i4b), dimension(1:,1:), intent(inout) :: ringphi
    integer(i4b),                   intent(inout) :: ngr

    integer(i4b) :: i, j, k, kk, ngr_out, diff, iphi, i0
    real(dp), dimension(1:2*nsboost+1) :: phiw, phie
    real(dp) :: dd, dph, phic

  !=======================================================================
    if (nsboost <= 1) return

    do i=1, ngr ! loop on low-res rings
       i0 = ringphi(1,i) * nsboost - nsboost - irmin
       do k=-1,1,2 ! West and East side of disc
          kk = (k+5)/2 ! 2 or 3
222       continue
          iphi = ringphi(kk, i)
          if (ringphi(2,i) <= ringphi(3,i) .and. iphi >= 0) then
             call find_pixel_bounds(nside, nsboost, ringphi(1,i), iphi, phiw, phie)
             do j=1, 2*nsboost+1
                phic = (phie(i)+phiw(i))*0.5_dp ! pixel center
                dph  = (phie(i)-phiw(i))*0.5_dp + dphi(i0+j) ! pixel size + circle radius
                dd = abs(phi0 - phic)    ! distance from disc center to pixel border sample
                dd = min(dd, twopi - dd) ! in [0,Pi]
                if (dd <= dph) goto 1000 ! pixel touched by disc, move to next one
             enddo
             ringphi(kk, i)= iphi - k ! pixel not in disc, move edge pixel inwards
             goto 222 ! try next pixel inward
1000         continue
          endif
       enddo ! loop on side
    enddo ! loop on low-res rings

    ! remove empty rings
    ngr_out = 0
    do i=1,ngr
       diff = ringphi(3,i) - ringphi(2,i)
       if (ringphi(2,i) >=0 .and. ringphi(3,i) >=0 .and. diff /= -2 .and. diff /= -1) then
          ngr_out = ngr_out + 1
          ringphi(1:3, ngr_out) = ringphi(1:3, i)
       endif
    enddo
    ! set empty rings to -1
    do i=ngr_out+1, ngr
       ringphi(2:3, i) = -1
    enddo

    ngr = ngr_out
    return
  end subroutine check_edge_pixels


  !=======================================================================
  subroutine find_pixel_bounds (nside, nsboost, iring, iphi, phiw, phie)
    !=======================================================================
    integer(i4b),               intent(in)  :: nside, nsboost, iring, iphi
    real(dp),     dimension(1:2*nsboost+1), intent(out) :: phiw, phie
    
    real(dp),     dimension(1:2*nsboost+1) :: f, f1, phiw_t, phie_t
    real(dp) :: c0, quad, phie1, phie2, phiw1, phiw2, cv
    integer(i4b) :: npr, kshift, nq, ip, i
    logical(lgt) :: transition
  !=======================================================================

    call pixels_per_ring(nside, iring, npr, kshift)
    f = ((/ (i,i=0,2*nsboost) /) - nsboost) / nsboost

    nq = npr/4 ! number of pixels on current ring in [0,Pi/2] (quadrant)
    transition = (iring == nside .or. iring == nside*3)

    if (nq == nside .or. transition) then ! equatorial region (and transition rings)

       f1 = (1.0_dp-abs(f))*0.5_dp    ! triangle of height 1/2
       f1 = halfpi * f1 / nq
       c0 = halfpi * (iphi + kshift*0.5_dp) / nq
       phiw = c0 - f1
       phie = c0 + f1
       if (transition) then ! store for future use
          phiw_t = phiw
          phie_t = phie
       endif
    endif

    if (nq < nside .or. transition) then ! polar regions and transition rings
       ip = mod(iphi,nq) ! in [0,nq-1]
       quad = iphi / nq ! quadrant in [0,3]
       if (iring <= nside*2) then
          f1 = halfpi / (nq + f) 
       else
          f1 = halfpi / (nq - f)! swap sign for South pole
       endif
       do i=1, 2*nsboost+1
          cv = f1(i)
          phiw1 = min(cv *     ip,    halfpi)
          phie1 = min(cv * (   ip+1), halfpi)
          phiw2 = min(cv * (nq-ip-1), halfpi)
          phie2 = min(cv * (nq-ip),   halfpi)
          phiw(i) = max(phiw1, halfpi - phie2) + (quad * halfpi)
          phie(i) = min(phie1, halfpi - phiw2) + (quad * halfpi)
       enddo
    endif

    if (transition) then 
       if (iring == nside) then ! transition in N hemisphere
          phiw(nsboost+2:2*nsboost+1) = phiw_t(nsboost+2:2*nsboost+1) 
          phie(nsboost+2:2*nsboost+1) = phie_t(nsboost+2:2*nsboost+1)
       else ! transition in S hemisphere
          phiw(1:nsboost+1) = phiw_t(1:nsboost+1)
          phie(1:nsboost+1) = phie_t(1:nsboost+1)
       endif
    endif

    return
  end subroutine find_pixel_bounds


end module
