function angular_dipdir_hist,mctheta,mcphi,nside=nside,fname=fname,map=map,bandmap=bandmap,thin=thin,phin=phin,nrad=nrad,rcount=rcount,bandrcount=bandrcount,rad=rad,rval=rval

  if n_params() eq 0 then print, 'Usage: rcount=angular_dipdir_hist(mctheta,mcphi,nside=nside,fname=fname,map=map,bandmap=bandmap,thin=thin,phin=phin,nrad=nrad,rcount=rcount,bandrcount=bandrcount,rad=rad)'


  if not keyword_set(phin) then  phin=225d0*!dtor
  if not keyword_set(thin) then  thin=89d0*!dtor
  if not keyword_set(nrad) then  nrad=20

  nband = n_elements(mctheta[*,0])
  nsim = n_elements(mctheta[0,*])

  if not keyword_set(nside) then nside=8l

  map = intarr(nside2npix(nside))
  bandmap = intarr(nside2npix(nside),nband)

  for ii=0,nband-1 do begin
     ang2pix_ring, nside,mctheta[ii,*], mcphi[ii,*], ipring

     for jj=0,nsim-1 do begin
        map[ipring[jj]] = map[ipring[jj]] + 1
        bandmap[ipring[jj],ii] = bandmap[ipring[jj]] + 1
     endfor

  endfor


  val=2.

  if not keyword_set(rval) then rval=linspace(0,90,nrad)

  rcount=intarr(nrad)
  bandrcount=intarr(nrad,nband)
  rad=fltarr(nrad)

  for i=0,nrad-1 do begin
     r=rval[i]

     mapx=make_disk_to_mollview(thin,phin,val,r,nside,bad_value=bad_value)
     ind = where(mapx gt  bad_value+1)
     ;; print, 'r, count',round(r), total(map[ind])

     rcount[i]=total(map[ind])
     bandrcount[i,*]=total(bandmap[ind,*],1)
     ;;subtract points that fall in the previous rings
     if i ne 0 then begin
        rcount[i]=rcount[i]-total(rcount[0:i-1])
        bandrcount[i,*]=bandrcount[i,*] - total(bandrcount[0:i-1,*],1)
     endif

     rad[i]=r
  endfor

  return,rcount

end
