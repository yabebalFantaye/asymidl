pro pixel_pairs, nside,map,map_out,theta=theta,phi=phi,nside_out=nside_out,order=order,mnval=mnval,msval=msval,$
                 mval=mval,ratio=ratio,ppair=ppair,prank=prank,irank=irank,verbose=verbose

if n_params() lt 1 then begin
   print, 'Usage: '
   print, 'pixel_pairs, nside,map,map_out,theta=theta,phi=phi,nside_out=nside_out,order=order,mpair=mpair,$'
   print, '    mval=mval,diff=diff,ppair=ppair,prank=prank,irank=irank,verbose=verbose'
   return
endif

  !quiet=1


if not keyword_set(nside_out) then nside_out = nside
Nside_out = nside_out


pixels = lindgen(nside2npix(Nside_out))

;print, 'nside_out= ',nside_out,nside
if not keyword_set(order) then order='RING'

  if keyword_set(verbose) then print, ' generate the pixel-antipixel correspondence in the ring scheme'
  ;; (since the generated maps are in this scheme)
  if (order eq 'NESTED') then begin
     ;; pix2vec_nest, Nside_out, pixels, vectors
     ;; anti_vectors=-vectors
     ;; vec2pix_nest, Nside_out, anti_vectors, anti_pixels
     print,'ERROR: pixel pair for nested ordering needs attention!'
     on_error,2
  endif else begin
     if (order eq 'RING') then begin
        pix2vec_ring, Nside_out, pixels, vectors
        anti_vectors=-vectors
        vec2pix_ring, Nside_out, anti_vectors, anti_pixels
     endif
  endelse

npix=nside2npix(Nside_out)
upixels = pixels(0:npix/2-1)
anti_upixels = anti_pixels(0:npix/2-1)

ppair=[[upixels],[anti_upixels]]


prank = upixels

if n_params() gt 1 then begin

   if keyword_set(verbose) then print, 'prank contains ranking of abs of opposite pixel differences in descending order!'

  ;;read_fits_map, filename, map
  map_deg = map
  if (nside ne nside_out) then ud_grade, map, map_deg, nside_out=Nside_out, order_in=order,order_out=order

  if not keyword_set(ratio) then begin
     ;;print, ' generating map of differences for opposite pixels'
     map_out=(map_deg-map_deg[anti_pixels])
  endif else begin
     ;;print, ' generating map of differences-ratio for opposite pixels'
     map_out=(map_deg-map_deg[anti_pixels])*100d0/map_deg
     ;;map_out=map_deg/map_deg[anti_pixels]   ;*100d0/map_deg
  endelse

  prank=reverse(sort(abs(map_out[upixels]))) ;;now map_out is sorted in decreasing order


  ;;if theta and phi are given compute irank
  ;;else use irank keyword or set irank=1 to get map_value at irank

  if keyword_set(theta) and keyword_set(phi) then begin
     ang2pix_ring, nside_out, theta, phi, ipix
     irank=ipix
     for ii=0,n_elements(ipix)-1 do begin
        if ipix[ii] lt npix/2 then irank[ii] = where(ipix[ii] eq upixels[prank])+1 
        if ipix[ii] gt npix/2 then irank[ii] = where(ipix[ii] eq anti_upixels[prank])+1 
     endfor
    if keyword_set(verbose) then  print, 'npix, ipix, irank',npix, ipix,irank
  endif else begin
     if not keyword_set(irank) then irank=1
  endelse


  mval=abs(map_out[upixels[prank[irank-1,*]]]) ;;the ith rank value fo the map
  mnval = map[upixels[prank[irank-1,*]]]
  msval = map[anti_upixels[prank[irank-1,*]]]

endif

end
