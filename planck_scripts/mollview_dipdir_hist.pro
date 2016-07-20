pro mollview_dipdir_hist,mctheta,mcphi,nside=nside,fname=fname,fband=fband,map=map,bandmap=bandmap

nband = n_elements(mctheta[*,0])
nsim = n_elements(mctheta[0,*])

if not keyword_set(nside) then nside=8l

map = dblarr(nside2npix(nside))
bandmap = dblarr(nside2npix(nside),nband)

for ii=0,nband-1 do begin
      ang2pix_ring, nside,mctheta[ii,*], mcphi[ii,*], ipring

     for jj=0,nsim-1 do begin
        map[ipring[jj]] = map[ipring[jj]] + 1
        bandmap[ipring[jj],ii] = bandmap[ipring[jj]] + 1
     endfor

endfor


;;plot dipole directions 

if keyword_set(fname) then begin
   print, 'saving dipole direction histogram plot in:' 
   print, fname
endif

plot_mollview,map,min(map),max(map),file_ps=fname,ctitle='Dipole directions count',title=' ',/show_dipdir

if keyword_set(fband) and keyword_set(fname) then begin
   print, 'bin: saving dipole direction histogram plot in:'
   print, fband
   for ii=0,nband-1 do begin
      fname = fband+'_bin_'+strn(ii)+'.ps'
      map = bandmap[*,ii]
      plot_mollview,map,min(map),max(map),file_ps=fname,ctitle='Dipole directions count',title=' ',/show_dipdir
   endfor
endif

end
