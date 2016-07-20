pro get_clratio,res,mcres,droot=droot,nband=nband,xth=xth,xph=xph,$
                 out_struct=out_struct,nrep=nrep,iout=iout,mres=mres,dres=dres,umcres=umcres,ures=ures,$
                 nnbin=nnbin,nproc=nproc,sdisk=sdisk,overwrite=overwrite,nosum=nosum,dth=dth,dph=dph

;;read spectra from data:
  if not keyword_set(iout) then iout=0
  if not keyword_set(nrep) then nrep=1
  if not keyword_set(nband) then nband=15
  if not keyword_set(nnbin) then nnbin=93
  if not keyword_set(nproc) then nproc=40
  if not keyword_set(sdisk) then sdisk=90d0

  lastb=nband-1

  nspots = n_elements(res[0,*])
  nbins = n_elements(res[*,0])
  nsim = 1
  help, mcres

  if keyword_set(mcres) then begin
     nsim = n_elements(mcres[0,0,*])
     mcres2=mcres
     ;;a second mcres to use for obtaining dipole directions for the simulation
     if keyword_set(umcres) then begin
        print, '**** using umcres ***'
        mcres2=umcres
     endif
  endif


  bins = lindgen(nbins)*16l+2l

  ;;a second res to use for obtaining dipole directions for the data
  res2=res
  if keyword_set(ures) then begin
     print, '**** using ures ***'
     res2=ures
  endif

  if not (keyword_set(xph) or keyword_set(xth)) and keyword_set(mcres2) then begin

     if not (keyword_set(mres) and keyword_set(dres)) then begin 

        nspots2 = n_elements(res2[0,*])
        nbins2 = n_elements(res2[*,0])
        if keyword_set(mcres) then nsim = n_elements(mcres[0,0,*])
        mres=dblarr(nbins2,nspots2)
        dres=mres
        for i=0,nbins2-1 do begin
           if i mod 20 eq 0 then print,i
           for j=0,nspots2-1 do begin
              mres(i,j)=mean(mcres(i,j,0:nsim-1))
              dres(i,j)=stdev(mcres(i,j,0:nsim-1))
           endfor
        endfor
     endif

     if keyword_set(droot) then spawn, 'mkdir -p '+droot
     if not keyword_set(droot) then   overwrite=1
     if not keyword_set(droot) then   droot='./'

     fth=droot+'theta1_nnbin'+strn(nnbin)+'.unf'
     fph=droot+'phi1_nnbin'+strn(nnbin)+'.unf'

     if not file_test(fth,/regular) or keyword_set(overwrite) then begin 
        get_cross_dipdir,bins,res2,dres,mcres2,mres,$
                         nband=1,nnbin=nnbin,nproc=nproc,$
                         theta=xth,phi=xph,bres=bres,dipamp=dipamp
        wunf, xth,fth
        wunf, xph,fph
     endif else begin
        xth = dblarr(nsim+1)
        xph=xth

        runf, xth,fth
        runf, xph,fph
     endelse

     print, 'data (th, ph) = (',xth[0]*!rtod,xph[0]*!rtod,')'
     print, 'mcres 2 (th, ph) = (',xth[3]*!rtod,xph[3]*!rtod,')'
 
  endif


  nsim=n_elements(xth)-1

  n16bins=nbins

  ddx9 = dblarr(nbins,3)
  bddx9 = dblarr(nband,3)
 
  bmcres=dblarr(nband,3,nsim)
  mbmcres=dblarr(nband,3)
  dbmcres=dblarr(nband,3)

  
  ddx9res=dblarr(n16bins,nspots)


  if (keyword_set(dth) and keyword_set(dph)) then begin
     x_thph = make_thph_mat(dth,dph)
  endif else begin
     x_thph = make_thph_mat(xth[0],xph[0])
  endelse


  print, 'data thph: ', x_thph*!rtod


  nside_map=npix2nside(nspots)

  if keyword_set(nosum) then begin
     ang2pix_ring, npix2nside(nspots), x_thph[0,0], x_thph[0,1], ind
     ang2pix_ring, npix2nside(nspots), x_thph[1,0], x_thph[1,1], indn
  endif else begin
     map = make_disk_to_mollview(x_thph[0,0], x_thph[0,1],[2,2],sdisk,nside_map, bad_value=bad_value)
     ind = where(map gt  bad_value+1)
     
     ;;negative direction
     map = make_disk_to_mollview(x_thph[1,0], x_thph[1,1],[2,2],sdisk,nside_map, bad_value=bad_value)
     indn = where(map gt  bad_value+1)
  endelse

   weight_ind = 1d0+0d0/dres[0:n16bins-1,ind]^2
   weight_indn = 1d0+0d0/dres[0:n16bins-1,indn]^2

  print, 'positive direction: averaging over nspots=',n_elements(ind) ;ind
  print, '---------------------'
  print, 'negative direction: averaging over nspots=',n_elements(indn) ;indn

;;get dx9 power spectra averaged on two hemisphere (disks)
 if n_elements(ind) gt 1 then  ddx9[0:n16bins-1,0] = total(res[0:n16bins-1,ind]*weight_ind,2,/double)/total(weight_ind,2,/double) else  ddx9[0:n16bins-1,0] = res[0:n16bins-1,ind]
 if n_elements(indn) gt 1 then ddx9[0:n16bins-1,1] = total(res[0:n16bins-1,indn]*weight_indn,2,/double)/total(weight_ind,2,/double) else ddx9[0:n16bins-1,1] = res[0:n16bins-1,indn]

 nbin2bin=6l

;;binning data cl and construct ratio
  bin_res, bins, ddx9[*,*], xx,bbins=bbins,nnbin=nbin2bin,nband=nband
  bddx9[*,*] = xx
  xx = bddx9[*,0]
  yy = bddx9[*,1]
  bddx9[*,2] = (xx-yy)/xx

  if size(out_struct,/type) eq 8 then begin
     print,'Out_struct is defined: saving data clratio to the old struc'
     out_struct[iout].bddx9 = bddx9
  endif


print, '****** nsim, n_elements(xth)', nsim, n_elements(xth)
  if nsim gt 10 then  begin
     for iisim=0,nsim-1 do begin
        
        

        if n_elements(xth) eq nsim+1 then x_thph = make_thph_mat(xth[iisim+1],xph[iisim+1])
        ;;print, iisim,x_thph

        if keyword_set(nosum) then begin
           ang2pix_ring, npix2nside(nspots), x_thph[0,0], x_thph[0,1], ind
           ang2pix_ring, npix2nside(nspots), x_thph[1,0], x_thph[1,1], indn
        endif else begin
           map = make_disk_to_mollview(x_thph[0,0], x_thph[0,1],[2,2],sdisk,nside_map, bad_value=bad_value)
           ind = where(map gt  bad_value+1)
           ;;negative direction
           map = make_disk_to_mollview(x_thph[1,0], x_thph[1,1],[2,2],sdisk,nside_map, bad_value=bad_value)
           indn = where(map gt  bad_value+1)
        endelse

        ;;inverse variance weighting 
        weight_ind = 1d0+0d0/dres[0:n16bins-1,ind]^2
        weight_indn = 1d0+0d0/dres[0:n16bins-1,indn]^2

        if n_elements(ind) gt 1 then  bin_res, bins, total(mcres[*,ind,iisim]*weight_ind,2,/double)/total(weight_ind,2,/double), xx,bbins=bbins,nnbin=nbin2bin,nband=nband $
        else  bin_res, bins, mcres[*,ind,iisim], xx,bbins=bbins,nnbin=nbin2bin,nband=nband

        bmcres[*,0,iisim] = xx

        ;;negative direction
        if n_elements(ind) gt 1 then  bin_res, bins, total(mcres[*,indn,iisim]*weight_indn,2,/double)/total(weight_indn,2,/double), yy,bbins=bbins,nnbin=nbin2bin,nband=nband $
        else  bin_res, bins, mcres[*,indn,iisim], yy,bbins=bbins,nnbin=nbin2bin,nband=nband

        bmcres[*,1,iisim] = yy
        
        bmcres[*,2,iisim] = (xx-yy)/xx
     endfor

;;get the mean and variance of mcres_hemisphere
     for i=0,nband-1 do begin
        
        for j=0,2 do begin
           mbmcres[i,j]=mean(bmcres[i,j,*])
           dbmcres[i,j]=stdev(bmcres[i,j,*])
        endfor
        
     endfor

     if size(out_struct,/type) eq 8 then begin
        print,'Out_struct is defined, and nsim>10: save to the old struc'
        out_struct[iout].bmcres = bmcres
        out_struct[iout].mbmcres=mbmcres
        out_struct[iout].dbmcres=dbmcres
        out_struct[iout].xth=xth
        out_struct[iout].xph=xph
     endif


  endif


  if size(out_struct,/type) ne 8 and (nsim gt 10) then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim>10'
     out_struct = {bmcres:bmcres,mbmcres:mbmcres,dbmcres:dbmcres,xth:xth,xph:xph,bddx9:bddx9}
     out_struct = replicate(out_struct,nrep)
  endif 

  if size(out_struct,/type) ne 8 and n_params() eq 1 then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim<10'
     out_struct = {bddx9:bddx9}
     out_struct = replicate(out_struct,nrep)
  endif 



end
