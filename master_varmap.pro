

;;;IDL extension to HEALPix [TESTVERSION]: Written by Frode K. Hansen, Nov. 2006
function master_varmap,kern,varmap,lmax=lmax_in,bins=bins_in,beam1=beam_in,fwhm1=fwhm_in,silent=silent,beam2=beam2_in,$
                       fwhm2=fwhm2_in,nproc=nproc,mask=mask_in,ncl=ncl

  dir='./temp/'
  direxe='./'
  bindir='/mn/owl1/d3/yabebalf/planck/bin/'
  spawn,'mkdir -p '+dir

  unit=99
  status = FSTAT(unit)
  if (status.open ne 0) then begin
     WHILE(status.open ne 0) do begin
        unit=unit-1
        status = FSTAT(unit)
     ENDWHILE
  endif

  routine = 'master_varmap'
  uroutine = strupcase(routine)
  if (n_params() ne 2) then begin
     PRINT, 'Wrong number of arguments in ' + uroutine
     print,'Syntax : '
     print,uroutine+', kern, map '
     print,'              [LMAX=, BEAM1=, FWHM1=, BINS= ,BEAM2= ,FWHM2= ,MASK= ,/SILENT]'
     return,0
  endif

  if ((KEYWORD_SET(fwhm_in)) or (KEYWORD_SET(beam_in))) then begin
     if (NOT (KEYWORD_SET(silent))) then begin
        print,'ASSUMING that FWHM and BEAM was NOT set in make_kern! If FWHM or BEAM was already specified in make_kern, it SHOULD NOT be repeated here, this would result in a double deconvolution!'
     endif
  endif


  lmax=20
  if (KEYWORD_SET(lmax_in)) then lmax=LONG(lmax_in)
  polar=1l

  nsim=n_elements(varmap[0,*])
  npix=n_elements(varmap[*,0])
  nside=npix2nside(npix)

  ;;binning
  if (KEYWORD_SET(bins_in)) then begin
     nbins=LONG(n_elements(bins_in(*,0)))
  endif else begin
     nbins=LONG(lmax)-1
  endelse

  ;;binning kernel
  binkern=0l
  if (KEYWORD_SET(bins_in)) then begin
     binkern=1l
     sbins=SIZE(bins_in)
     sbins=LONG(sbins[0])
     if (sbins ne 1l) then begin
        print,'bins must be on the form bins=lonarr(nbins)'
        return,0
     endif

     allbins=dblarr(nbins)
     allbins(*)=bins_in(0l:nbins-1l)

     openw,unit,dir+'nbins.unf',/f77
     writeu,unit,LONG(nbins)
     close,unit

     openw,unit,dir+'binvec.unf',/f77
     writeu,unit,LONG(allbins)
     close,unit
  endif



  if ((n_elements(reform(kern(*,0))) ne nbins) or (n_elements(reform(kern(*,0))) ne n_elements(reform(kern(0,*))))) then begin
     print,'kern must be of the form kern=dblarr(nbins,nbins) or (lmax-1,lmax-1)'
     return,0
  endif



  fwhm0=0.
  if (KEYWORD_SET(fwhm_in)) then fwhm0=FLOAT(fwhm_in)
  beam0=dblarr(lmax+1)
  beam0(*)=1d0

  if (KEYWORD_SET(beam_in)) then beam0(*)=beam_in(0:lmax)
  sigma=fwhm0/60d0/180d0*!pi/sqrt(8d0*alog(2d0))

  if not keyword_set(ncl) then ncl=dblarr(lmax+1)

  openw,unit,dir+'ncl.unf',/f77
  writeu,unit,double(ncl)
  close,unit

  openw,unit,dir+'lmax.unf',/f77
  writeu,unit,LONG(lmax)
  close,unit

  openw,unit,dir+'binkern.unf',/f77
  writeu,unit,binkern
  close,unit

  help, kern
  openw,unit,dir+'kern.unf',/f77
  writeu,unit,double(kern(0l:nbins-1l,0l:nbins-1l))
  close,unit

  openw,unit,dir+'nsim.unf',/f77
  writeu,unit,LONG(nsim)
  close,unit

  openw,unit,dir+'nside.unf',/f77
  writeu,unit,LONG(nside)
  close,unit

  help, varmap
  openw,unit,dir+'varmap.unf',/f77
  writeu,unit,double(varmap)
  close,unit


  ;;define mask
  if (keyword_set(mask_in)) then begin
     mask=mask_in
     npx=n_elements(mask)
     if (npx ne nside^2*12l) then begin
        print,'size of mask does not agree with nside!'
        return,0
     endif
  endif else begin
     mask=fltarr(nside^2*12l)
     mask(*)=1.
  endelse

  help, mask
  openw,unit,dir+'mask.unf',/f77
  writeu,unit,float(mask)
  close,unit


  if not keyword_set(nproc) then nproc=min([125,nsim])

  ;;run code
  f90code=bindir+'master_varmap'
  str='cd '+direxe+'; time mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME  '+f90code
  print, '************************* spawning **********************'
  print, '** ', str, ' **'
  print, '*********************************************************'


  ;;load res 
  res=dblarr(nbins,nsim) 
  openr,unit,dir+'res.unf',/f77
  readu,unit,res
  close,unit


  return,res

end
