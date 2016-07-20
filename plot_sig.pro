pro plot_sig, sigmat, fname=fname,ires=ires,colvec=colvec,$
              xra=xra, yra=yra, nx=nx,nmx=nmx,nmy=nmy,ny=ny,nflmax=nflmax,nflmin=nflmin,legtxt=legtxt,lminmax=lminmax,shift=shift,nothph=nothph,xval=xval,circol=circol,cirsz=cirsz,_extra=extra


nothph=1
flmax=1
flmin=1
if keyword_set(nflmin) then flmin=0
if keyword_set(nflmax) then flmax=0

 if not keyword_set(cirsz) then cirsz=0.4
  if not keyword_set(xra) then xra = [0,1500]
  ;;if not keyword_set(yra) then yra = [-10,15]
  if not keyword_set(ny) then ny = 5
  if not keyword_set(nx) then nx = 3
  if not keyword_set(nmx) then nmx=5
  if not keyword_set(nmy) then nmy=5
  if not keyword_set(ires) then ires=0
  
print, 'ny,nym,nx, nxm',ny,nmy,nx, nmx

if keyword_set(flmin) and keyword_set(flmax) then lminlmax='min/max' & ltxt=['lmin-1500','2-lmax']
if keyword_set(flmin) and not keyword_set(flmax) then lminlmax='min' & ltxt=['lmin-1500','']
if not keyword_set(flmin) and keyword_set(flmax) then lminlmax='max' & ltxt=['2-lmax','']

if not keyword_set(legtxt) then legtxt=ltxt
if keyword_set(lminmax) then lminlmax=lminmax

llmax=xra[1]

  nband = n_elements(sigmat[*,0,0])
  nsim = n_elements(sigmat[0,0,*])
  ncase = n_elements(sigmat[0,*,0])

  nbin2bin=6l
  dl=100
  dl_cent=50
  if nband gt 90l then begin
     nbin2bin=1l
     dl=16
     dl_cent=8
  endif
  if (nband eq 1) then nbin2bin=floor(llmax/(16l*nband))
  lnew = indgen(nband)*dl+dl_cent
  val=lnew

if keyword_set(xval) then lnew=xval

if not keyword_set(shift) then dl=0

;;sig is in percentage
  sig=sigmat

print, 'nband, ncase, nsim: ',nband,ncase,nsim

  if nsim gt 1 then begin
     msig= dblarr(nband,ncase)
     dsig=msig
     for i=0,nband-2 do begin
        for j=0,3 do begin
           msig(i,j)=mean(sig(i,j,*))
           dsig(i,j)=stdev(sig(i,j,*))
        endfor
     endfor
  endif else begin
     msig=sig
     dsig=sig
  endelse


;;======================
;;------ Significances plot --
;;===================


  planck_ct,cvec

  if not keyword_set(colvec) then colvec=cvec


if not keyword_set(yra) then begin
  if nsim gt 1 then begin
     yra = round([min(msig)-max(dsig),max(msig)+max(dsig)])
     yra = [floor(-yra[1]/(ny-1.)),(ny-1)*ceil(yra[1]/(ny-1.))]
  endif else begin
     yra = round([min(msig)-1.2*min(msig),max(msig)+1.2*max(msig)])
     yra = [floor(-yra[1]/(ny-1.)),(ny-1)*ceil(yra[1]/(ny-1.))]
  endelse
endif


one_percent=1
if max( msig[0:nband-1,0]) gt 50 then one_percent=99

;;-------------------------------
;;plot theta sum 
;;-------------------------------
  if keyword_set(fname) then fcase=fname+'clustering.ps'

  planck_plot,lnew[0:nband-1]+dl,msig[0:nband-1,0],color=0,fname=fcase, $
              xra=xra, yra=yra, thick=2,xt='xttl', yt='yttl', $ ;, XTTL_DY=xtdy
              /noeras,/nodata, nx=nx,nmx=nmx,nmy=nmy,ny=ny,xs=xsize,ys=ysize,_extra=extra


  ls_circle,sz=0.4,color=circol

  if nsim gt 1 then begin
     if keyword_set(flmax) then oploterror, lnew[0:nband-1]+dl, msig[0:nband-1,0], dsig[0:nband-1,0],color=colvec[0],thick=2,psym=-8,errcolor=colvec[0]
     if keyword_set(flmin) then oploterror, lnew[0:nband-1]+dl, msig[0:nband-1,1],dsig[0:nband-1,1], color=colvec[1],thick=2,psym=-8,errcolor=colvec[1]
  endif else begin
     if keyword_set(flmax) then oplot, lnew[0:nband-1]+dl, msig[0:nband-1,0],color=colvec[0],thick=2,psym=-8
     if keyword_set(flmin) then oplot, lnew[0:nband-1]+dl, msig[0:nband-1,1], color=colvec[1],thick=2,psym=-8
  endelse
  oplot,indgen(llmax),0*indgen(llmax)+one_percent,line=2,thick=1,color=0


  if keyword_set(flmin) and keyword_set(flmax) then begin
     if max( msig[0:nband-1,0]) gt 50 then begin
        legend,legtxt,color=colvec[[0,1]],textcolor=[0,0],/bottom,/right,box=0,$
               line=[0,0],pspacing=1,thick=4,charsize=lcharsize
     endif else begin
        legend,legtxt,color=colvec[[0,1]],textcolor=[0,0],/top,/right,box=0,$
               line=[0,0],pspacing=1,thick=4,charsize=lcharsize
     endelse
  endif

  ax=['xttl','yttl']
  ay=['$\ell_{\rm '+lminlmax+'}$','Significances of clustering $[\%]$']

  if keyword_set(fname) then ps_end_planck, /png,feps=fcase,/latex,/cl,axisx=ax,axisy=ay,xsize=xsize,ysize=ysize

;;-------------------------------
;;plot cos(theta) sum 
;;-------------------------------

if not keyword_set(nothph) then begin

  if keyword_set(fname)  then fcase=fname+'cluster_wrt_thph.ps'

  planck_plot,lnew[0:nband-1],msig[0:nband-1,0],color=0,fname=fcase, $
              xra=xra, yra=yra, thick=2,xt='xttl', yt='yttl', $ ;, XTTL_DY=xtdy
              /noeras,/nodata, nx=nx,nmx=nmx,nmy=nmy,ny=ny,xs=xsize,ys=ysize,_extra=extra


  LS_circle,sz=cirsz,color=circol

  if nsim gt 1 then begin
      if keyword_set(flmax) then oploterror, lnew[0:nband-1]+dl, msig[0:nband-1,2], dsig[0:nband-1,2],color=colvec[0],thick=2,psym=-8,errcolor=colvec[0]
      if keyword_set(flmin) then oploterror, lnew[0:nband-1]-dl, msig[0:nband-1,3],dsig[0:nband-1,3], color=colvec[1],thick=2,psym=-8,errcolor=colvec[1]
  endif else begin
      if keyword_set(flmax) then oplot, lnew[0:nband-1]+dl, msig[0:nband-1,2],color=colvec[0],thick=2,psym=-8
      if keyword_set(flmin) then oplot, lnew[0:nband-1]-dl, msig[0:nband-1,3], color=colvec[1],thick=2,psym=-8
  endelse
  oplot,indgen(llmax),0*indgen(llmax)+one_percent,line=2,thick=1,color=0

  if keyword_set(flmin) and keyword_set(flmax) then begin
     if max( msig[0:nband-1,0]) gt 50 then begin
        legend,legtxt,color=colvec[[0,1]],textcolor=[0,0],/bottom,/right,box=0,$
               line=[0,0],pspacing=1,thick=4,charsize=lcharsize
     endif else begin
        legend,legtxt,color=colvec[[0,1]],textcolor=[0,0],/top,/right,box=0,$
               line=[0,0],pspacing=1,thick=4,charsize=lcharsize
     endelse
  endif

  ax=['xttl','yttl']
  ay=['$\ell_{\rm '+lminlmax+'}$','Significances of dispersion $[\%]$']

  if keyword_set(fname) then ps_end_planck, /png,feps=fcase,/latex,/cl,axisx=ax,axisy=ay,xsize=xsize,ysize=ysize
endif

end
