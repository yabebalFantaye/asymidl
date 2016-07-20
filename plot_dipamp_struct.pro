pro plot_dipamp,ddipamp,dipamp,fname=fname, legtxt=legtxt,issim=issim,val=val,ay=ay,$
  xra=xra, yra=yra, nx=nx,nmx=nmx,nmy=nmy,ny=ny,_extra=extra,dshift=dshift,txtshift=txtshift,xvec=xvec,nolabel=nolabel,DECMODX=DECMODX, DECMODY=DECMODY,format=format,xdipamp=xdipamp

if not keyword_set(legtxt) then legtxt=['Data ', 'Isotropic']
if not keyword_set(xra) then xra = [0,1500]
if not keyword_set(ny) then ny = 4
if not keyword_set(nx) then nx = 4
if not keyword_set(nmx) then nmx=1
if not keyword_set(nmy) then nmy=1

if not keyword_set(dshift) then dshift=10
if not keyword_set(txtshift) then txtshift=13

;;xvec is the x value
;;val is value at the xvec 

llmax=xra[1]


if not keyword_set(yra) then begin
   yra = [0,ceil(1.2*max(dipamp))]
   yra = [floor(-yra[1]/(ny-1.)),(ny-1)*ceil(yra[1]/(ny-1.))]
endif


if keyword_set(fname) then print, 'saving file to '+fname

nband = n_elements(dipamp[*,0])
nsim = n_elements(dipamp[0,*])
nsim2 = n_elements(ddipamp[0,*])



psym2=-2
thick2=1.5
symsize2=0.5
if nsim2 gt 1 then begin
   psym2=2
   thick2=0.5
   symsize2=0.2
endif
nbin2bin=6l
dl=100
dl_cent=50
if nband eq 93l then begin
   nbin2bin=1l
   dl=16
   dl_cent=8
endif
if (nband eq 1) then nbin2bin=floor(llmax/(16l*nband))
lnew = indgen(nband)*dl+dl_cent


if keyword_set(val) then lnew=val


;; name=replicate('xyz',nband)
;; for i=0,nband-1 do begin
;;    name[i]=strn(val[i])
;; endfor

;; print, name

planck_plot,lnew[0:nband-1],ddipamp[0:nband-1,0],color=0,fname=fname, $
            xra=xra, yra=yra, thick=2,xt='xttl',yt='yttl', DECMODX=DECMODX, DECMODY=DECMODY,format=format,$ 
            nx=nx,nmx=nmx,nmy=nmy,ny=ny,xs=xsize,ys=ysize,/noerase,/nodata,psym=-2,symsize=0, _extra=extra ;,XTickV=lnew, XTickname=name


for ijk=0,nsim-1 do begin
   oplot,lnew[0:nband-1],dipamp[0:nband-1,ijk],color=20,thick=0.5,psym=2,symsize=0.2
endfor


;;compute significance based on amplitude
sig_amp=dblarr(nband)
for ijk=0,nband-1 do begin

   x = dipamp[ijk,*] & y = ddipamp[ijk]
   if keyword_set(issim) then  x = dipamp[ijk,*] & y = ddipamp[ijk,*]


   sig_amp[ijk]=total( x[0,*] gt y[0] )  
   if keyword_set(issim) then sig_amp[ijk]=total(y[0,*] lt max(x[0,*]))

  ;; print,'i, sig_amp: ', ijk, sig_amp[ijk]
   if not keyword_set(nolabel) then   cgtext,lnew[ijk]+txtshift,mean(x), strn(round(sig_amp[ijk])),col=0,charsize=0.5,charthick=0.5
endfor

;;plot anisotropic simulation dipole amplitudes
for ijk=0,nsim2-1 do begin
   oplot,lnew[0:nband-1]+dshift,ddipamp[0:nband-1,ijk],color=11,thick=thick2,psym=psym2,symsize=symsize2
endfor
;if nsim2 gt 10 then oplot,lnew[0:nband-1],min(ddipamp[0:nband-1,*],dim=2),color=11,thick=2,psym=-2,symsize=0.5


if keyword_set(xdipamp) then begin
   nsim3=n_elements(xdipamp[0,*])
   for ijk=0,nsim3-1 do begin
      oplot,lnew[0:nband-1]+2*dshift,xdipamp[0:nband-1,ijk],color=12,thick=thick,psym=psym,symsize=symsize
   endfor
   legend,legtxt,color=[11,12,20],textcolor=[0,0,0],/top,/right,box=0,$
          line=[0,0,0],pspacing=1,thick=4,charsize=lcharsize
endif else begin
   legend,legtxt,color=[11,20],textcolor=[0,0],/top,/right,box=0,$
          line=[0,0],pspacing=1,thick=4,charsize=lcharsize
endelse

ax=['xttl','yttl']
if not keyword_set(ay) then ay=['$\ell$','Dipole Amplitude']
if keyword_set(fname) then ps_end_planck, /png,feps=fname,/latex,/cl,axisx=ax,axisy=ay,xsize=xsize,ysize=ysize

end
