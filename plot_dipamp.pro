pro plot_dipamp,ddipamp,dipamp,fname=fname, legtxt=legtxt,issim=issim,val=val,ay=ay,sig_amp=sig_amp,$
                xra=xra, yra=yra, nx=nx,nmx=nmx,nmy=nmy,ny=ny,_extra=extra,dshift=dshift,txtshift=txtshift,$
                nolabel=nolabel,DECMODX=DECMODX, DECMODY=DECMODY,format=format,xdipamp=xdipamp,$
                xtickname=xtickname,xtickv=xtickv,xvec=xvec

if not keyword_set(legtxt) then legtxt=['Data ', 'Isotropic']

if not keyword_set(ny) then ny = 3
if not keyword_set(nx) then nx = 4
if not keyword_set(nmx) then nmx=1
if not keyword_set(nmy) then nmy=1

if not keyword_set(dshift) then dshift=0
if not keyword_set(txtshift) then txtshift=13

;;xvec is the x value
;;val is value at the xvec 


ymax = 1.1*max([max(ddipamp),max(dipamp)])
ymax = (ny-1)*float(sigfig(ymax/(ny-1),2))

if not keyword_set(yra) then begin
   yra = [0,ymax]
endif


if keyword_set(fname) then print, 'saving file to '+fname

nband = n_elements(dipamp[*,0])
nsim = n_elements(dipamp[0,*])
nsim2 = n_elements(ddipamp[0,*])



psym2=-1.5
thick2=1.3
symsize2=0.3
if nsim2 gt 1 then begin
   psym2=2
   thick2=0.5
   symsize2=0.2
endif


lnew = indgen(nband)+1
xra = [min(lnew)-1, max(lnew)+1]
xtickv=[lnew[0]-1,lnew,lnew[nband-1]+1]

;;if not keyword_set(xra) then 
;;if not keyword_set(xtickv) then 

if not keyword_set(val) then val=lnew


if not keyword_set(xtickname) then begin
   xtickname=replicate(' ',nband+2)
   for i=1,nband-1 do begin
      xtickname[i]=strn(round(val[i]))
   endfor
endif
nx=n_elements(xtickv)-1

print, 'xra: ',xra
print, 'xgrid plot: ',lnew
print, 'xtickv: ',xtickv
print, 'xtickname: ',xtickname

planck_plot,lnew[0:nband-1],ddipamp[0:nband-1,0],color=0,fname=fname, $
            xra=xra, yra=yra, thick=2,xt='xttl',yt='yttl', DECMODX=DECMODX, DECMODY=DECMODY,format=format,$ 
            nx=nx,nmx=nmx,nmy=nmy,ny=ny,xs=xsize,ys=ysize,/noerase,/nodata,psym=-2,symsize=0,$
            XTickV=xtickv, XTickname=xtickname, _extra=extra 


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
      oplot,lnew[0:nband-1]+2*dshift,xdipamp[0:nband-1,ijk],color=13,thick=0.5,psym=2,symsize=0.2
   endfor
   legend,legtxt,color=[11,13,20],textcolor=[0,0,0],/top,/left,box=0,$
          line=[0,0,0],pspacing=1,thick=4,charsize=lcharsize
endif else begin
   legend,legtxt,color=[11,20],textcolor=[0,0],/top,/left,box=0,$
          line=[0,0],pspacing=1,thick=4,charsize=lcharsize
endelse

if not keyword_set(ax) then ax=['xttl','yttl']
if not keyword_set(ay) then ay=['$\ell$','Dipole Amplitude']
if keyword_set(fname) then ps_end_planck, /png,feps=fname,/latex,/cl,axisx=ax,axisy=ay,xsize=xsize,ysize=ysize,/nops

end
