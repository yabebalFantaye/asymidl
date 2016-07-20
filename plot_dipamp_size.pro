pro plot_dipamp_size,ddipamp,dipamp,fname=fname, legtxt=legtxt,issim=issim,val=val,ay=ay,sig_amp=sig_amp,$
  xra=xra, yra=yra, nx=nx,nmx=nmx,nmy=nmy,ny=ny,_extra=extra,dshift=dshift,txtshift=txtshift,xvec=xvec,nolabel=nolabel,DECMODX=DECMODX, DECMODY=DECMODY,format=format,xdipamp=xdipamp

if not keyword_set(legtxt) then legtxt=['Data ', 'Isotropic']
if not keyword_set(xra) then xra = [0,1500]
if not keyword_set(ny) then ny = 4
if not keyword_set(nx) then nx = 4
if not keyword_set(nmx) then nmx=1
if not keyword_set(nmy) then nmy=1

if not keyword_set(dshift) then dshift=0
if not keyword_set(txtshift) then txtshift=13

;;xvec is the x value
;;val is value at the xvec 

llmax=xra[1]

ymax = (ny-1)*float(sigfig(max(dipamp)/(ny-1),2))

if not keyword_set(yra) then begin
   yra = [0,ymax]
endif


if keyword_set(fname) then print, 'saving file to '+fname

nband = n_elements(dipamp[0,*,0])
nsim = n_elements(dipamp[0,0,*])
nradi = n_elements(ddipamp[*,0,0])



psym2=-1.5
thick2=1.3
symsize2=0.3
if nradi gt 1 then begin
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

planck_plot,lnew[0:nband-1],ddipamp[0,0:nband-1,0],color=0,fname=fname, $
            xra=xra, yra=yra, thick=2,xt='xttl',yt='yttl', DECMODX=DECMODX, DECMODY=DECMODY,format=format,$ 
            nx=nx,nmx=nmx,nmy=nmy,ny=ny,xs=xsize,ys=ysize,/noerase,/nodata,psym=-2,symsize=0, _extra=extra ,XTickV=lnew, XTickname=name


;; for ijk=0,nsim-1 do begin
;;    oplot,lnew[0:nband-1],dipamp[0:nband-1,ijk],color=20,thick=0.5,psym=2,symsize=0.2
;; endfor


;;compute significance based on amplitude
sig_amp=dblarr(nradi,nband)
sym_size=sig_amp
for tt=0,nradi-1 do begin
   for ijk=0,nband-1 do begin
      x = reform(dipamp[tt,ijk,*]) & y = reform(ddipamp[tt,ijk])
    
      sig_amp[tt,ijk]=total(  x  gt y )  
      sym_size[tt,ijk] = round(nsim*min([sig_amp[tt,ijk]/double(nsim),0.01])) ;;p-value converted 0-10
   endfor
   
endfor

; Add a colorbar
ncolors=10
loadct,39, NColors=ncolors
loadct,39, NColors=ncolors,rgb_table=rgb

planck_ct
print, sig_amp
print, sym_size

;;plot anisotropic simulation dipole amplitudes
mysim=['o','s','h','tu','S','X','+']
max_tt = nradi-1 ;min([nradi-1,n_elements(mysim)])
for tt=0,max_tt do begin

   LS_circle,sz=1+float(tt)/float(max_tt)
   oplot,lnew[0:nband-1],reform(ddipamp[tt, 0:nband-1]),psym=-8, color=11+tt,thick=2.-1.5*tt/max_tt

   ;; p=plot(lnew[0:nband-1],reform(ddipamp[tt, 0:nband-1]),symbol='o',rgb_table=rgb, sym_thick=2.-1.5*tt/max_tt,$
   ;;             vert_colors=reform(sym_size[tt,0:nband-1]),/current)
endfor

;c = colorbar(taper=3,tickname=strn(indgen(ncolors)),rgb_table=39, ORIENTATION=0, /normal, POSITION=[0.2,0.1,0.8,0.3], TITLE='p-value')


;; if keyword_set(xdipamp) then begin
;;    nsim3=n_elements(xdipamp[0,*])
;;    for ijk=0,nsim3-1 do begin
;;       oplot,lnew[0:nband-1]+2*dshift,xdipamp[0:nband-1,ijk],color=13,thick=0.5,psym=2,symsize=0.2
;;    endfor
;;    legend,legtxt,color=[11,13,20],textcolor=[0,0,0],/top,/left,box=0,$
;;           line=[0,0,0],pspacing=1,thick=4,charsize=lcharsize
;; endif else begin
;;    legend,legtxt,color=[11,20],textcolor=[0,0],/top,/left,box=0,$
;;           line=[0,0],pspacing=1,thick=4,charsize=lcharsize
;; endelse

if not keyword_set(ax) then ax=['xttl','yttl']
if not keyword_set(ay) then ay=['$\ell$','Dipole Amplitude']
if keyword_set(fname) then ps_end_planck, /png,feps=fname,/latex,/cl,axisx=ax,axisy=ay,xsize=xsize,ysize=ysize,/nops

end
