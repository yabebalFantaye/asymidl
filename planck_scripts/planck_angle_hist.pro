pro planck_angle_hist,mm,fname=fname,ny=ny,nmy=nmy,yra=yra,nx=nx,$
                      nmx=nmx,xra=xra,nbins=nbins,ax=ax,ay=ay,legtxt=legtxt,$
                      help=help,stat=stat,pos=pos
;;Plot histogram of the dipole angle distributions or amplitudes

if keyword_set(help) then begin
   print, 'Usage: '
   print, 'planck_angle_hist,mm,fname=fname,ny=ny,nmy=nmy,yra=yra,nx=nx,$'
   print, '                      nxm=nxm,xra=xra,nbins=nbins,ax=ax,ay=ay,legtxt=legtxt,$'
   print, '                      help=help'
endif

if keyword_set(fname) then begin
   print, 'Histogram plot will be save in: ',fname
   ps_start_planck,file=fname,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize ;,/large
endif

planck_ct,cvec  ;;load Planck style colors for CR, NILC,SEVEM, SMICA, sim

vcol=cvec
tcvec = [0,0,0,0,0,0]
lcharsize=1


if not keyword_set(pos) then pos = [0.2,0.9]

if not keyword_set(xra) then xra = [min(mm),max(mm)]
if not keyword_set(nx) then nx=5
if not keyword_set(nmx) then nmx = 0
mmin = xra[0]
mmax = xra[1]

resolve_routine,'hfi_plot',/either, /compile_full_file

nsim = n_elements(mm[*,0])-1
nmask = n_elements(mm[0,*])

mnangle = mm[1:nsim,0]
mndat = mm[0,0]

if not keyword_set(nbins) then nbins=30

stat=dblarr(3,nmask)
for iii=0,nmask-1 do begin            
   mnangle = mm[1:nsim,iii]
   mndat = mm[0,iii]
   medang=median(mnangle)
   stdang=stdev(mnangle)
   meanang=mean(mnangle)

   stat[*,iii]=[meanang,stdang,medang]
endfor

hist = Histogram(mnangle, min=mmin, max=mmax,nbins=nbins,location=xbins)

if not keyword_set(yra) then yra = [0,max([10*(max(hist)/10),10])]
if not keyword_set(ny) then ny=4
if not keyword_set(nmy) then nmy = 1


print,'hist yra = ',yra

ytickv = intarr(ny+1)
ytickname = strarr(ny+1)
for jj=0,ny do begin
   yval = yra[0] + jj*(yra[1]-yra[0])/ny
   ytickv[jj] = yval
   ytickname[jj] = strn(yval)
endfor


plot,xbins,hist,color=0,xr=xra,yr=yra,/xs,/ys, BACKGROUND=255,charsize=0.8,thick=2,$
         YTITLE='yttl', XTITLE='xttl', /noeras,xticks=nx,XMINOR=nmx,yminor=nmy,yticks=ny,$
     /nodata ,ytickv=ytickv,ytickname=ytickname




 for iii=0,nmask-1 do begin            
   mnangle = mm[1:nsim,iii]
   mndat = mm[0,iii]

   hist = Histogram(mnangle, min=mmin, max=mmax,nbins=nbins,location=xbins)

   oplot, xbins, hist, PSYM = 10, col=vcol[iii],thick=1.5
   oplot,mndat+0.*indgen(yra[1]), indgen(max(yra[1])+1),thick=3,col=vcol[iii]
endfor

if keyword_set(legtxt) then legend,legtxt[0:(nmask-1)],color=cvec[0:(nmask-1)],textcolor=tcvec[0:(nmask-1)],pos=pos,/norm,box=0,$
       line=tcvec[0:(nmask-1)],pspacing=1,thick=4,charsize=lcharsize

if not keyword_set(ax) then ax=['xttl','yttl']
if not keyword_set(ay) then ay = ['Dipole dispersion angles (deg)','Frequency']
if keyword_set(fname) then ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=ax,axisy=ay

end
