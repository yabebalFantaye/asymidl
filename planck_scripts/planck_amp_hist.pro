pro planck_amp_hist,xy_sim, d, llmin=llmin,llmax=llmax,fname=fname,ny=ny,nym=nym,$
                    yra=yra,nx=nx,nxm=nxm,xra=xra,nbins=nbins,ax=ax,ay=ay,legtxt=legtxt,$
                    help=help

if keyword_set(help) then begin
   print, 'usage: '
   print, 'planck_amp_hist,xy_sim, d, llmin=llmin,llmax=llmax,fname=fname,ny=ny,nym=nym,$'
   print,'             yra=yra,nx=nx,nxm=nxm,xra=xra,nbins=nbins,ax=ax,ay=ay,legtxt=legtxt,$'
   print, '            help=help'
   return
endif


if keyword_set(fname) then  begin
   print, 'Clratio figure will be saved in: ',fname
   ps_start_planck,file=fname,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize ;,/large
endif

planck_ct,cvec  ;;load Planck style colors for CR, NILC,SEVEM, SMICA, sim

nband = n_elements(xy_sim[*,0,0])
nsim = n_elements(xy_sim[0,*,0])
nmask = n_elements(xy_sim[0,0,*])

if not keyword_set(nbins) then nbins=30


bmin=0
lastb = nband-1

if keyword_set(llmin) then bmin=min([0,(llmin-100)/100])
if keyword_set(llmax) then lastb = llmax/100 - 1




if not keyword_set(xra) then xra = [-2,4]
if not keyword_set(nx) then nx=4
if not keyword_set(nxm) then nxm = 0
mmin = xra[0]
mmax = xra[1]

print, lastb
zzz = reform(mean(xy_sim[0:lastb,*,0]*100,dimension=1,/double))
hist = Histogram(zzz, min=mmin, max=mmax,location=xbins)

if not keyword_set(yra) then yra = [0,max([10*(max(hist)/10),10])]
if not keyword_set(ny) then ny=4
if not keyword_set(nym) then nym = 1

ytickv = intarr(ny+1)
ytickname = strarr(ny+1)
for jj=0,ny do begin
   yval = yra[0] + jj*(yra[1]-yra[0])/ny
   ytickv[jj] = yval
   ytickname[jj] = strn(yval)
endfor





hfi_plot,xbins,hist,color=0,$
         xr=xra,yr=yra,/xs,/ys, BACKGROUND=255,charsize=0.8,thick=2,$
         YTITLE='yttl', XTITLE='xttl', XMARGIN=XMAR, YMARGIN=YMAR,$ ; YTTL_DX = ytdx,
         /noeras,xticks=nx,XMINOR=nxm,yminor=nym,yticks=ny,X_DX=0,/nodata ,ytickv=ytickv,ytickname=ytickname



for ii=0,nmask-1 do begin
   zzz = reform(mean(xy_sim[bmin:lastb,*,ii]*100,dimension=1,/double))
   hist = Histogram(zzz, min=-3, max=3,nbins=nbins, locations=xbins)


   sig_ind = (where(mean(d[bmin:lastb,ii]*100) le zzz))
   sig_num=0
   if total(sig_ind) gt 0 then begin
      sig_num = n_elements(sig_ind) 
      print, 'simulations with mean cl/cl gt data: ',zzz[sig_ind]
   endif
   
   print, 'mean, sig/'+strn(nsim)+', sig(%): ', mean(d[bmin:lastb,ii]*100),sig_num,sig_num*100d0/(nsim*1d0)

   oplot,xbins,hist,color=cvec[ii],thick=1.5,psym=10
   oplot,mean(d[bmin:lastb,ii]*100)+0.*indgen(max(hist)+200), indgen(max(hist)+200),thick=3,color=cvec[ii] ;,line=ii
endfor

if keyword_set(legtxt) then legend,legtxt[0:(nmask-1)],color=cvec[0:(nmask-1)],textcolor=tcvec[0:(nmask-1)],/top,/right,box=0,$
       line=tcvec[0:(nmask-1)],pspacing=1,thick=4,charsize=lcharsize

if not keyword_set(ax) then ax=['xttl','yttl']
if not keyword_set(ay) then ay = ['$ \langle \Delta C_\ell/C_\ell \rangle [\%]$','Frequency']
if keyword_set(fname) then ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=ax,axisy=ay


end
