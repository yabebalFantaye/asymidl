pro planck_hist,mm,fname=fname,ny=ny,nmy=nmy,yra=yra,nx=nx,rot=rot,$
                nmx=nmx,xra=xra,nbins=nbins,ax=ax,ay=ay,legtxt=legtxt,$
                help=help,stat=stat,pos=pos,deftxt=deftxt,defpos=defpos,$
                charsize=charsize,lcharsize=lcharsize,noplot=noplot, norm=norm
;;Plot histogram of the dipole angle distributions or amplitudes

  if keyword_set(help) then begin
     print, 'Usage: '
     print, 'planck_angle_hist,mm,fname=fname,ny=ny,nmy=nmy,yra=yra,nx=nx,$'
     print, '                      nxm=nxm,xra=xra,nbins=nbins,ax=ax,ay=ay,legtxt=legtxt,$'
     print, '                      help=help'
  endif




  if not keyword_set(lcharsize) then lcharsize=1
  if not keyword_set(charsize) then charsize=1


  if not keyword_set(pos) then pos = [0.7,0.9]


  nnx=4

  nsim = n_elements(mm[*,0])
  mdata=reform(mm[0,*])
  mdatamin = min(mdata)
  mdatamax = max(mdata)

  if not keyword_set(xra) then begin
     boundsim = transpose(mm[1:nsim-1,*])
     pvalvec=[0.003,0.997]
     simconflevel,boundsim,y1,y2,mmat=mmat,pvalue=pvalvec
     
     msimmin = min(y1)
     msimmax = max(y2)
     
     xmin = min([mdatamin,msimmin])
     xmax = max([mdatamax,msimmax])
     if abs(xmax) gt 1 or abs(xmin) gt 1 then begin
        xmin=floor(xmin)
        xmax=ceil(xmax)
     endif
     dx = (xmax-xmin)/(nnx+1)
     xmax = xmin+(nnx+1)*dx
     xra = [xmin,xmax]
  endif

  if not keyword_set(nx) then nx=nnx+1
  if not keyword_set(nmx) then nmx = 1
  mmin = xra[0]
  mmax = xra[1]

  resolve_routine,'hfi_plot',/either, /compile_full_file

;;here just make nsim to be only for sims
  nsim = n_elements(mm[*,0])-1
  nmask = n_elements(mm[0,*])

  mnangle = mm[1:nsim,0]
  mndat = mm[0,0]

  if not keyword_set(nbins) then nbins=30

  stat=dblarr(nmask,5)
  for iii=0,nmask-1 do begin            
     mnangle = mm[1:nsim,iii]
     mndat = mm[0,iii]

     xfrac=total(mndat lt mnangle) ;*100.0/double(nsim-1)
     print, 'iii, data, minsim, maxsim, pvalue',iii,mndat,min(mnangle),max(mnangle), reform(xfrac)
     medang=median(mnangle)
     stdang=stdev(mnangle)
     meanang=mean(mnangle)

     stat[iii,*]=[xfrac,mndat,meanang,stdang,xfrac*100d0/double(nsim)]
  endfor

  hist = Histogram(mnangle, min=mmin, max=mmax,nbins=nbins,location=xbins)

  nny=2
  ymax = max(hist)


  if keyword_set(norm) then begin
     ymax=1
     hist = float(hist)/max(float(hist))
  endif


  if not keyword_set(yra) then yra = [0, ymax]
  if not keyword_set(ny) then ny=nny+1
  if not keyword_set(nmy) then nmy = 1


  print,'**** hist yra = ',yra

  ytickv = intarr(ny+1)
  ytickname = strarr(ny+1)
  for jj=0,ny do begin
     yval = yra[0] + jj*(yra[1]-yra[0])/ny
     ytickv[jj] = yval
     ytickname[jj] = strn(yval)
  endfor
  format = '(g0)'  

  ;; yra = smart_yra(yra[0], yra[1],ytickv=ytickv,ytickname=ytickname,format=format,/chatty)
  ;; ny = n_elements(ytickv)

  if keyword_set(rot) then yformat = '(A1)' 
  if not keyword_set(rot) then yformat = '(g0)' 

  print,'**** hist xra = ',xra

  xtickv = intarr(nx+1)
  xtickname = strarr(nx+1)
  for jj=0,nx do begin
     xval = xra[0] + jj*(xra[1]-xra[0])/nx
     xtickv[jj] = xval
     xtickname[jj] = strn(xval)
  endfor

  xtt='xttl'
  ytt='yttl'
  if not keyword_set(ax) then ax=[xtt,ytt]

  ;;===============================

  if not keyword_set(noplot) then begin

     if keyword_set(fname) then begin
        print, 'Histogram plot will be save in: ',fname
        ps_start_planck,file=fname,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize ;,/large
     endif


     pstyle_plot,xbins,hist,color=0,xr=xra,yr=yra,BACKGROUND=255,charsize=charsize,/xs,/ys,$
          YTITLE=ax[1]+'!C', XTITLE=ax[0],xticks=nx,XMINOR=5,yticks=ny,$
          /nodata, xcharsize=1.2, ycharsize=1,xmar=[10,4],ymar=[4,2]
;;,xtickv=xtickv,xtickname=xtickname, 
     

     planck_ct,cvec  ;;load Planck style colors for CR, NILC,SEVEM, SMICA, sim
     
     vcol=cvec
     tcvec = [0,0,0,0,0,0]
     
     for iii=0,nmask-1 do begin            
        mnangle = mm[1:nsim,iii]
        mndat = mm[0,iii]
       
        hist = Histogram(mnangle, min=mmin, max=mmax,nbins=nbins,location=xbins)

        print, 'min(hist), max(hist)',min(hist), max(hist)
        print, 'hist before ',hist

        hist = float(hist)/max(float(hist))

        print, 'min(hist), max(hist)',min(hist), max(hist)
        print, 'after before ',hist

        oplot, xbins, hist, PSYM = 10, col=vcol[iii],thick=1.5
        oplot,mndat+0.*indgen(max([100,hist])+1), indgen(max([100,hist])+1),thick=3,col=vcol[iii]
     endfor

     if keyword_set(legtxt) then legend,legtxt[0:(nmask-1)],color=cvec[0:(nmask-1)],textcolor=tcvec[0:(nmask-1)],pos=pos,/norm,box=0,$
                                        line=tcvec[0:(nmask-1)],pspacing=1,thick=4,charsize=lcharsize

     if keyword_set(deftxt) then ndef=n_elements(deftxt)
     if keyword_set(deftxt) then legend,deftxt[0:(ndef-1)],color=0*cvec[0:(ndef-1)],textcolor=0*tcvec[0:(ndef-1)],pos=defpos,/norm,box=0,$
                                        pspacing=1,thick=4,charsize=lcharsize


     if not keyword_set(ay) then ay = ['Dipole dispersion angles (deg)','{\rm Probability}']
     if keyword_set(fname) then ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=ax,axisy=ay,charsize=replicate(1.5,n_elements(ay))


  endif

end
