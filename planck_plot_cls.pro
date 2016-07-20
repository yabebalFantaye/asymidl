pro planck_plot_cls,dxy,xy_sim,mxy_sim,dxy_sim,fname=fname,psend=psend,scale=scale,ay=ay,ax=ax,$
                 xra=xra, yra=yra, nx=nx,nmx=nmx,nmy=nmy,ny=ny,colvec=colvec,legtxt=legtxt,$
                    lvec=lvec,rmll=rmll,multll=multll, txtshift=txtshift,siglabel=siglabel,sig_amp=sig_amp,nolabel=nolabel,$
                    _extra=extra

  if n_params() lt 1 then begin 
     print, 'usage:'
     print, 'planck_plot_cls,dxy,xy_sim,mxy_sim,dxy_sim,fname=fname,psend=psend,scale=scale,$'
     print, '                 xra=xra, yra=yra, nx=nx,nmx=nmx,nmy=nmy,ny=ny,colvec=colvec,legtxt=legtxt,$'
     print, '                 lvec=lvec,_extra=extra'
     return
  endif



   zzz=dxy
   if size(dxy,/type) eq 8 then zzz=dxy[0].xy

   nsimd = n_elements(zzz[0,*])
   nsim=nsimd

   
   if n_params() gt 1 then    begin
      nsim = n_elements(xy_sim[0,*])
      mxy_sim=mean(xy_sim,dimension=2)
      dxy_sim=stddev(xy_sim,dimension=2)
   endif


   nband = n_elements(zzz[*,0])

   lnew = indgen(nband)+1
   if keyword_set(lvec) then lnew=lvec
   ll=lnew*(lnew+1)/(2d0*!pi)

   if not keyword_set(scale) then scale = lnew*0 + 1.
   if keyword_set(rmll) then scale=scale/ll
   if keyword_set(multll) then scale=scale*ll


   if not keyword_set(ny) then ny = 5
   if not keyword_set(nx) then nx = 3
   if not keyword_set(nmx) then nmx=5
   if not keyword_set(nmy) then nmy=5

   if not keyword_set(xra) then xra = [min(lnew),max(lnew)]
   if not keyword_set(yra) then yra = [floor(min(xy_sim)) ,(ny+1)*(ceil(max(zzz))-floor(min(zzz)))/ny]
   


   planck_ct, cvec


   planck_plot,lnew,zzz[*,0]*scale,color=0,fname=fname, $
               xra=xra, yra=yra, thick=2,xt='xttl', yt='yttl', $ ;, XTTL_DY=xtdy
               /noeras,/nodata, nx=nx,nmx=nmx,nmy=nmy,ny=ny,xs=xsize,ys=ysize,_extra=extra

   if not keyword_set(colvec) then colvec=cvec

   LS_circle,sz=0.4

   ;;simulation case
   if n_params() gt 1 then begin
      for kk=0,nsim-1 do begin
         oplot,lnew, xy_sim[*,kk]*scale,color=20,symsize=0.2,psym=8
      endfor
      if n_params() gt 3 then begin 
         sigma = dxy_sim
         oband, lnew,  (mxy_sim-sigma)*scale, (mxy_sim+sigma)*scale,color=19,border=colvec[0],symsize=0.2,psym=8
      endif
   endif

   ;;data case
   for ij = 0, nsimd-1 do begin
      if size(dxy,/type) eq 8 then begin
         ddipamp = dxy[0].xy
         for kk=0,n_elements(dxy)-1 do begin
            oplot, lnew, dxy[kk].xy[*,ij]*scale, color=colvec[kk],thick=2,psym=-8
         endfor
      endif else  begin
         ddipamp=dxy
         oplot, lnew, dxy[*,ij]*scale, color=colvec[0],thick=2,psym=-8
      endelse
   endfor

   ;;oplot,indgen(max(lnew)+10),
   ;;indgen(max(lnew)+10)*0,color=0,line=0,thick=1

   ;;compute significance based on amplitude

   if keyword_set(siglabel) then begin

      if not keyword_set(txtshift) then txtshift=0

      dipamp=xy_sim
      sig_amp=dblarr(nband)
      for ijk=0,nband-1 do begin
      
         x = dipamp[ijk,*]*scale[ijk] & y = ddipamp[ijk]*scale[ijk]         
         ;print, max(x[0,*]), y[0]

         sig_amp[ijk]=total( x[0,*] gt y[0] )  

         print, 'sig: ',sig_amp[ijk]

         x = reform(x)
         y=reform(y)

         ;;print, 'sig yval: ',(max(x)+(max(x)-mean(x))/10.)

         if not keyword_set(nolabel) then   cgtext,lnew[ijk]+txtshift,mean(x), strn(round(sig_amp[ijk])),col=0,charsize=0.7,charthick=0.5
      endfor

   endif


   if keyword_set(legtxt) then begin

      legend,legtxt,color=colvec[0:n_elements(legtxt)-1],textcolor=colvec[0:n_elements(legtxt)-1],/top,/right,box=0,$
             pspacing=1,thick=4,charsize=lcharsize ;, line=0*indgen(n_elements(legtxt))
   endif



   if keyword_set(ay)  and not keyword_set(ax) then axisx=['xttl','yttl']

   if keyword_set(fname) then  ps_end_planck, /png,feps=fname,/latex,/cl ,xsize=xsize,ysize=ysize,axisy=ay,axisx=ax



end
