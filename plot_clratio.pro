pro plot_clratio,dxy,xy_sim,mxy_sim,dxy_sim,fname=fname,psend=psend,$
                 xra=xra, yra=yra, nx=nx,nmx=nmx,nmy=nmy,ny=ny,colvec=colvec,legtxt=legtxt



   if not keyword_set(xra) then xra = [0,1500]
   if not keyword_set(yra) then yra = [-10,15]
   if not keyword_set(ny) then ny = 5
   if not keyword_set(nx) then nx = 3
   if not keyword_set(nmx) then nmx=5
   if not keyword_set(nmy) then nmy=5
   

   nsim = n_elements(xy_sim[0,*])
   nband = n_elements(xy_sim[*,0])

   lnew = indgen(nband)*100l+50l

   planck_ct, cvec

   planck_plot,lnew,mxy_sim[*,0]*100,color=0,fname=fname, $
               xra=xra, yra=yra, thick=2,xt='xttl', yt='yttl2', $ ;, XTTL_DY=xtdy
               /noeras,/nodata, nx=nx,nmx=nmx,nmy=nmy,ny=ny,xs=xsize,ys=ysize

   if not keyword_set(colvec) then colvec=cvec

   LS_circle,sz=0.4

   

   for kk=0,nsim-1 do begin
      oplot,lnew, xy_sim[*,kk]*100,color=20
   endfor
   sigma = dxy_sim
   oband, lnew,  (mxy_sim-sigma)*100, (mxy_sim+sigma)*100,color=19,border=colvec[0]
   
   if size(dxy,/type) eq 8 then begin
      for kk=0,n_elements(dxy)-1 do begin
         oplot, lnew, dxy[kk].xy*100, color=colvec[kk],thick=2,psym=-8
      endfor
   endif else    oplot, lnew, dxy*100, color=colvec[0],thick=2,psym=-8

   oplot,indgen(2500), indgen(2500)*0,color=0,line=0,thick=1


   if keyword_set(legtxt) then begin

      legend,legtxt,color=colvec,textcolor=colvec[0:n_elements(legtxt)-1],/top,/right,box=0,$
             line=0*indgen(n_elements(legtxt)),pspacing=1,thick=4,charsize=lcharsize
   endif

   if keyword_set(fname) then  ps_end_planck, /png,feps=fname,/latex,/cl ,xsize=xsize,ysize=ysize



end
