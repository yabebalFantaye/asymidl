pro plot_mollview,map,minval,maxval,large=large,medium=medium,fopen=fopen,show_dipdir=show_dipdir,mtheta=mtheta,mphi=mphi,$
                  file_ps=file_ps,ctitle=ctitle,title=title,cl=cl,ax=ax,ay=ay,fclose=fclose,help=help,nopng=nopng,$
                  midval=midval,tkey=tkey,colt=colt,hist=hist,ellcol=ellcol,ellsize=ellsize,method_vec=method_vec,$
                  fthph_out=fthph_out,thph_out=thph_out,ecp=ecp,lowell=lowell,val=val,txpx=txpx,ftxpx=ftxpx,fftxpx=fftxpx,$
                  stxpx=stxpx,xsc=xsc,ysc=ysc,txrad=txrad,ddir=ddir,bwhite=bwhite,fwhm=fwhm,lmax_smooth=lmax_smooth,mask=mask,$
                  _extra=extra,mollview=mollview,fmtbar=fmtbar, ncolors=ncolors,grls=grls,nocgbar=nocgbar

;;stxpx is the [color,charsize, chartick] property of txpx

;;thph_out (fthph_out) and txpx (ftxpx) serve the same purpose, draw
;;arbitrary point in the map with fthph name. txpx can be an array


  if keyword_set(help) or n_params() eq 0 then begin
     print, 'USAGE: '
     print, 'plot_mollview,map,minval,maxval,large=large,medium=medium,fopen=fopen,show_dipdir=show_dipdir,mtheta=mtheta,mphi=mphi,'
     print, '    file_ps=file_ps,ctitle=ctitle,title=title,cl=cl,ax=ax,ay=ay,fclose=fclose,help=help,$'
     print, '       midval=midval,tkey=tkey,colt=colt,hist=hist,ellcol=ellcol,ellsize=ellsize,method_vec=method_vec,$'
     print, '             fthph_out=fthph_out,thph_out=thph_out,ecp=ecp,lowell=lowell'
  endif

  png=1
  if keyword_set(nopng) then png=0

;;plot dipole directions 
  if not keyword_set(ctitle) then ctitle=' '
  if not keyword_set(title) then title=' '
  if not keyword_set(midval) then midval = 0
  if not keyword_set(colt) then colt = 41
  if not keyword_set(ellcol) then ellcol=0
  if not keyword_set(ellsize) then ellsize=0.4

  nobar=1
  if keyword_set(nocgbar) then nobar=0

  if not keyword_set(grls) then begin
     grls=1
     glsize=0.5
  endif
  if grls lt 0 then begin
     grls=0
     glsize=0
  endif

  white=1
  if keyword_set(bwhite) then white=0

  if midval lt 1 then midval=0

  SZ = 88d
  scale_sz=1
  if keyword_set(medium) then begin
     SZ=120d
     scale_sz=1.5
  endif
  if keyword_set(large) then begin
     SZ=180d
     scale_sz=2
  endif

  ysize=SZ/10d/1.6d
  xsize=SZ/10d

  ;;check if smoothing is set
  if keyword_set(fwhm) then begin
     if not keyword_set(lmax_smooth) then lmax_smooth = 2*npix2nside(n_elements(map))
     print, '*****************************'
     print, 'plot_mollview: map will be smoothed by fwhm arcmin and lmax: ',fwhm,lmax_smooth
     print, '*****************************'

     if keyword_set(mask) then map=map*mask
     ismoothing,map, map1,fwhm_arcmin=fwhm,lmax=lmax_smooth,/ring,/silent,regression=2
     map=map1
     minval = min(map)
     maxval = max(map)     
  endif

  if keyword_set(mask) then begin
     ipmask=where(mask lt 0.5)
     init_healpix
     map[ipmask] = !healpix.bad_value
  endif
  
  ;;=======================
  mct=2
  mcs=scale_sz*1.2
  cgcs=scale_sz*1.2
  cgct=2
  !P.CHARTHICK=2

  tcgcs=scale_sz*1.2
  tcgct=1.5

  mnep='+'
  melse='.'
  mbiposh='x'

  times="264B
  plus="53B
  dot="267B
  star="52B
  diamond="340B
  diamondFilled="250B
  nabla="104B
  grad="321B
  box="200B
  otimes="304B
  oplus="305B
  
  mwmap='!9' + String(star) + '!X'
  mlowl='!9' + String(dot) + '!X'
  mbiposh='!9' + String(times) + '!X'


  mdipole='!9' + String(otimes) + '!X'
  mnep='!9' + String(plus) + '!X'

  ;;=======================
  
  ;; mwmap=14
  ;; mdipole=15
  ;; mlowl=16
  ;; mnep=1
  ;; melse=2

if not keyword_set(minval) then minval=!NULL
if not keyword_set(maxval) then maxval=!NULL




  if not keyword_set(fclose) then begin
     planck_mollview_window,out,CTDIR,CTFILE,GR,W,PX,large=large,medium=medium,NColors=ncolors ,bottom=1
     ;;no graticule
     if glsize eq 0 then GR=0

     if keyword_set(file_ps) then begin
        ls_mollview, map, COLT=colt, CTDIR=CTDIR, CTFILE=CTFILE, MIN=minval, MAX=maxval, $
                     CHARSIZE=mcs,charthick=mct,GRATICULE=GR, GLSIZE=glsize,HXSIZE=W, FLIP=0, hist=hist,bwhite=bwhite,$
                     GRMIN=[-179, -89], GRMAX=[179,89], GRLS =grls, TITLE=title, PXSIZE=PX, CBLBL=ctitle, $
                     nobar=nobar,PS=file_ps,/keep_file_open  ,OUTLINE=out,_extra=extra ;execute='dipole_directions' ;;
     endif else begin
        ls_mollview, map, COLT=colt, CTDIR=CTDIR, CTFILE=CTFILE, MIN=minval, MAX=maxval, CHARSIZE=mcs,$
                     charthick=mct,HXSIZE=W, hist=hist,bwhite=bwhite,$
                     GRMIN=[-179, -89], GRMAX=[179,89], TITLE=title, PXSIZE=PX,nobar=nobar,PS=file_ps,CBLBL=ctitle, $ 
                     /keep_file_open  ,OUTLINE=out,_extra=extra ;execute='dipole_directions'
     endelse
  endif


if not keyword_set(minval) then minval=min(map)
if not keyword_set(maxval) then maxval=max(map)

nnval=5
if keyword_set(val) then nnval=n_elements(val)
if not keyword_set(val) then val=linspace(minval,maxval,nnval)
if keyword_set(ncolors) then discrete=1

if not keyword_set(fmtbar) then fmtbar='(g0.1)'
if nobar eq 1 then cgColorbar, NColors=ncolors,  Bottom=1, discrete=discrete,Ticklen=0.25,Range=fix([minval,maxval]), divisions=nnval-1,$ 
            position=[0.2,0.12,0.8,0.15],title=ctitle,TCHARSIZE=0.8,tlocation='bottom',charsize=0.8,TICKNAMES=strtrim(string(val,format=fmtbar),2)


;; cgdcbar, NColors=ncolors, labels=strtrim(string(val),2), spacing=0.75,Bottom=1,$ 
;;             position=[0.2,0.12,0.8,0.15],title=ctitle,TCHARSIZE=0.8,charsize=0.8


     tek_color

  if keyword_set(thph_out) then begin

     if not keyword_set(stxpx) then stxpx=[0,1,1]

     ;;arbitrary point to plot together
     theta_out = thph_out[0]
     phi_out = thph_out[1]
     
     ;;wmap9 mean value
     ang2vec,theta_out,phi_out,vec
     vec2moll,vec,u,v
     cgtext,u,v,'+',col=white,charsize=cgcs,charthick=cgct,align=0.5 
     if keyword_set(fthph_out) then xyouts,u*0.95,v*1.2,fthph_out,col=stxpx[0],charsize=tcgcs,charthick=stxpx[2]*1
  endif
  

  if keyword_set(show_dipdir) then begin


     mmsize=[40,40,40,40]
     ys=[0.14,0.14,0.14,0.14]
     ccsize=[2,2,2,2]
     msize=5
     csize=1.2


     dd = ['*','*','*','*']

;;C-R    green        #008000 =>  0;128;0
;;NILC   deepskyblue  #00BFFF =>  0;191;255
;;SEVEM  red          #FF0000 =>  255;0;0
;;SMICA  orange       #FFA500 =>  255;165;0
     tvlct, 0B,   128B, 0B,  11
     tvlct, 0B,   191B, 255B,12
     tvlct, 255B, 0B,   0B,  13
     tvlct, 255B, 165B, 0B,  14

     cc = [11,12,13,14]

     ;;wmap9 total mean theta and phi upto lmax=600
     ;theta2 = [2.0366871]
     ;phi2 = [3.9453986]

     ;;mean of the individual 16l angles
     theta2 = [1.98968]
     phi2 = [3.64774]

     ;;wmap9 mean value
     ang2vec,theta2[0],phi2[0],vec
     vec2moll,vec,u,v
     cgtext,u,v,mwmap,col=white,charsize=cgcs,charthick=cgct,align=0.5 
     xyouts,u*0.9,v*1.2,'wmapix',col=white,charsize=tcgcs,charthick=tcgct

     ;oplot,u,v,psym=symcat(wmap),col=white,charsize=cgcs,charthick=cgct  ;;,align=0.5
;;sss = SYMBOL(u, v, '*', /data,SYM_COLOR='Blue', SYM_SIZE=2,LABEL_FONT_STYLE='bi',LABEL_FONT_SIZE=2,LABEL_FONT_BACKGROUND=1, $
;;  SYM_THICK=3, LABEL_STRING='wmapix',/SYM_FILLED)

     ;xyouts,u,v-0.03,'x',col=0,charsize=scale_sz*csize,charthick=msize


;;dipole modulation
     dm_th = -15.0947870995432 
     dm_ph = 227.460937500000 
     theta_dm=dm_th             ;*!dtor
     phi_dm=dm_ph               ;*!dtor
     ang2vec,theta_dm,phi_dm,vec ,/astro
     vec2moll,vec,u,v
     cgtext,u,v,mlowl,col=white,charsize=cgcs,charthick=cgct,align=0.5
     xyouts,u*0.9,v*0.65,'lowell',col=white,charsize=tcgcs,charthick=tcgct

     ;xyouts,u-0.01,v-0.03,'o',col=cc[0],charsize=scale_sz*csize,charthick=msize

;;CMB dipole direction
;;http://www.sciencedirect.com/science/article/pii/S1387647306001990
     theta = 90d0-41.75         ;*!dtor
     phi = 263.85               ;*!dtor
     ang2vec,theta,phi,vec,/astro
     vec2moll,vec,u,v
     cgtext,u,v,mdipole,col=white,charsize=cgcs,charthick=cgct,align=0.5
     cgtext,u*0.8,v*1.15,'Dipole',col=white,charsize=tcgcs,charthick=tcgct
     

;;Eclipticl poles
     ang2vec,-29.81, 276.38,vecsep,/astro
     vec2moll,vecsep,u,v
     xyouts,u,v,mnep,col=0,align=0.5,charsize=cgcs,charthick=5
     xyouts,u-0.01,v+0.05,'SEP',col=0,align=0.5 ,charsize=tcgcs,charthick=tcgct
     
     ang2vec,29.81, 96.38,vecnep,/astro
     vec2moll,vecnep,u,v
     xyouts,u,v,mnep,col=0,align=0.5 ,charsize=cgcs,charthick=5
     xyouts,u-0.01,v+0.05,'NEP',col=0,align=0.5 ,charsize=tcgcs,charthick=tcgct

     ;;Planck, C-R, NILC, SEVEM and SMICA dipole directions
     if keyword_set(method_vec) then begin
        planck_th=[-0.51,1.41,3.58,3.01]
        planck_ph=[212.61, 209.48, 210.33, 211.13] 
        name=['C-R', 'NILC', 'SEVEM','SMICA']
        planck_ct,vcol
        xshift=[0.9,0.9,0.9,0.9]
        yshift=[0.95,1.3,0.95,0.95]

        for jj=0,n_elements(method_vec)-1 do begin
           ii = method_vec[jj]
           ang2vec,planck_th[ii],planck_ph[ii],vec ,/astro
           vec2moll,vec,u,v
           xyouts,u,v,'.',col=vcol[ii],align=0.5 ,charsize=scale_sz*5,charthick=5
           ;xyouts,u*xshift[ii],v*yshift[ii],name[ii],col=vcol[ii],align=0.5 ,charsize=scale_sz*0.3,charthick=1
        endfor
     endif

  endif 

  if keyword_set(txpx) then begin
     ;;wmap1 mean theta and phi upto lmax=600
     if not keyword_set(stxpx) then stxpx=[0,1,1,0,2]

     print, 'stxpx = ',stxpx

     ;;xsc = [0.96, 0.9, 0.9, 1.0, 0.96, 0.9, 0.99]
     ;;ysc=  [0.7,  1.4, 0.8, 0.89,  1.3,  0.7, 1.2]

     if not keyword_set(xsc) then xsc = [1, 3, 1, 1, 1,  1, 0]
     if not keyword_set(ysc) then  ysc=  [-1, -1, 0, 0, 0, -1, -1]
     if not keyword_set(tprad) then  tprad = [3,3,3,3,3,4,4]-0.5

    if not keyword_set(fftxpx) then  fftxpx=replicate('txpx1',n_elements(txpx[0,*]))

     for i=0,n_elements(txpx[0,*])-1 do begin
        theta2 = [txpx[0,i]*!dtor]
        phi2 = [txpx[1,i]*!dtor]
        
        ;; mean value
        ang2vec,theta2[0],phi2[0],vec
        vec2moll,vec,u,v

        if not stxpx[4] eq 0 then cgtext,u,v,'.',col=stxpx[3],charsize=scale_sz*stxpx[4],charthick=10,align=0.5
        if not keyword_set(fftxpx) then fftxpx[i]='txpx'+strn(100*i)

        ang2vec,theta2[0]+ysc(i)*tprad(i)*!dtor,phi2[0]+xsc(i)*tprad(i)*!dtor,vec
        vec2moll,vec,u2,v2
        
        cgtext,u2,v2,mbiposh,col=white,charsize=cgcs,charthick=cgct,align=0.5
        xyouts,u2*0.95,v2*1.15,fftxpx[i],col=stxpx[0],charsize=scale_sz*stxpx[1]*0.5,charthick=stxpx[2]*1
     endfor
  endif

  if keyword_set(lowell) then begin
     ;;wmap1 mean theta and phi upto lmax=600
     theta2 = [117d0*!dtor]
     phi2 = [225d0*!dtor]

     ;;wmap1 mean value
     ang2vec,theta2[0],phi2[0],vec
     vec2moll,vec,u,v
     ;;cgtext,u,v,mlowl,col=white,charsize=cgcs,charthick=cgct,align=0.5
     plot,u,v,psym=symcat(lowl),col=white,charsize=cgcs,charthick=cgct ;;,align=0.5
     xyouts,u*0.96,v*0.825,'lowell',col=white,charsize=tcgcs,charthick=tcgct



     ;;CMB dipole direction
     ;;http://www.sciencedirect.com/science/article/pii/S1387647306001990
     theta = 90d0-41.75         ;*!dtor
     phi = 263.85               ;*!dtor
     ang2vec,theta,phi,vec,/astro
     vec2moll,vec,u,v
     ;;cgtext,u,v,mdipole,col=white,charsize=cgcs,charthick=cgct,align=0.5
     plot,u,v,psym=symcat(dipole),col=white,charsize=cgcs,charthick=cgct ;;,align=0.5
     cgtext,u*0.95,v*0.9,'Dipole',col=white,charsize=tcgcs,charthick=tcgct
     

     ;;Eclipticl poles
     ang2vec,-29.81, 276.38,vecsep,/astro
     vec2moll,vecsep,u,v
     ;;xyouts,u,v,mnep,col=0,align=0.5 ,charsize=cgcs,charthick=5
     plot,u,v,psym=symcat(mnep),col=white,charsize=cgcs,charthick=cgct ;;,align=0.5
     xyouts,u-0.01,v+0.05,'SEP',col=0,align=0.5 ,charsize=tcgcs,charthick=tcgct
     
     ang2vec,29.81, 96.38,vecnep,/astro
     vec2moll,vecnep,u,v
     ;;xyouts,u,v,mnep,col=0,align=0.5 ,charsize=cgcs,charthick=5
     plot,u,v,psym=symcat(mnep),col=white,charsize=cgcs,charthick=cgct ;;,align=0.5
     xyouts,u-0.01,v+0.05,'NEP',col=0,align=0.5 ,charsize=tcgcs,charthick=tcgct
  endif

  if keyword_set(ecp) then begin
     ;;CMB dipole direction
     ;;http://www.sciencedirect.com/science/article/pii/S1387647306001990
     theta = 90d0-41.75         ;*!dtor
     phi = 263.85               ;*!dtor
     ang2vec,theta,phi,vec,/astro
     vec2moll,vec,u,v
     ;;cgtext,u,v,mdipole,col=white,charsize=cgcs,charthick=cgct
     plot,u,v,psym=symcat(dipole),col=white,charsize=cgcs,charthick=cgct ;;,align=0.5
     cgtext,u*0.95,v*0.9,'Dipole',col=white,charsize=tcgcs,charthick=tcgct
     

     ;;Eclipticl poles
     ang2vec,-29.81, 276.38,vecsep,/astro
     vec2moll,vecsep,u,v
     ;;xyouts,u,v,mnep,col=0,align=0.5 ,charsize=cgcs,charthick=5
     plot,u,v,psym=symcat(mnep),col=white,charsize=cgcs,charthick=cgct ;;,align=0.5
     xyouts,u-0.01,v+0.05,'SEP',col=0,align=0.5 ,charsize=tcgcs,charthick=tcgct
     
     ang2vec,29.81, 96.38,vecnep,/astro
     vec2moll,vecnep,u,v
     ;;xyouts,u,v,mnep,col=0,align=0.5 ,charsize=cgcs,charthick=5
     plot,u,v,psym=symcat(mnep),col=white,charsize=cgcs,charthick=cgct ;;,align=0.5
     xyouts,u-0.01,v+0.05,'NEP',col=0,align=0.5 ,charsize=tcgcs,charthick=tcgct
  endif


  if keyword_set(mtheta) and keyword_set(mphi) then begin
     nval = n_elements(mtheta)
     if not keyword_set(val) then begin
       val = indgen(nval)
        for ii=0,nval-1 do begin
           val[ii] = max([minval, midval+ii*100+minval])
        endfor
     endif

     ;xdir_str=strn(val[0:nval-1])
     ;ydir_str = xdir_str

   for ii=0,nval-1 do begin
      theta = mtheta[ii]
      phi = mphi[ii]
      ang2vec,theta,phi,vec
      vec2moll,vec,u,v
      ell_val = val[ii]
      ;xdir_str[ii]='ddir'+strn(ii)

      csize=1
      ecol=ellcol
      esize=ellsize
      if keyword_set(ddir) then begin
         if ii lt 6 then begin
            csize=1.4
            esize=1.0*ellsize
            ;ecol=4
         endif else begin
            csize=1.4
            esize=1.0*ellsize
         endelse
      endif

      cgtext,u-0.1,v,strn(ell_val),col=ecol,charsize=scale_sz*esize,charthick=csize
   endfor
  endif

  if keyword_set(tkey) then begin
     xyouts,1.5,0.88,tkey,col=2,align=0.4 ,charsize=cgcs,charthick=5
  endif

if keyword_set(file_ps) then begin

  if not keyword_set(fopen) then begin

     axisx = ['-120','-60','0','60','120','-45','45','Central Multipole','Dipole','lowell','wmapix','microk','txpx','SEP','NEP']
     axisy = ['240$^\circ$','300$^\circ$','0','60$^\circ$','120$^\circ$','-45$^\circ$','45$^\circ$','$\ell_{\rm central}$',$
              '{\bf \rm CMB dipole}','{\bf \rm low-}$\ell$','{\bf \rm {\it WMAP}-9}','$\mu$K','Input','{\bf \rm SEP}','{\bf \rm NEP}']

     ;;charsz = [replicate(1.8,8),2.2,2.2,2.0,2.2,2.2,2.2,2.2]
     charsz = replicate(2.5,n_elements(axisx))



     if keyword_set(ax) then    axisx = [axisx[*],ax[*]]
     if keyword_set(ay) then    axisy = [axisy[*],ay[*]]

     if keyword_set(ax) then begin
        exsz = replicate(3,n_elements(ax))
        charsz = [charsz,scale_sz*exsz]
     endif

     if keyword_set(ftxpx) then begin
        axisx=[axisx[*],fftxpx[*]]
        axisy=[axisy[*],ftxpx[*]]
        exsz = replicate(stxpx[1],n_elements(ftxpx))
        charsz = [charsz,scale_sz*exsz]
     endif

     if keyword_set(fthph_out) then begin
        axisx=[axisx[*],fthph_out]
        axisy=[axisy[*],fthph_out]
        exsz = replicate(stxpx[1,0],n_elements(fthph_out))
        charsz = [charsz,scale_sz*exsz]
     endif

     print, axisx
     print, axisy
     print, 'caharsize = ',charsz


     ps_end_planck,png=png,feps=file_ps,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy,charsize=charsz,cl=cl,/nops ;,mollview=mollview
  endif

endif

end
