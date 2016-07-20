pro plot_gnomview,map,minval,maxval,large=large,medium=medium,fopen=fopen,show_dipdir=show_dipdir,mtheta=mtheta,mphi=mphi,$
                  file_ps=file_ps,ctitle=ctitle,title=title,cl=cl,ax=ax,ay=ay,fclose=fclose,help=help,$
                  midval=midval,tkey=tkey,colt=colt,hist=hist,ellcol=ellcol,ellsize=ellsize,method_vec=method_vec

  if keyword_set(help) then begin
     print, 'USAGE: '
     print, 'plot_mollview,map,minval,maxval,large=large,medium=medium,fopen=fopen,file_ps=file_ps,ctitle=ctitle,title=title,cl=cl,ax=ax,ay=ay,fclose=fclose'
  endif

  resolve_routine,'ls_mollview',/COMPILE_FULL_FILE,/either

;;plot dipole directions 
  if not keyword_set(ctitle) then ctitle='Central Multipole'
  if not keyword_set(title) then title=' '
  if not keyword_set(midval) then midval = 50
  if not keyword_set(colt) then colt = 41
  if not keyword_set(ellcol) then ellcol=0
  if not keyword_set(ellsize) then ellsize=0.4

  if midval lt 1 then midval=0

  SZ = 88d
  if keyword_set(medium) then SZ=120d
  if keyword_set(large) then SZ=180d

  ysize=SZ/10d/1.6d
  xsize=SZ/10d

  if not keyword_set(fclose) then begin
     planck_mollview_window,out,CTDIR,CTFILE,GR,W,PX,large=large,medium=medium

     triangle= create_struct('coord','G','ra',[0,80,0],'dec',[40,45,65])

     LS_gnomview, map, res=10,rot=[212,86],graticule=[45,30],COLT=colt, CTDIR=CTDIR, CTFILE=CTFILE, MIN=minval, MAX=maxval, CHARSIZE=0.8,charthick=3, hist=hist,$
                  GRMIN=[-179, -89], GRMAX=[179,89], GRLS =1, PXSIZE=PX, HXSIZE=W,TITLE=title, CBLBL=ctitle,PS=file_ps,/keep_file_open  ,OUTLINE=out ;execute='dipole_directions'
  endif


  if keyword_set(show_dipdir) then begin

     tek_color

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

     ;;wmap9 mean theta and phi upto lmax=600
     theta2 = [2.0366871]
     phi2 = [3.9453986]

     ;;wmap9 mean value
     ang2vec,theta2[0],phi2[0],vec
     vec2moll,vec,u,v
     cgtext,u,v,'.',col=1,charsize=2,charthick=2
     xyouts,u*0.95,v*1.2,'WMAP9',col=1,charsize=0.5,charthick=1

     ;xyouts,u,v-0.03,'x',col=0,charsize=csize,charthick=msize


;;dipole modulation
     dm_th = -15.0947870995432 
     dm_ph = 227.460937500000 
     theta_dm=dm_th             ;*!dtor
     phi_dm=dm_ph               ;*!dtor
     ang2vec,theta_dm,phi_dm,vec ,/astro
     vec2moll,vec,u,v
     cgtext,u,v,'.',col=1,charsize=2,charthick=2
     xyouts,u*0.96,v*0.86,'lowell',col=1,charsize=0.5,charthick=1

     ;xyouts,u-0.01,v-0.03,'o',col=cc[0],charsize=csize,charthick=msize

;;CMB dipole direction
;;http://www.sciencedirect.com/science/article/pii/S1387647306001990
     theta = 90d0-41.75         ;*!dtor
     phi = 263.85               ;*!dtor
     ang2vec,theta,phi,vec,/astro
     vec2moll,vec,u,v
     cgtext,u,v,'.',col=1,charsize=2,charthick=2
     cgtext,u*0.95,v*0.9,'dipole',col=1,charsize=0.5,charthick=1
     

;;Eclipticl poles
     ang2vec,-29.81, 276.38,vecsep,/astro
     vec2moll,vecsep,u,v
     xyouts,u,v,'.',col=0,align=0.5 ,charsize=2,charthick=5
     xyouts,u-0.01,v+0.05,'SEP',col=0,align=0.5 ,charsize=0.5,charthick=1
     
     ang2vec,29.81, 96.38,vecnep,/astro
     vec2moll,vecnep,u,v
     xyouts,u,v,'.',col=0,align=0.5 ,charsize=2,charthick=5
     xyouts,u-0.01,v+0.05,'NEP',col=0,align=0.5 ,charsize=0.5,charthick=1

     ;;Planck, C-R, NILC, SEVEM and SMICA dipole directions
     if keyword_set(method_vec) then begin
        planck_th=[-9.,3.,4.,4.]
        planck_ph=[201.,210.,211.,212.]
        name=['C-R', 'NILC', 'SEVEM','SMICA']
        planck_ct,vcol
        xshift=[0.9,0.9,0.9,0.9]
        yshift=[0.95,1.3,0.95,0.95]

        for jj=0,n_elements(method_vec)-1 do begin
           ii = method_vec[jj]
           ang2vec,planck_th[ii],planck_ph[ii],vec ,/astro
           vec2moll,vec,u,v
           xyouts,u,v,'.',col=vcol[ii],align=0.5 ,charsize=5,charthick=5
           ;xyouts,u*xshift[ii],v*yshift[ii],name[ii],col=vcol[ii],align=0.5 ,charsize=0.3,charthick=1
        endfor
     endif

  endif


  if keyword_set(mtheta) and keyword_set(mphi) then begin
     nval = n_elements(mtheta)
   for ii=0,nval-1 do begin
      theta = mtheta[ii]
      phi = mphi[ii]
      ang2vec,theta,phi,vec
      vec2moll,vec,u,v
      ell_val = max([minval, midval+ii*100])
      cgtext,u,v,strn(ell_val),col=ellcol,charsize=ellsize,charthick=1   
   endfor
  endif

  if keyword_set(tkey) then begin
     xyouts,1.5,0.88,tkey,col=2,align=0.4 ,charsize=2,charthick=5
  endif

  if not keyword_set(fopen) then begin
     exsz = [3]
     axisx = ['-120','-60','0','60','120','-45','45','Central Multipole','lowell']
     axisy = ['240$^\circ$','300$^\circ$','0$^\circ$','60$^\circ$','120$^\circ$','-45$^\circ$','45$^\circ$','$\ell_{\rm central}$','low-$\ell$']
     csdefault = replicate(1.5,n_elements(axisx))

     if keyword_set(ax) then    axisx = ['-120','-60','0','60','120','-45','45','lowell',ax[*]]
     if keyword_set(ay) then    axisy = ['240$^\circ$','300$^\circ$','0$^\circ$','60$^\circ$','120$^\circ$','-45$^\circ$','45$^\circ$','low-$\ell$',ay[*]]
     if keyword_set(ax) then exsz = replicate(3,n_elements(ax))

     print, axisx
     print, axisy


     charsz = [csdefault,exsz]
     ps_end_planck,/png,feps=file_ps,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy,charsize=charsz,cl=cl
  endif



end
