pro planck_plot,xin,yin,fname=fname,cl=cl,ax=ax,ay=ay,xra=xra,yra=yra,nodata=nodata,psend=psend,xtitle=xt,ytitle=yt,xs=xs,ys=ys,pmulti=pmulti,$
                xgrid=xgrid,ygrid=ygrid,ny=ny,nmy=nmy,nx=nx,nmx=nmx,color=color,noeras=noeras,yover=yover,dops=dops,line=line, legtxt=legtxt,$
                thick=thick,_extra=extra,sxyra=sxyra,angle=angle,format=format,medium=medium,large=large,DECMODX=DECMODX, DECMODY=DECMODY,hline=hline

resolve_routine,'hfi_plot',/COMPILE_FULL_FILE,/either


;; !P.CHARTHICK = 1d
;; !P.CHARSIZE=1                   ;	Set the charactersize to not be scaled from that above.
;; !X.CHARSIZE=1                   ;	Set the X-label the same as the main figure text.
;; !Y.CHARSIZE=1                   ;	Set the Y-label the same as the main figure text.
;; !p.thick = 1.0d                 ;	Set the lines a bit thicker the nthe minimum of 1 pt
;; !x.thick = !P.thick             ;	Set x-axis lines the same as others within the plot
;; !y.thick = !P.thick             ;	The same for y-axis lines

if n_elements(angle) eq 0 then angle=90d
;;print, 'yaxis_label_angle=',angle
;;if not keyword_set(format) then format='(f3.1)'

if keyword_set(fname) then  ps_start_planck,file=fname,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xs,ysize=ys,medium=medium,large=large ;, /large ;,/medium ;large

planck_ct,vcol


xmat=xin
ymat=yin
x=xin
y=yin

sizexmat=size(x)
sizeymat=size(y)
ncase=n_elements(ymat[0,*]) 

;;make xmat and ymat to have same rank
if (sizexmat[0] ne sizeymat[0]) and sizexmat[0] ne 1 then on_error,2
if (sizexmat[0] eq 1) and (sizeymat[0] ne 1) then begin
   xmat=ymat
   for ii=0,ncase-1 do begin
      xmat[*,ii]=x
   endfor
endif


x = xmat[*,0]
y=ymat[*,0]



if keyword_set(ax) and not keyword_set(xt) then    xt = ax[0]
if keyword_set(ax) and not keyword_set(yt) then    yt = ax[1]

if not keyword_set(xt) then   xt = 'y'
if not keyword_set(yt) then   yt = 'x'


if keyword_set(cl) then begin
   xt = 'xttl'
   yt = 'yttl'
endif

 if not keyword_set(yra) then yra = [(min(y)),(max(y))]
;    if not keyword_set(ny) then ny=4
;    if not keyword_set(nmy) then nmy = 4    

if not keyword_set(xra) then  xra = [(min(x)),(max(x))]
;if not keyword_set(nx) then nx=4
;if not keyword_set(nmx) then nmx = 4



yformat=''
if max(y) gt 1e3 or max(y) lt 1e-3 then yformat='(f5.3)'
if max(x) gt 1e3 or max(x) lt 1e-3 then xformat='exponent'
if keyword_set(ylog) then yformat='exponent'

;;multiple plots supported
nplots=1
if keyword_set(pmulti) then begin
   !P.multi = pmulti
   nplots = pmulti[1]*pmulti[2]

   pstr = create_struct('p',!P,'x',!X,'y',!Y)
   pmulti=replicate(pstr, nplots)

endif

for ijk=0,nplots-1 do begin ;;multiplot loop refer http://www.idlcoyote.com/tips/oplot_pmulti.html

   xytpm=''

   if keyword_set(pmulti) then xytpm = strn(ijk+1)
   if keyword_set(sxyra) then begin
      yra=sxyra[ijk].yra & xra=sxyra[ijk].xra & format=sxyra[ijk].format & ny=sxyra[ijk].ny
   endif
   ;;print, 'format = ',format

   ;;x corrdinate
   if keyword_set(xlog) then begin
      xtickv = LogLevels(xra)
      xtickname=string(xtickv)
      nx = N_Elements(xtickv)-1
   endif 

;;y log case
   if keyword_set(ylog) then begin
      ytickv = LogLevels(yra)
      ytickname=string(ticksy)
      ny = N_Elements(ticksy)-1
   endif 



colorvec=[0,2,3,4,5,6,7]
if keyword_set(color) then begin
   if n_elements(color) eq 1 then colorvec[0]=color
   if n_elements(color) gt 1 then colorvec=color
endif

linevec=[0,0,0,0,0,0,0]
if keyword_set(line) then begin
   if n_elements(line) eq 1 then linevec[0]=line
   if n_elements(line) gt 1 then linevec=line
endif

   pstyle_plot,x,y,color=0,hline=hline, $
            /xs,xr=xra, yr=yra, /ys, BACKGROUND=255,thick=2,$
            YTITLE=yt+xytpm, XTITLE=xt+xytpm,$ , XMARGIN=xmar, YMARGIN=ymar ,$
            line=0,nodata=nodata, noeras=neras, _extra=extra,/YNOZERO,$
            xticks=nx,XMINOR=nmx,$ 
            yminor=nmy,yticks=ny 


   ;; yra = smart_yra(yramin, yramax,ytickv=yaxval,ytickname=yaxstr,format=format,/chatty)
   ;; axlabel,yaxval,orientation=90,/yaxis,format=format,ytype=yaxstr[0]
   
   if not keyword_set(nodata) then begin
      for ii=0,ncase-1 do begin
         oplot,xmat[*,ii],ymat[*,ii],color=colorvec[ii],line=linevec[ii]
      endfor
   endif
;,ytickv=ytickv,ytickname=ytickname,YTicklen=ygrid, YGridStyle=1,$
;            xtickv=xtickv,xtickname=xtickname 


;;            DECMODX=DECMODX, DECMODY=DECMODY,YTTL_DX = ytdx, XTTL_DY=xtdy,Y_DY=0,X_DY=0,angle=angle

;,xtickv=xtickv,xtickname=xtickname

;;, $
;;xtickv=indgen(14),xtickname=['','1','2','3','4','5','6','7','8','9','10','11','12',''],XTickFormat='(A1)'

   if keyword_set(yover) then oplot, x,yover, color=colorvec[1],line=linevec[1]


   if keyword_set(pmulti) then begin
      pmulti[ijk].p = !P &  pmulti[ijk].x=!X &  pmulti[ijk].y=!Y
   endif


endfor

if keyword_set(legtxt) then  legend,legtxt,/top,textcolor=colorvec, right=right,box=0,thick=4 ,pspacing=1,delimiter=''

if keyword_set(fname) then begin
   if keyword_set(psend) then begin
      ps_end_planck, /png,feps=fname,/latex,xsize=xs,ysize=ys ,axisx=ax,axisy=ay,cl=cl,dops=dops
   endif
endif


end
