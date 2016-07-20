pro chain_plot_dipole_angle,mm,fname,nsim=nsim

version='ddx9'

figdir = '../output/'+version+'/figures/'
dirout = '../output/'+version+'/data/'
datadir = '../output/data/'

spawn,'mkdir -p '+figdir


method = ['sevem','smica','ruler','nilc']
mask_kind = ['_sevem_','_smica_','_ruler_','_nilc_','_common_']

dcases = [2,3,0,1,1]
mcases = [4,4,4,4,2]

legtxt = ['C-R','NILC','SEVEM','SMICA','SMICA**']


if not keyword_set(nsim) then nsim = 500l

if size(mm,/type) eq 7 then begin
   ;;first row is for the data and the rest are for the nsim simulations
   mm = dblarr(nsim+1,n_elements(dcases))
   readcol,fdata, rul_c, nilc_c,sev_c, smi_c, smi_s ; rul_c,nil_c,smi_s
   mm[*,0]=rul_c
   mm[*,1]=nilc_c
   mm[*,2]=sev_c
   mm[*,3]=smi_c
   mm[*,4]=smi_s
endif

f=dblarr(nsim)
x = indgen(nsim)+1

tek_color

    

ps_start_planck,file=fname,xmar=xmar,ymar=ymar,ytdx=ytdx,xtdy=xtdy,xsize=xsize,ysize=ysize

;;C-R    green        #008000 =>  0;128;0
;;NILC   deepskyblue  #00BFFF =>  0;191;255
;;SEVEM  red          #FF0000 =>  255;0;0
;;SMICA  orange       #FFA500 =>  255;165;0
tvlct, 0B,   128B, 0B,  11
tvlct, 0B,   191B, 255B,12
tvlct, 255B, 0B,   0B,  13
tvlct, 255B, 165B, 0B,  14

vcol = [11,12,13,14]
vline = 0*indgen(n_elements(dcases))    

!P.CHARTHICK = 1d
!P.CHARSIZE=1                   ;	Set the charactersize to not be scaled from that above.
!X.CHARSIZE=1                   ;	Set the X-label the same as the main figure text.
!Y.CHARSIZE=1                   ;	Set the Y-label the same as the main figure text.
!p.thick = 1.0d                 ;	Set the lines a bit thicker the nthe minimum of 1 pt
!x.thick = !P.thick             ;	Set x-axis lines the same as others within the plot
!y.thick = !P.thick             ;	The same for y-axis lines
!P.font=2

ytit = textoidl('Dipoles dispersion angle (deg)')
xtit = textoidl('FFp6 simulation number')


axisx = ['xlbl','yttl']
axisy = ['Simulation Number',"Dipole's dispersion angle (deg)"]

hfi_plot,x,mm(1:nsim,4),color=0, charsize=0.8,thick=2,$
      xtitle='xlbl',/xs, xr=[1,nsim],yr=[30,100], /ys, BACKGROUND=255, $
      YTITLE='yttl', XMARGIN=xmar+0.1, YMARGIN=ymar,yticks=7 ;, YTTL_DX = ytdx, XTTL_DY=xtdy,yticks=6

oplot,x,mm(1:nsim,4),color=0
f(*)=mm(0,4)
oplot,f ,col=0,thick=2

for iii=0,n_elements(dcases)-2 do begin            
   oplot,x,mm(1:nsim,iii),color=vcol[iii]
   f(*)=mm(0,iii)
   oplot,f ,col=vcol[iii],thick=2
endfor

vcol = [11,12,13,14,0]
legend,legtxt,color=vcol,textcolor=vline,/bottom,box=0,line=vline,pspacing=1,thick=4


ps_end_planck, /png,feps=fname,/latex,xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy


end
