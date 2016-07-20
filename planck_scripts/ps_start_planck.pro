pro ps_start_planck,file=file,medium=medium,large=large,mollview=mollview,ytdx=ytdx,xtdy=xtdy,xmar=xmar,ymar=ymar,xsize=xsize,ysize=ysize,font=font

  X1 = 88d/25.4d                ; 88 mm by 25.m mm/inch to get the size in inches
  AR = 1.5d                     ; an aspect ratio of 1.5
  Y1 = X1/AR
                                ;
  X2 = 120d/25.4d               ; 120 mm
  Y2 = X2/AR
                                ;
  X3 = 180d/25.4d               ; 
  Y3 = X3/AR
                                ;
  X4 = X1
  Y4 = Y1*2                     ; this will be a 2 x 1 plot.
                                ;
  X5 = X2
  Y5 = Y2*2                     ; this will be a 2 x 1 plot.
                                ;
  X6 = X3
  Y6 = Y3*2
                                ;
  colours                       ;	This routine is included in this file (see bottom) and sets up a few default colours.
  Dname = !D.NAME               ;	This is different for windows and linux (and likely MAC)

;;  XMAR_ = [12,2.25]           ; Set this to minimize white space on left and right sides of figure.                                ;
  XMAR_ = [8.15,2.25]           ; Set this to minimize white space on left and right sides of figure.
  XMAR_ = [6,2.25]
  YMAR_ = [3,0.5]               ; Set this to minimize white space on bottom and top of figure.

                                ;
  XMAR1 = XMAR_
  XMAR2 = XMAR_*120d/88d
  XMAR3 = XMAR_*180d/88d
                                ;
  YMAR1 = YMAR_
  YMAR2 = YMAR_*120d/88d
  YMAR3 = YMAR_*180d/88d
                                ;
  YTTL_DXs = [-50d, -50d*88d/120d, -50d*88d/180d]
  XTTL_DYs = [-100d, -100d*88d/120d, -100d*88d/180d]
                                ;
  FNTsz = 8                     ; Set this to the desired font size (pt)
                                ;
  if not keyword_set(font) then font=FNTsz
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   I have the sizes for 5 figures mapped out:
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  small, medium, and large single figures, and a small and large 1x2 figure pair.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  Do the first three as a for loop, and the last two as custom runs.
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                                ;
                                ;
  SZstrs = ['88','120','180']
  Xs = [X1,X2,X3]
  Ys = [Y1,Y2,Y3]
                                ;
  jj=0
  if keyword_set(medium) then jj=1
  if keyword_set(large) then jj=2
  Xsz = Xs[jj]
  Ysz = Ys[jj]

  CASE jj OF
     0: XMAR = XMAR1
     1: XMAR = XMAR2
     2: XMAR = XMAR3
  ENDCASE
                                ;
  CASE jj OF
     0: YMAR = YMAR1
     1: YMAR = YMAR2
     2: YMAR = YMAR3
  ENDCASE
                                ;
  SZstr = SZstrs[jj]

xmar=XMAR
ymar=YMAR
xtdy= XTTL_DYs[jj]
ytdx= YTTL_DXs[jj]

  set_plot, "ps"
  !P.font = 0
  if not keyword_set(mollview) then  begin
     xsize=Xsz
     ysize=Ysz
     device, FILENAME=file,xsize=xsize,ysize=ysize,/inches,/color, /ENCAPSULATED, /HELVETICA, FONT_size=FNTsz
  endif else begin
     hxsize = 15.
     xsize = Xsz
     ysize = Ysz
     papersize = 'a4'

     du_dv = 2.                 ; aspect ratio
     fudge = 1.02               ; spare some space around the Mollweide egg
     xc = 0.5*(xsize-1) & delta_x = (xsize-1 - xc)
     yc = 0.5*(ysize-1) & delta_y = (ysize-1 - yc) 
     ;; x and y range of egg
     umin = - du_dv * fudge & umax = du_dv * fudge
     vmin = - fudge         & vmax =         fudge
     ;; position of the egg in the final window
     w_xll = 0.0 & w_xur = 1.0 & w_dx = w_xur - w_xll
     w_yll = 0.1 & w_yur = 0.9 & w_dy = w_yur - w_yll
     w_dx_dy = w_dx / w_dy      ; 1./.8
     ;; color bar, position, dimension
                                ;cbar_dx = 1./3.
                                ;cbar_dy = 1./70.
     cbar_dx = 2./3.
     cbar_dy = 1./32.
     cbar_xll = (1. - cbar_dx)/2.
     cbar_xur = (1. + cbar_dx)/2.
     cbar_yur = w_yll - cbar_dy
     cbar_yll = cbar_yur - cbar_dy
     ;; polarisation color ring, position, dimension
     cring_dx = 1./10.
     cring_dy = 1./10.
     cring_xll = .025

     cring_yll = .025
     ;; location of pol vector scale
     vscal_x = 0.05
     vscal_y = 0.02
     ;; location of title and subtitle
     x_title = 0.5 & y_title = 0.95
     x_subtl = 0.5 & y_subtl = 0.905

     xsize=hxsize
     ysize=hxsize/du_dv*w_dx_dy

     yoffset = (papersize eq 'a4') ? 2 : 1
     DEVICE, FILENAME=file, XSIZE=hxsize, YSIZE=ysize, XOFFSET=4, YOFFSET=hxsize+yoffset, /ENCAPSULATED, /HELVETICA, FONT_size=font
  endelse
  ;; ;
  ;; !P.CHARTHICK = 1d
  ;; !P.CHARSIZE=1						;	Set the charactersize to not be scaled from that above.
  ;; !X.CHARSIZE=1						;	Set the X-label the same as the main figure text.
  ;; !Y.CHARSIZE=1						;	Set the Y-label the same as the main figure text.
  ;; !p.thick = 1.0d						;	Set the lines a bit thicker the nthe minimum of 1 pt
  ;; !x.thick = !P.thick 					;	Set x-axis lines the same as others within the plot
  ;; !y.thick = !P.thick 					;	The same for y-axis lines

end
