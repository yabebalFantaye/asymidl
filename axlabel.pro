;-------------------------------------------------------------
; $Id$
;+
; NAME:
;        AXLABEL
;
; PURPOSE:
;        Put previously calculated axis labels onto the screen
;        at proper position. This routine was designed to work 
;        together with LOGLEVELS to produce fancy log plots.
;        It involves several coordinate transformations in order
;        to be device independent and take into account the 
;        character size. The user can specify a label format
;        and use 'external' formatting functions similar to
;        the [XYZ]TICKFORMAT keyword of PLOT.
;
; CATEGORY:
;        Plotting
;
; CALLING SEQUENCE:
;        AXLABEL,Value [,/XAxis] [,keywords]
;
; INPUTS:
;        VALUE -> A vector with the values to be labelled on the 
;             axis.
;
; KEYWORD PARAMETERS:
;        /XAxis -> If set, the labels are placed on the X achis
;             rather than on the Y axis
;
;        /YAxis -> Place the labels on the Y axis (this is the
;        default,
;             and this keyword is there for purely aesthetic reasons)
;
;        CHARSIZE -> The character size of the label
;
;        FORMAT -> An IDL format string (used as argument to the
;              STRING function) or the name of a function that returns
;              formatted labels. This function must accept three
;              arguments, the third of which is the current value
;              (see the online help to [XYZ]TICKFORMAT for more
;              details).
;              AXLABEL always passes 0 to the first two arguments.
;
;        _EXTRA  keywords are passed on to XYOUTS (e.g. COLOR or
;              ORIENTATION). Note that the ALIGN keyword value is 
;              determined automatically.
;
; OUTPUTS:
;        Axis labels without fuss.
;
; SUBROUTINES:
;        None.
;
; REQUIREMENTS:
;        A DATA coordinate system must be established by a previous
;        PLOT command.
;
; NOTES:
;        AXLABEL currently operates only on the left and bottom axes.
;
; EXAMPLE:
;          xrange = [0.3,3.0]   ; define axis range
;          yrange = [0.3,3.0]
;          plot,[1],xr=xrange,yr=yrange, $   ; do the plot
;          title='Logarithmic X axis, Logarithmic Y axis',$
;          xtickf='(a1)',ytickf='(a1)', /ylog,/xlog
;          ; important: turn the tick labeling off with
;          ?tickformat='(A1)'
;          xlblv = loglevels(xrange)   ; get nice label values (0.5,
;          1., 2.)
;          ylblv = loglevels(yrange)
;          axlabel,xlblv, /xaxis       ; plot the labels
;          axlabel,ylblv, /yaxis
;
; MODIFICATION HISTORY:
;        mgs, 10 Sep 1999: VERSION 1.00
;        mgs, 23 Sep 1999: - bug fix for log-log plots
;
;-
; Copyright (C) 1999, Martin Schultz, Max-Planck-Institut
; f. Meteorologie
; This software is provided as is without any warranty
; whatsoever. It may be freely used, copied or distributed
; for non-commercial purposes. This copyright notice must be
; kept with any copy of this software. If this software shall
; be used commercially or sold as part of a larger package,
; please contact the author.
; Bugs and comments should be directed to martin.schultz@dkrz.de
; with subject "IDL routine axlabel"
;-------------------------------------------------------------


pro axlabel,value,Charsize=Charsize,XAxis=XAxis,YAxis=YAxis, orientation=orientation, charscale=charscale, $
            Format=Format,_EXTRA=e
  
  rot=0
  if keyword_set(orientation) then rot=orientation
  
                                ; Error catching
  if (N_Elements(VALUE) eq 0) then begin
     message,'Must supply at least one label value to AXLABEL!'
  endif
  
                                ; Set default for CHARSIZE and FORMAT
  if (not keyword_set(CHARSIZE)) then $
     CHARSIZE = 1.
  if (not keyword_set(CHARSCALE)) then $
     CHARSCALE = 1.0
  if (not keyword_set(FORMAT)) then $
     FORMAT = '(g0.1)'


                                ; Format VALUE according to format string. If this string
                                ; does not begin with '(', it is assumed that the user has passed
                                ; a formatting function as for [XYZ]TICKFORMAT
                                ; However, only the third (NUMBER)
                                ; argument of this function is used
  
  if (STRPOS(FORMAT,'(') ne 0) then begin
     ValS = STRARR(N_Elements(VALUE))
     for j=0,N_Elements(VALUE)-1 do $
        ValS[j] = CALL_FUNCTION(FORMAT,0,0,VALUE[j])
  endif else $                  ; apply format string directly
     ValS = STRING(VALUE,format=FORMAT)
  
  ValS = STRTRIM(ValS,2)

  print, 'ValS nchar: ',strlen(ValS)
  
  if (keyword_set(XAxis)) then begin
     
                                ; Get y position for label
                                ; Subtract one character size
     PY = !Y.Window[0] 
     PYOFF = CONVERT_COORD(1,!D.Y_CH_SIZE*CHARSIZE,/DEVICE,/TO_NORMAL)
     PY = PY - 1.05*PYOFF[1]
     ;;print,'X:PY:',py
     PY = REPLICATE(PY,N_Elements(VALUE))
     
                                ; Convert data values to normalized x coordinates
     Y0 = !Y.CRANGE[0]
     if (!Y.TYPE eq 1) then $
        Y0 = 10.^Y0
     PX = CONVERT_COORD(VALUE,REPLICATE(Y0,N_Elements(VALUE)), $
                        /DATA,/TO_NORMAL)
     PX = PX[0,*]
     ;;print,'X:PX=',px
     
  endif else begin              ; Y axis label (default)
     
     numticks=N_Elements(VALUE)-1

     ;; Get x position for label
     PX = !X.Window[0] - 0.01
     PX = REPLICATE(PX,N_Elements(VALUE))


     ;;-----------------
     ;; Get y position for label

     ;; ;;NCHARW=CHARSIZE*strlen(ValS)
     ;; zzz = CONVERT_COORD(0,!D.Y_CH_SIZE,/DEVICE,/TO_NORMAL)
     ;; NCHARW = zzz[1]*CHARSIZE*strlen(ValS)


     ;; print, 'ntick,ValS',numticks, ': ', ValS
     ;; print, 'strlen(ValS): ',strlen(ValS)
     ;; print, 'zzz[1],CHARSIZE, 2.6*NCHARW[numticks]',zzz[1], CHARSIZE, 2.6*NCHARW[numticks]
     ;; print, 'NCHARW',NCHARW
     ;; print, '!Y.Window: ',!Y.Window

     ;; ymin = !Y.Window[0] - 0.2*NCHARW[0]
     ;; ymax = !Y.Window[1] - 2.8*NCHARW[numticks] ;


     ;; zzz = CONVERT_COORD(0,value,/DATA,/TO_NORMAL)
     ;; PYY = reform(zzz[1,*])

     ;; if rot gt 0 then begin
     ;;    ydiff=(ymax-ymin)/(numticks+1)        
     ;;    PY = PYY        
     ;;    PY[0]=ymin
     ;;    PY[numticks]=ymax
     ;;    for i=1,numticks-1,1 do begin
     ;;       coef = 0.5
     ;;       if i eq 3 then coef = 0.6
     ;;       if i eq 2 then coef = 0.45
     ;;       if i eq 1 then coef = 0.22
     ;;       PY[i]=PY[0] + i*ydiff + coef*NCHARW[i] 
     ;;    endfor

     ;; endif else PY=PYY

     ;; ; Get x position for label
     ;; PX = !X.OWindow[0] - 0.010
     ;; PX = REPLICATE(PX,N_Elements(VALUE))

     ;;-----------------

     
     ;; ; Convert data values to normalized coordinates and
     ;; ; subtract half the character size
     PYOFF = CONVERT_COORD(0,!D.Y_CH_SIZE*CHARSIZE,/DEVICE,/TO_NORMAL)
     PY =CONVERT_COORD(REPLICATE(!X.CRANGE[0],N_Elements(VALUE)),VALUE,  $
                        /DATA,/TO_NORMAL)
     PY = PY[1,*]-0.5*PYOFF[1]




  endelse
  
  
  XYOUTS,PX,PY,ValS,/NORMAL,orientation=rot,charsize=CHARSIZE,_EXTRA=e
  
  return
end
 
