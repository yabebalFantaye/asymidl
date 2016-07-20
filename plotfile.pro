pro PLOTFILE, file, a, _EXTRA=extra, $
NODATA=nodata, XLOG=pXLOG,YLOG=pYLOG, $
XCOL=xcol,YCOL=ycol, KIND=kind, FREQUENCY=frequency,$
HISTO=histo,NBINS=nbins, GAUSSFIT=gauss, SHIFTVAL=shiftval,HISDATA=HISDATA, $
LEGEND=LEGEND, AXISCOLOR=AXISCOLOR, $
EXECUTE = pexecute, Y2AXIS=y2axis, X2AXIS=x2axis, ERRPLOT=pERRPLOT, $
NDATA=ndata,VALNDATA=valndata, MESHCOL=meshcol,MESHKIND=meshkind,$
MESHPOINTS=meshpoints,$
XRANGE=XRANGE,YRANGE=YRANGE,TITLE=TITLE,XTITLE=XTITLE,YTITLE=YTITLE
;+
; NAME:
;	PLOTFILE
; PURPOSE:
;	to plot a simple N column file with IDL
; CATEGORY:
;
; CALLING SEQUENCE:
;	plotfile,'myfile' [,idl_var]
; INPUTS:
;	myfile  name of the file with data (between quotes) or
;               name of an IDL variable (without quotes) or
;		idl "handle" or tree of handles. In the case that
;		input handle is a list, overplots all the  data.
; OPTIONAL INPUT PARAMETERS:
;	idl_var an idl named variable where to store the data
; KEYWORD PARAMETERS:
;		Other keywords than those described here will be passed
;		to the plot (or surface or contour when MESH plot is set)
;		routines. See the idl doc for a complete description of
;		possible keywords.
;
;		XCOL = column for x (default=1)
;		YCOL = column for y (default=2)
;		KIND = kind of plot:
;			1 normal (default)
;			2 overplots the previous plot
;			3 lineal logaritmic plot
;			4 log lin plot
;			5 log log plot
;		TITLE='top_title' title to be written at the top
;		XTITLE='x_title'  title of x axis
;		YTITLE='y_title'           y
;		XRANGE=[xmin,xmax] range of the x variable for the plot
;		YRANGE=[ymin,ymax]              y
;               HISTO when set to 1 perform the histogram with the fits
;			column of the input data
;               NBINS = number of bins for the histogram
;		HISDATA = fltarr(2,nppints) where to optionally store
;			the histogram data (only works is HISTO is set)
;               GAUSSFIT when set to 1 perform a gauusian fit of the data
;		NODATA = when set to 1 do not display any graphical output
;			(useful for storing the file data in IDL)
;		AXISCOLOR = color index for the axis (not working if
;			NDATA is set)
;		FREQUENCY   1 (default) read all the points
;			    2 read alternatively half of the points
;			    n read 1/n oner the whole set of points
;		XLOG,YLOG   as for other IDL plot procedures
;		SHIFTVAL	    a 2-dim vector containing the values of 
;				the shiftval of the graphic (H anv V 
;				respectively)
;		LEGEND	    String with legend to be instert. Uses the LEGEND
;				procedure of F K Knight (knight@ll.mit.edu)
;			    The LEGEND string must contain the arguments 
;				and keywords to be passed to LEGEND.
;			    For example:
;			    plotfile,'file',legend="['1','2'],psym=[4,5]"
;		X2AXIS = [min,max] a 2-elements vector with the limits
;				of an optional second X-axis (top axis).
;		Y2AXIS = [min,max] a 2-elements vector with the limits
;				of an optional second Y-axis (right axis).
;		ERRPLOT=[col1,col2] If set, overplots error bars with 
;			bottom(top) values in the column number col1 (col2).
;
;		The following keywords refer to 3D plot:
;		To set the MESH (3D) view, set one of the following:
;               NDATA = Sets 3D multi-set plot. NDATA is the number of data 
;			sets in the input file. When this keyword is set, 
;			additional keywords AZ, ZAXIS,ZRANGE and ZLOG 
;			(same meaning as in SURFACE) are allowed.
;               MESHCOL = Number of the column that holds the label values 
;			 for the 3D plot. If MESHCOL is negative, these
;			 values are output in the plot, otherwise they
;			 are labelled as 0,1,2...
;               MESHPOINTS = Number of the points for each individual graph.
;		The following keywords are valid if one of the previous
;		keywords (for setting 3D view) is set:
;               VALNDATA =  VALNDATA is the array with the third axis.
;		MESHKIND =	0 (default) lines
;				1 surface
;				2 contour
;
;		Obsolete keywords that are still valid but inheritated and
;		not explicitely set:
;
;		PSYM = symbol of the plot (0=cont. line, default):
;                       1 Plus sign
;                       2 Asterisk
;                       3 Period
;                       4 Diamond
;                       5 Triangle
;                       6 Square
;                       7 X
;                       8 (See IDL users guide, D-15)
;                       9 Undefined
;                      10 Histogram mode
;		THICK	    as for other IDL plot procedures
;               LINESTYLE = 0 Solid (default)
;                           1 Dotted
;                           2 Dashed
;                           3 Dash Dot
;                           4 Dash Dot Dot Dot
;                           5 Long Dashes
;		COLOR = color index for the plot 
;		TICKLEN	    as for other IDL plot procedures
;		CHARSIZE    as for other IDL plot procedures
;
; OUTPUTS:
;	a plot
; OPTIONAL OUTPUT PARAMETERS:
;	the data of the file  in the variable selected (idl_var)
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Straightforward.
;
; EXAMPLE OF USE OF 3D MULTISET PLOT (keywords ndata, meshcol and valndata)
;	
;	; imagine we have two curves with points (1,1), (2,6) and (3,3) 
;	; for the first one and (1,3), (2,7) and (3,0) for the second one. 
;	; The first curve is labelled with the value 0.5 and the second is 
;	;  labelled with 1.0.
;	; 1) creation of idl array holding the data:
;	tmp = fltarr(3,6)
;	tmp(0,*) = [1,2,3,1,2,3]		; col 1
;	tmp(1,*) = [0.5,0.5,0.5,1.0,1.0,1.0]	; col 2
;	tmp(2,*) = [1,6,3,3,7,0]		; col 3
;	; the following command makes the 3D plot
;	plotfile,tmp,ycol=3,meshcol=2
;	; same plot is obtained:
;	plotfile,tmp,ycol=3,ndata=2
;	; same plot is obtained:
;	plotfile,tmp,ycol=3,meshpoints=3
;	;
;	; if you want to, in addition, label the plots with the 
;	; right values use either:
;	plotfile,tmp,ycol=3,ndata=2,valndata=tmp(1,*)
;	; or
;	plotfile,tmp,ycol=3,meshcol=-2
;	; We see, in conclusion that the same results can be obtained 
;	; using either the NDATA, MESHCOL or MESHNPOINTS keywords, but the 
;	; use of MESHCOL or MESHNPOINTS is easier.
;	
; MODIFICATION HISTORY:
;	by M. Sanchez del Rio. ESRF. Grenoble, Oct 1991
;	21-Nov-1991 include PSYM keyword
;	16-Dec-1991 include XYRANGE keyword
;	07-Sep-1992 include HISTO option, few other modifications.
;	18-Nov-1992 MSR includes NDATA option
;	16-May-1993 MSR introduces COLOR option.
;	11-Jun-1993 MSR correct a bug with color in the 'ps' device option
;	29-Jun-1993 MSR removes NCOL keyword by using the functions
;			RASCII and STRPARSE of Mati Meron
;	30-Jul-1993 MSR adds linestyle keyword
;	17-Nov-1993 MSR adds frequency keyword
;	15-Dec-1993 MSR adds ticklen charsize and thick keywords
;	17-Dec-1993 MSR modify [xy]range keyword to act though ![xy].range
;	19-Dec-1993 MSR adds the shiftval keyword
;	27-May-1994 MSR adds the HISDATA keyword
;	10-Oct-1994 MSR allows the use of "handles"
;	01-Aug-1995 MSR adds /xlog and /ylog keywords (idl4.0)
;	08-Aug-1995 MSR adds legend and axiscolor keywords.
;	17-Aug-1995 MSR sort colors, cosmetics.
;	01-Sep-1995 MSR adds [x,y]2axis keywords.
;	13-Sep-1995 MSR adds ERRPLOT (error bars)
;	23-May-1996 MSR adds AZ,ZAXIS,ZRANGE and ZLOG keywords.
;	26-Aug-1997 MSR adds MESHCOL keyword and update doc.
;	28-Aug-1997 MSR Removes the [xyz]title, psym, thick, linestyle,
;		az,zaxia,zrange and zlog EXPLICIT keywords. They are
;		always available because we use since now the keyword
;		inheritance mechanism. Error managing using catch.
;	29-Aug-1997 MSR in MESH option change sur by transpose(sur) for
;		a correct position of titles.
;	01-Set-1997 MSR xyrange keyword not longer support. Inheritate
;		kwywords [xy]range, ticklen, charsize and color
;	03-Set-1997 MSR fix some problems wuth errplot, redefine explicit 
;		kws: [xy]range, [ xy]title
;	17-Feb-1998 MSR adds !err_string when an error is caught.
;	10-Mar-1998 MSR fix bug when using simultaneously /NODATA and /HISTO
;	25-Aug-1998 MSR adds "shade" option in meshkind.
;	17-Set-1998 MSR cosmetics with shiftval kw.
;	03-Mar-1999 MSR adds printer device together ps.
;
;-	
catch, error_status
if error_status ne 0 then begin
   message,/info,'error caught: '+!err_string
   if sdep(/w) then itmp = Dialog_Message(/Error,$
	'PLOTFILE: error caught: '+!err_string)
   catch, /cancel
   on_error,2
   goto,out
endif
npar = n_params()
if npar EQ 0 then begin
  message,/info,'Usage: plotfile,file'
  return
endif

;
; ========================== set initial values ===============
;
;
; color
;
if keyword_set(title) then begin
 title_old = !p.title
 !p.title=title
endif
if keyword_set(xtitle) then begin
 xtitle_old = !x.title
 !x.title=xtitle
endif
if keyword_set(ytitle) then begin
 ytitle_old = !y.title
 !y.title=ytitle
endif
if keyword_set(xrange) then begin
 xrange_old = !x.range
 !x.range=xrange
endif
if keyword_set(yrange) then begin
 yrange_old = !y.range
 !y.range=yrange
endif
if !d.name EQ 'PS' OR !d.name EQ 'PRINTER'  then begin
  if keyword_extra(extra,'color') then  begin
    if keyword_extra(extra,'color',/return_value) EQ 255 then $
      extra.color = 0
  endif
  if keyword_set(axiscolor) then  begin
    if axiscolor EQ 255 then axiscolor=0
  endif
endif

;
; log axis def
;
if not(keyword_set(kind)) then kind = 1
case kind of
  1: begin
	xlog = 0 & ylog = 0
     end
  2: begin 
        axiscolor=0
	xlog = 0 & ylog = 0
     end
  3: begin
	xlog = 0 & ylog = 1
     end
  4: begin
	xlog = 1 & ylog = 0
     end
  5: begin
	xlog = 1 & ylog = 1
     end
  else: message,'kind set to a non valid value'
endcase
if keyword_set(pxlog) then xlog=1
if keyword_set(pylog) then ylog=1

if (keyword_set(NDATA) OR keyword_set(MESHCOL) OR keyword_set(MESHPOINTS)) $
  then imesh = 1 else imesh = 0
if not(keyword_set(xcol)) then xcol = 1
if not(keyword_set(ycol)) then ycol = 2
;if not(keyword_set(shiftval)) then shiftval = [0.0, 0.0]


input_type = (size(file)) ((size(file)) (0)+1)
;
; ========================== end setting initial values ===============
;
;
;	store data  
;
imore = 0
if (input_type eq 7) then begin
  a=RASCII(file)
endif else if (input_type eq 3) then begin
  if ((n_elements(input_type) EQ 1) and (handle_info(file))) then begin
    ;print,'PLOTFILE: handle pointer used as input'
    nchild = handle_info(file, /NUM_CHILDREN)
     ;print,'PLOTFILE: ',nchild,' children
     if nchild EQ 0 then handle_value,file,a else begin
       handle=handle_info(file,/FIRST_CHILD)
       handle_value,handle,a
       imore = 1
     endelse
  endif
endif else a=file
;
if keyword_set(frequency) then begin
  frequency = fix(frequency)
  nxx = n_elements(a(0,*)) /frequency
  nyy = n_elements(a(*,0))
  ind = lindgen(nxx) * frequency
  a_old = a
  a = fltarr(nyy,nxx)
  for i=0,n_elements(a(*,0))-1 do begin
   a(i,*) = a_old(i,ind) 
  endfor
endif

;
; return (no plot) if nodata is set
;
if not(keyword_set(histo)) and keyword_set(nodata) then goto,out
;
; prepare the arrays for X and Y
;
if n_elements(a(0,*)) GT 1 then arr1=a(xcol-1,*) else arr1=reform(a)
;arr1=a(xcol-1,*)
if not(keyword_set(histo)) then arr2=a(ycol-1,*)
;
; prepare data for multiple data plot (3D or MESH plots)
;
if imesh then begin
  if keyword_set(MESHCOL) then begin
    tmp = where( a(abs(MESHCOL)-1,*) NE a(abs(MESHCOL)-1,0))
    npointset=tmp(0)
  endif
  if keyword_set(MESHPOINTS) then begin
    npointset = MESHPOINTS
  endif
  if n_elements(npointset) EQ 0 then $
   npointset = n_elements(a(0,*))/ndata
  if npointset EQ 1 then $
    message,'Error in column with multiple data values'
  ndata = n_elements(a(0,*))/npointset
endif
;
; make histograms (if selected)
;
if keyword_set(histo) then begin
  if not(keyword_set(nbins)) then nbins = 25
  if not(keyword_extra(extra,'xrange')) then $
	xrange1=[min(arr1),max(arr1)] else xrange1=extra.xrange
  w = arr1*0 + 1
  binsize = (xrange1(1)-xrange1(0))/float(nbins)
  hy = histogramw(arr1,w,binsize=binsize,min=xrange1(0),max=xrange1(1))
  hx = fltarr(nbins)
  for j=0,nbins-1 do hx(j)=xrange1(0) + binsize/2 + binsize*j
  arr1 = hx
  arr2 = hy
  ;psym = 10
  hisdata = fltarr(2,n_elements(hx))
  hisdata(0,*) = hx  &  hisdata(1,*) = hy
  if keyword_set(nodata) then goto,out
endif
;
; set shifts
;
if keyword_set(shiftval) then begin
  if n_elements(shiftval) NE 2 then begin
    itmp = Dialog_Message(/Error,Dialog_parent=group,$
	'PLOTFILE: ShiftVal has incorrect dimensions. Not used.')
  endif else begin
	arr1 = arr1 + shiftval[0]
	arr2 = arr2 + shiftval[1]
  endelse
endif
;
; plot for single set case
;
IF IMESH NE 1 THEN BEGIN
  if keyword_set(axiscolor) then begin  
    ; plot axis
    if keyword_extra(extra,'color') then begin
      tmp = extra
      tmp.color=axiscolor
      plot,arr1,arr2,xlog=xlog,ylog=ylog,$
	/nodata,_EXTRA=tmp
    endif else plot,arr1,arr2,xlog=xlog,ylog=ylog,$
       /nodata,_EXTRA=extra,color=axiscolor
    ; plot lines
    oplot,arr1,arr2,_EXTRA=extra
  endif else begin
    case kind of
      2: oplot,arr1,arr2,_EXTRA=extra
      else: plot,arr1,arr2,xlog=xlog,ylog=ylog,$
	_EXTRA=extra
    endcase
  endelse
ENDIF
;
; multiset case
;
IF IMESH THEN BEGIN
  sur = fltarr(npointset,ndata)
  for j=0,ndata-1 do begin
    sur(*,j) = arr2(j*npointset:j*npointset+npointset-1)
  endfor
  message,/info,strcompress(ndata,/rem)+' sets of '+$
    strcompress(npointset,/rem)+' points.'
  ; calculate arr3 with labels
  if keyword_set(valndata) then begin
    if (n_elements(valndata) ne ndata) then begin
      message,/info,$
	'Inconsistent number of elements in VALNDATA. Using internal labels'
      arr3 = indgen(ndata)
    endif else begin
      arr3 = valndata
    endelse
  endif 
  if keyword_set(MESHCOL) then begin
    if MESHCOL LT 0 then begin
      arr3 = a(abs(MESHCOL)-1,indgen(ndata)*npointset)
    endif
  endif
  if n_elements(arr3) EQ 0 then arr3 = indgen(ndata)

  if not(keyword_set(meshkind)) then meshkind = 0
  case meshkind of 
    0: begin
        if keyword_set(axiscolor) then begin
          ; plot axis
          if keyword_extra(extra,'color') then begin
            tmp = extra
            tmp.color=axiscolor
            surface,sur,arr1(0:npointset-1),arr3,/save,/nodata,$
             xlog=xlog,ylog=ylog,_EXTRA=tmp
          endif else surface,sur,arr1(0:npointset-1),arr3,/save,/nodata,$
             xlog=xlog,ylog=ylog,_EXTRA=extra,color=axiscolor
        endif else surface,sur,arr1(0:npointset-1),arr3,/save,/nodata,$
             xlog=xlog,ylog=ylog,_EXTRA=extra
        ; plot lines
        for j=0,ndata-1 do begin
          aa = sur(*,j)
          plots,arr1(0:npointset-1),arr3(j),aa,/data,/t3d,_EXTRA=extra
        endfor
      end
    1: begin
        if keyword_set(axiscolor) then begin
          if keyword_extra(extra,'color') then begin
            tmp = extra
            tmp.color=axiscolor
            surface,sur,arr1(0:npointset-1),arr3,$
             xlog=xlog,ylog=ylog,_EXTRA=tmp,/nodata
            surface,sur,arr1(0:npointset-1),arr3,$
             xlog=xlog,ylog=ylog,_EXTRA=extra,/noerase,$
	     xstyle=4,ystyle=4,zstyle=4
	  endif else begin
	     surface,sur,arr1(0:npointset-1),arr3,$
             xlog=xlog,ylog=ylog,_EXTRA=extra,COLOR=axiscolor,/nodata
            surface,sur,arr1(0:npointset-1),arr3,$
             xlog=xlog,ylog=ylog,_EXTRA=extra,/noerase,$
	     xstyle=4,ystyle=4,zstyle=4
	  endelse
        endif else surface,sur,arr1(0:npointset-1),arr3,$
           xlog=xlog,ylog=ylog,_EXTRA=extra 
       end
    2: begin
        contour,sur,arr1(0:npointset-1),arr3,_EXTRA=extra
       end
    3: begin
        shade_surf,sur,arr1(0:npointset-1),arr3,_EXTRA=extra
       end
    else: 
  endcase
ENDIF
;
; multiset when pointer (handle) input type
;
; bug: makes the plot with /nodata is set. Corrected MSR 95-02-13
if keyword_set(imore) then begin
   for i=1,nchild-1 do begin
     handle = handle_info(handle, /SIBLING)
     handle_value,handle,tmpset
     if n_elements(tmpset) gt 1 then oplot,tmpset(xcol-1,*),tmpset(ycol-1,*),$
	_EXTRA=extra
   endfor
endif
;
; gaussian fit
;
if keyword_set(gauss) then begin
  yfit = gauss_fit(arr1,arr2,afit)
  if not(keyword_set(nodata)) then oplot, arr1,yfit,_EXTRA=extra
  print,' Sigma of the gaussian is : ',afit(2)
  print,' FWHM of the gaussian is : ',2.345*afit(2)
endif

;
; legend
;
if keyword_set(legend) then begin
  if strcompress(legend,/rem) NE '' then begin
    command='legend,'+legend
    print,command & tmp=execute(command)
  endif
endif

if keyword_set(pexecute) then $
  if n_elements(pexecute) EQ 1 then tmp = execute(pexecute) else $
	for i=0l,n_elements(pexecute)-1 do tmp = execute(pexecute(i))
;

;
;overplot second axes (if selected)
;
if keyword_set(y2axis) then if (y2axis(0) NE y2axis(1)) then begin
  if keyword_set(axiscolor) then $
    axis,/YAXIS,YRANGE=y2axis,/YSTYLE,YTICKLEN=0.000001,color=axiscolor $
    else axis,/YAXIS,YRANGE=y2axis,/YSTYLE,YTICKLEN=0.000001
endif
if keyword_set(x2axis) then if (x2axis(0) NE x2axis(1)) then begin
  if keyword_set(axiscolor) then $
    axis,/XAXIS,XRANGE=X2axis,/XSTYLE,XTICKLEN=0.000001,color=axiscolor $
    else axis,/XAXIS,XRANGE=X2axis,/XSTYLE,XTICKLEN=0.000001
endif
;
; overplots error columns (if selected)
;
if  keyword_set(perrplot) then $
if perrplot(0) NE perrplot(1) then begin
  if max ([perrplot(0),perrplot(1)]) GT n_elements(a[*,0]) then begin
    message,/info,'Error column number does not exist.'
  endif else begin
    errplot,arr1,a(perrplot(0)-1,*),a(perrplot(1)-1,*)
  endelse
endif
out:
if n_elements(title_old) GT 0 then !p.title=title_old
if n_elements(xtitle_old) GT 0 then !x.title=xtitle_old
if n_elements(ytitle_old) GT 0 then !y.title=ytitle_old
if n_elements(xrange_old) GT 0 then !x.range=xrange_old
if n_elements(yrange_old) GT 0 then !y.range=yrange_old
end
