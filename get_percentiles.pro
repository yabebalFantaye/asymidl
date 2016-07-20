;-------------------------------------------------------------
;+
; NAME:
;        PERCENTILES
;
; PURPOSE:
;        compute percentiles of a data array
;
; CATEGORY:
;        statistical function
;
; CALLING SEQUENCE:
;        Y = PERCENTILES(DATA [,VALUE=value-array])
;
; INPUTS:
;        DATA --> the vector containing the data
;
; KEYWORD PARAMETERS:
;        VALUE --> compute specified percentiles
;        default is a standard set of min, 25%, median (=50%), 75%, and max 
;        which can be used for box- and whisker plots.
;        The values in the VALUE array must lie between 0. and 1. !
;
; OUTPUTS:
;        The function returns an array with the percentile values or
;        -1 if no data was passed or value contains invalid numbers.
;
; SUBROUTINES:
;
; REQUIREMENTS:
;
; NOTES:
;
; EXAMPLE:
;      x = (findgen(31)-15.)*0.2     ; create sample data
;      y = exp(-x^2)/3.14159         ; compute some Gauss distribution
;      p = percentiles(y,value=[0.05,0.1,0.9,0.95])
;      print,p
;
;      IDL prints :  3.92826e-05  0.000125309     0.305829     0.318310
 
;
; MODIFICATION HISTORY:
;        mgs, 03 Aug 1997: VERSION 1.00
;        mgs, 20 Feb 1998: - improved speed and memory usage 
;                (after tip from Stein Vidar on newsgroup)
;        antunes, 16 May 2007: changed 'fix' to 'long' so this works
;                on data larger than 128x128
;
;-
; Copyright (C) 1997, Martin Schultz, Harvard University
; This software is provided as is without any warranty
; whatsoever. It may be freely used, copied or distributed
; for non-commercial purposes. This copyright notice must be
; kept with any copy of this software. If this software shall
; be used commercially or sold as part of a larger package,
; please contact the author to arrange payment.
; Bugs and comments should be directed to mgs@io.harvard.edu
; with subject "IDL routine percentiles"
;-------------------------------------------------------------


function get_percentiles,data,dvalue,ix=ix,iy=iy,log=log,sigma=sigma
  
  result = -1
  n = n_elements(data)
  if (n le 0) then return,result ; error : data not defined

;;get the shape and n_elements of data
  sm=size(data) 


;;sort the dvalue to which percentiles are to be found
  dvalue_temp = dvalue
  ;; if n_elements(dvalue) gt 1 then begin
  ;;    iy=sort(dvalue)
  ;;    dvalue_temp=dvalue[iy]
  ;; endif 



;;when data is just an array
  if (sm(0) EQ 1) then begin

     ;;print, 'case: data is vector and dvalue is scalar/vector'
     ;;create a temporary copy of the data and sort
     ix = sort(data)
     data_temp = data[ix]

     iz = value_locate(data_temp,dvalue_temp)
     if total(where(iz eq -1)) gt -1 then iz[where(iz eq -1)]=0 ;; values outside data are set as the extremem percentiles
     value=float(iz+1)/float(n_elements(data_temp)+1)

  endif 

;;when data is a matrix e.g. data[npix,nsim]
  if (sm(0) EQ 2) then begin


     ;;dvalue must be scalar or array of length sm[1]
     if n_elements(dvalue) eq 1 then dvalue=replicate(dvalue, sm[1])
     if n_elements(dvalue) ne sm[1] then on_error, 2 

     ;;print, 'case: data is matrix and dvalue is vector'

     ;;create a temporary copy of the data and sort
     data_temp = msort(data,/sort_now) ;;for mat[ncol,nrow], the row is sorted 

     value=dvalue_temp*0d0

     for ii=0,sm[1]-1 do begin

        iz = value_locate(reform(data_temp[ii,*]),dvalue_temp[ii])
        if iz eq -1 then iz=0 ;; values outside data are set as the extremem percentiles


        value[ii]=float(iz+1)/float(sm[2]+1)
     endfor

  endif



  if keyword_set(log) then begin
     value=-alog10(value)
  endif

  if keyword_set(sigma) then begin
     value=value*100
  endif
  
  return,value
end
 
