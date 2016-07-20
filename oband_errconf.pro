pro oband_errconf,x,y1,y2,cvec=cvec,mat=mat,pvalue=pvalue,_extra=extra


  if keyword_set(mat) then simconflevel,mat,y1,y2,mmat=mmat,pvalue=pvalue

  sy=size(y1)

  if n_params() lt 1 then x = lindgen(sy[1])

  if sy[0] gt 1 then begin

     ny = n_elements(y1[0,*])

     if not keyword_set(cvec) then begin
        loadct,54,ncolors=ny+6
        cvec=indgen(ny+6)
     endif

     for i=0,ny-1 do begin
        oband,x,reform(y1[*,i]),reform(y2[*,i]),color=cvec[i+2],_extra=extra
     endfor

  endif else begin

     if not keyword_set(cvec) then cvec=19
     oband,x,y1,y2,color=cvec,_extra=extra

  endelse

end
