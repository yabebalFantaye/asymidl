pro pstyle_plot, x, y, yrange=yrange,yticks=yticks, hline=hline, Charsize=Charsize, _extra=extra,$
                 xcharsize=xcharsize, ycharsize=ycharsize


resolve_routine,'hfi_plot',/COMPILE_FULL_FILE,/either

  ymax=max(y)
  ymin=min(y)
  if not keyword_set(yrange) then yrange=[ymin,ymax]
  


  yrange = smart_yra(yrange[0], yrange[1],ytickv=yval,ytickname=yname,format=format,ny=yticks,/chatty)

  print, 'strlen(ybame): ',strlen(yname)

  ;;print, 'pstyle_plot: yval',yval


  plot, x,y,ytickv=yval,ytickname=yname,yrange=yrange, yticks=n_elements(yval)-1,$
        _extra=extra,Charsize=Charsize , yTickFormat='(A1)',xcharsize=xcharsize, ycharsize=ycharsize
  

;;Position=[0.1, 0.1, 0.9, 0.9]
print, 'pstyle plot format: ',format

  axlabel,yval,orientation=90,/yaxis,format=format,Charsize=yCharsize
  
  ;;horizontal lines
  if keyword_set(hline) then begin
     for ii=0,n_elements(hline)-1 do begin
        oplot, x,y*0+hline[ii],color=0,line=2,thick=3
     endfor
  endif
  


end
