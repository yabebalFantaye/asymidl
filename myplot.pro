pro myplot, x,y,xra=xra,yra=yra,nx=nx, ny=ny, _extra=extra

if not keyword_set(nx) then nx=4
if not keyword_set(ny) then ny=4
if not keyword_set(xra) then xra=[min(x),max(x)]
if not keyword_set(yra) then yra=[min(y),max(y)]

if not keyword_set(ny) then ny=3
if not keyword_set(nmy) then nmy = 1
if not keyword_set(nx) then nx=4
if not keyword_set(nmx) then nmx = 1

xmax = 100*xra[1]

count_start=3
count=count_start
if xmax lt 1 then begin
   while xmax gt 1 do begin
      xmax = xra[1]*10^count
      count=count+1
   endwhile
endif


;xformat='()'

xlblv = linspace(xra[0],xra[1],nx)
ylblv = linspace(yra[0],yra[1],ny)

if keyword_set(xlog) then begin
   xtickv = LogLevels(xra)
   xtickname=string(xtickv)
   nx = N_Elements(xtickv)-1
endif else begin
   xtickv = intarr(nx+1)
   xtickname = strarr(nx+1)
   for jj=0,nx do begin
      xval = xra[0] + jj*(xra[1]-xra[0])/nx
      xtickv[jj] = xval
      xtickname[jj] = strn(xval)
   endfor
endelse

if keyword_set(ylog) then begin
   ytickv = LogLevels(yra)
   ytickname=string(ticksy)
   ny = N_Elements(ticksy)-1
endif else begin
   ytickv = intarr(ny+1)
   ytickname = strarr(ny+1)
   for jj=0,ny do begin
      yval = yra[0] + jj*(yra[1]-yra[0])/ny
      ytickv[jj] = yval
      ytickname[jj] = strn(yval)
   endfor   
endelse


plot,x,y,xr=xra,yr=yra, $ ; do the plot
     xtickf='(a1)',ytickf='(a1)',_extra=extra



axlabel,xlblv, /xaxis           ; plot the labels
axlabel,ylblv, /yaxis

end
