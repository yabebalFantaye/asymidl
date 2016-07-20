function get_xyrange, y, nxy=nxy, _extra=extra

nxy=4
if keyword_set(nxy) then ny=nxy


ymin = float(sigfig(min(y),2))
ymax = (ny-1)*float(sigfig(max(y-ymin)/(ny-1),2))

yra = [ymin,ymax]


end
