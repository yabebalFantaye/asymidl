pro sphere,r,bx,by,bz,upper=upper,lower=lower,npts=npts

; This creates the cartesian coordinates for a sphere
; which can be plotted with PLOTS.

  if not keyword_set(npts) then npts=100

  ; creating the bulge
  ;r = 5.
  nz = 21
  th0 = 0.
  dth = !dpi/(nz-1.)
  if keyword_set(upper) then begin
    th0 = 0.
    dth = 0.5*!dpi/(nz-1.)
  endif
  if keyword_set(lower) then begin
    th0 = !dpi/2.
    dth = 0.5*!dpi/(nz-1.)
  endif

  for i=0,nz-1 do begin
    z = r*cos(i*dth+th0)
    rho = r*sin(i*dth+th0)
    bxarr = rho*sin(findgen(npts))
    byarr = rho*cos(findgen(npts))
    bzarr = fltarr(npts)+z

    if not keyword_set(bx) then begin
      bx = bxarr
      by = byarr
      bz = bzarr
    endif else begin
      bx = [bx,bxarr]
      by = [by,byarr]
      bz = [bz,bzarr]
    endelse
  end

;stop

end
