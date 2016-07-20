function psym_to_mollview, th, phi,nside,psym,col=col,thick=thick,$
                           fill=fill,bad_value=bad_value,astro=astro, help=help,psize=psize
;;
;;This function makes disks on a healpix map of nside which subtends a
;;radious of rad and sets the value of it to be val
;;

;     PSYM -  The following integer values of PSYM will create the
;             corresponding plot symbols
;     0 - circle
;     1 - downward arrow (upper limit), base of arrow begins at plot
;         value             value
;     2 - upward arrow (lower limt)
;     3 - 5 pointed star
;     4 - triangle
;     5 - upside down triangle
;     6 - left pointing arrow
;     7 - right pointing arrow
;     8 - square
;;
  if keyword_set(help) then begin
     print, 'USAGE: '
     print, 'psym_to_mollview, th, phi,nside,psym,col=col,bad_value=bad_value,astro=astro, help=help'
     On_error,0
  endif

  if n_params() lt 2 then begin
     print, 'USAGE: '
     print, 'psym_to_mollview, th, phi,nside,psym,col=col,bad_value=bad_value,astro=astro, help=help'
     On_error,0
  endif

  nn = n_elements(th)
  if n_params() lt 4 then psym=2
  if n_params() lt 3 then nside=512l
  if not keyword_set(col) then col = findgen(nn) 
  val = fltarr(nn) 
  val[*] = col


  if N_elements(psym) LT 1 then begin
     print, 'psym_to_mollview, th, phi,nside,psym,col=col,bad_value=bad_value,astro=astro, help=help'
     print,'  PSYM values 0 - circle, 1 - down arrow, 2 - up arrow, 3 - star'
     print,'       4 - triangle, 5 - upside down triangle, 6 - left arrow'
     print,'       7 - right arrow, 8 - square'
     On_error,0
  endif

  map=fltarr(nside^2*12l)
  bad_value = float(-1d35)
  map(*)=bad_value

  if not keyword_set(psize) then psize = 0.001 
  if ~keyword_set(FILL) then fill = 0
  if ~keyword_set(thick) then thick=1


  for i=0,nn-1 do begin

     if keyword_set(astro) then begin
        ANG2VEC,  (!pi/2d0-th(i))*!rtod, phi(i)*!rtod, vec,/astro
     endif else begin
        ANG2VEC,  th(i), phi(i), vec
     endelse


     case psym of
        0:  begin               ;Circle
           ang = 2*!PI*findgen(49)/48. ;Get position every 5 deg
           xarr = psize*cos(ang)  &  yarr = psize*sin(ang)
        end
        1:  begin               ;Down arrow
           xarr = [0,0,.5,0,-.5]*psize
           yarr = [0,-2,-1.4,-2,-1.4]*psize
           fill = 0
        end
        2:  begin               ;Up arrow
           xarr = [0,0,.5,0,-.5]*psize
           yarr = [0,2,1.4,2,1.4]*psize
           fill = 0
        end
        3:  begin                                ;Star
           ang = (360. / 10 * findgen(21) + 90) / !RADEG ;star angles every 36 deg
           r = ang*0
           r[2*indgen(12)] = 1.
           cp5 = cos(!pi/10.)
           r1 = 2. * cp5 - 1. / cp5
           r[2*indgen(10)+1] = r1
           r = r * psize / sqrt(!pi/4.) * 2. / (1.+r1)
           xarr = r * cos(ang)   &   yarr = r * sin(ang)
        end
        4:  begin               ;Triangle
           xarr = [-1,0,1,-1]*psize
           yarr = [-1,1,-1,-1]*psize
        end
        5:  begin               ;Upside down triangle
           xarr = [-1, 0, 1, -1]*psize
           yarr = [ 1,-1, 1, 1]*psize
        end
        6:  begin               ;Left pointing arrow
           yarr = [0, 0, 0.5, 0, -.5]*psize
           xarr = [0,-2,-1.4,-2,-1.4]*psize
           fill = 0
        end
        7:  begin               ;Left pointing arrow
           yarr = [ 0, 0, 0.5, 0, -.5] * psize
           xarr = [ 0, 2, 1.4, 2, 1.4] * psize
           fill = 0
        end
        8:  begin               ;Square
           xarr = [-1,-1,1, 1,-1] * psize
           yarr = [-1, 1,1,-1,-1] * psize
        end
        else: message,'Unknown plotting symbol value of '+strtrim(psym,2)
     endcase

     xarr = xarr[0:n_elements(xarr)-2]
     yarr = yarr[0:n_elements(yarr)-2]


     vlist = dblarr(n_elements(xarr),3)
     vlist[*,0] = vec[0]+xarr 
     vlist[*,1] = vec[1]+yarr
     vlist[*,2] = vec[2]+0.*xarr 
 
     query_polygon,nside,vlist,listpix

     map(listpix)=val(i)    
  endfor

  return, map

end
