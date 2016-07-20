function make_disk_to_mollview, th, phi,va_inl,rrad_deg,nside,$
                                square=square,triangle=triangle,bad_value=bad_value,astro=astro
;;
;;This function makes disks on a healpix map of nside which subtends a
;;radious of rad and sets the value of it to be val
;;
;;

if n_params() lt 2 then begin
    print, 'USAGE: make_dist_to_mollview, th, phi,val,rad,nside'
    on_error, 0
endif

nn = n_elements(th)
if n_params() lt 2 then val = findgen(nn) else val=val_in ;bins=dindgen(50)*16d0+10d0
if n_params() lt 3 then rrad=3d0 else rrad=rrad_deg
if n_params() lt 4 then nside=512l

if n_elements(bad_value) eq 0 then bad_value = float(-1d35)

map=fltarr(nside^2*12l)

map(*)=bad_value

rad = rrad*!dtor

radvec = (fltarr(nn)+1)*rad[0]/float(nn)
if n_elements(rad) eq nn then radvec=rad
if (n_elements(rad) ne nn) and (n_elements(rad) ne 1)  then rad=rad[0]
if n_elements(rad) eq 1 then radvec=replicate(rad,nn)

;;if val is constant
if n_elements(val) eq 1 then val=replicate(val,nn)


for i=0,nn-1 do begin   
   rad = radvec[i]*1.4

    if not keyword_set(square) or not keyword_set(triangle) then begin
        ;print, th(i), phi(i)
       if keyword_set(astro) then begin
          if (th(i) lt !pi/9d0) and (th(i) gt -!pi/9d0) then rad=rad*1.2  ;;±20deg increase rad
          ANG2VEC,  (!pi/2d0-th(i))*!rtod, phi(i)*!rtod, vec,/astro
       endif else begin
          if (th(i) gt !pi/9d0) and (th(i) lt 11d0*!pi/18d0) then rad=rad*1.2 ;;±20deg increase rad
          ANG2VEC,  th(i), phi(i), vec
       endelse
        query_disc,nside,vec,rad,listpix
    endif

    if keyword_set(square) then begin
        vec = dblarr(4,3)

        rrad = rad
        if th(i) >0.9*!pi then rrad = sin(rrad)

        thnp = min([th(i) + rrad, !pi])
        thnn = max([th(i) - rrad, 0d0])

        phnp = phi(i) + rad ;/2d0 
        if phnp gt 2d0*!pi then phnp = phnp mod 2d0*!pi

        phnn = phi(i) - rad ;/2d0
        if phnn lt 0 then phnn = 2d0*!pi + phnn



        thn = thnp & phn = phnp
        if keyword_set(astro) then begin
           ANG2VEC,  (!pi/2-thn)*!rtod, phn*!rtod, vec,/astro
        endif else begin
           ANG2VEC,  thn, phn, temp
        endelse

        vec[0,*] = temp

        thn = thnn & phn = phnp
        if keyword_set(astro) then begin
           ANG2VEC,  (!pi/2-thn)*!rtod, phn*!rtod, vec,/astro
        endif else begin
           ANG2VEC,  thn, phn, temp
        endelse
        vec[1,*] = temp

        thn = thnn & phn = phnn
        if keyword_set(astro) then begin
           ANG2VEC,  (!pi/2-thn)*!rtod, phn*!rtod, vec,/astro
        endif else begin
           ANG2VEC,  thn, phn, temp
        endelse
        vec[2,*] = temp

        thn = thnp & phn = phnn
        if keyword_set(astro) then begin
           ANG2VEC,  (!pi/2-thn)*!rtod, phn*!rtod, vec,/astro
        endif else begin
           ANG2VEC,  thn, phn, temp
        endelse
        vec[3,*] = temp


        query_polygon,nside,vec,listpix
    endif

    if keyword_set(triangle) then begin
        vec = dblarr(3,3)

        rrad = rad
        if th(i) >0.9*!pi then rrad = sin(rrad)

        thnp = min([th(i) + rrad, !pi])
        thnn = max([th(i) - rrad, 0d0])

        phnp = phi(i) + rad ;/2d0 
        if phnp gt 2d0*!pi then phnp = phnp mod 2d0*!pi

        phnn = phi(i) - rad ;/2d0
        if phnn lt 0 then phnn = 2d0*!pi + phnn


        thn = thnp & phn = phi(i)+0.01
        if keyword_set(astro) then begin
           ANG2VEC,  (!pi/2-thn)*!rtod, phn*!rtod, vec,/astro
        endif else begin
           ANG2VEC,  thn, phn, temp
        endelse

        vec[0,*] = temp

        thn = thnn & phn = phnp
        if keyword_set(astro) then begin
           ANG2VEC,  (!pi/2-thn)*!rtod, phn*!rtod, vec,/astro
        endif else begin
           ANG2VEC,  thn, phn, temp
        endelse

        vec[1,*] = temp

        thn = thnp & phn = phnn
        if keyword_set(astro) then begin
           ANG2VEC,  (!pi/2-thn)*!rtod, phn*!rtod, vec,/astro
        endif else begin
           ANG2VEC,  thn, phn, temp
        endelse

        vec[2,*] = temp


        query_triangle,nside,vec[0,*],vec[1,*],vec[2,*], listpix
    endif

    map[listpix]=val[i]    
endfor

return, map

end
