function make_pix_disk, ipix,val_in,rrad_deg,nside,square=square,triangle=triangle,$
                        bad_value=bad_value,map=map,nsspot=nsspot,$
                        unique=unique, coadd=coadd
;;
;;This function makes disks on a healpix map of nside which subtends a
;;radious of rad and sets the value of it to be val
;;
;; The keywords unique and coadd set the behaviour of how val is put
;; in the output map. 
;; UNIQUE: makes sure that that disks which contains
;; previously visited pixels will not be overwriten, they keep the
;; previously set values.
;; assigned bad value 
;;
;; COADD: when setting a value to a disk, the values are coadded with
;; the previous set values (including bad_value)

if n_params() lt 1 then begin
    print, 'USAGE: make_dist_to_mollview, ipix,val,rad,nside'
    on_error, 0
endif

nn = n_elements(ipix)
if n_params() lt 2 then val = findgen(nn) else val=val_in ;bins=dindgen(50)*16d0+10d0
if n_elements(val) lt nn then val = val[0]+fltarr(nn) 

if n_params() lt 3 then rrad=3d0+fltarr(nn) else rrad=rrad_deg
if n_elements(rrad) lt nn then rrad = rrad[0]+fltarr(nn) 

if n_params() lt 4 then nside=512l

if not keyword_set(nsspot) then nsspot=nside
if n_elements(bad_value) eq 0 then bad_value = 0;float(-1d35)

if not keyword_set(map) then begin
   map=fltarr(long(nside)^2*12l)
   map(*)=bad_value
endif else nside=npix2nside(n_elements(map))


pix2ang_ring, nsspot, ipix, th, phi

for i=0,nn-1 do begin
   rad = rrad[i]*!dtor

    if not keyword_set(square) or not keyword_set(triangle) then begin
        ;print, th(i), phi(i)
       if keyword_set(astro) then begin
          ANG2VEC,  (!pi/2d0-th(i))*!rtod, phi(i)*!rtod, vec,/astro
       endif else begin
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

    if keyword_set(unique) then begin
       icheck = where(map[listpix] gt bad_value/2, count)
       ;;only set these pixels if all pixels in the disc have bad values
       if count eq 0 then map[listpix]=val[i]    
    endif else begin
       if keyword_set(coadd) then begin
          map[listpix]=map[listpix] + val[i] ;;for coadding 
       endif else map[listpix]=val[i] ;;for general     

    endelse

endfor

return, map

end
