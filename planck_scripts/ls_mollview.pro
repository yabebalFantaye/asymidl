;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2012  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.jpl.nasa.gov
;
; -----------------------------------------------------------------------------
pro LS_oplot_graticule, graticule, eul_mat, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic, flip = flip, _extra = oplot_kw, half_sky=half_sky, coordsys=coordsys, charsize=charsize, reso_rad=reso_rad, GRMIN=GRMIN, GRMAX=GRMAX
;+
; NAME:
;       OPLOT_GRATICULE
;
; PURPOSE:
;       overplots graticule (ie, spherical coordinates grid) on top
;       of existing map
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;       oplot_graticule, graticule, eul_mat, $
;          [projection=,mollweide=,gnomic=,cartesian=, $
;           orthographic=,flip=,half_sky=,coordsys=, reso_rad=, + all oplot keywords]
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;
;
; COMMON BLOCKS:
;
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;      calls oplot_sphere
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;         Feb 2003, corrected bug occuring when combining astrophysical coordinates
;         changes with arbitrary rotation
;         Jan. 2013, L. Spencer: added GRMIN and GRMAX keywords to prevent printing graticule labels on edges of map (i.e. -180, +180, -90, etc.)
;-
IF N_ELEMENTS(GRMIN) EQ 1 THEN GRMIN = DBLARR(2) + GRMIN
IF N_ELEMENTS(GRMAX) EQ 1 THEN GRMAX = DBLARR(2) + GRMAX
;
identify_projection, projtype, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian

if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)

do_ig = (n_elements(coordsys) eq 2) 

; define default graticule spacing, in Degrees
if (projtype eq 2) then begin ; gnom
    dlong = 5. 
    dlat  = 5. 
    gratmin = 0.
endif
if (projtype eq 1) then begin ; moll
    dlong = 45.
    dlat  = 45.
    gratmin = 10.
endif
if (projtype eq 3) then begin ; cart
    dlong = 45. 
    dlat  = 45. 
    gratmin = 2.
endif
if (projtype eq 4) then begin ; ortho
    dlong = 45.
    dlat  = 45.
    gratmin = 10.
endif

; read in user defined grid spacings
if (n_elements(graticule) eq 2 and min(graticule) ge gratmin) then begin
    dlong = float(graticule[0]) & dlat = float(graticule[1])
endif else begin
    if (graticule[0] gt gratmin) then begin 
        dlong = float(graticule[0]) & dlat = dlong
    endif
endelse

fsgrat = (defined(reso_rad)) ? !dtor/(reso_rad+1.d-11) : 1. ; number of pixels / degree
fsgrat = long(fsgrat) > 1 < 5
; define variables
epsilon = 1.d-5
nmerid = fix(181./dlong)
nparal = fix(90./dlat)-1
nv = (projtype eq 2 || projtype eq 3) ? 721*fsgrat : 361 ; more points for partial projections
vector = DINDGEN(nv)/(nv-1.d0) * (1.d0-2.d0*epsilon) + epsilon ; from epsilon to 1-epsilon
bounds = [[-nmerid,nmerid-1],[-nparal,nparal]]

do_rot = (n_elements(eul_mat) eq 9)
form = '(i4)'
if (abs(round(dlong)-dlong) gt 1.e-2 || abs(round(dlat)-dlat) gt 1.e-2) then form = '(f6.1)'

case projtype of
2: begin  ; gnomic : straightforward
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
;             if (jg eq 0) then ang2vec, vector*!pi,        replicate(i*dlong*!DtoR,nv), vv ; meridians
;             if (jg eq 1) then ang2vec, replicate((90.-i*dlat)*!DtoR,nv), vector*(2.*!pi), vv ; parallels
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            k = where(vv(*,0) gt 0, nk)
            if (nk gt 0) then begin
                u = vv(k,1)/vv(k,0)
                v = vv(k,2)/vv(k,0)
                good = where(abs(u) lt !x.crange[1]*1.1 and abs(v) lt !y.crange[1]*1.1 ,ng)
                ; reorder points to have one continuous segment across the plot
                bad = where(good-shift(good,1) ne 1, nbad)
                if (nbad gt 1) then good = shift(good, bad[1])
;                oplot, flipconv * u, v, _extra = oplot_kw
               if (ng gt 1) then oplot_sphere, flipconv *u[good], v[good], _extra = oplot_kw, linelabel=linelabel,/flush, charsize=charsize
            endif
        endfor
    endfor
end
1: begin  ; mollweide : deal with boundaries
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                IF N_ELEMENTS(GRMIN) GT 0 THEN BEGIN
                   IF mylong LT GRMIN[0] THEN linelabel = '' ELSE linelabel = strtrim(string(mylong,form=form),2);+'!Xo'
                   IF N_ELEMENTS(GRMAX) GT 0 THEN IF mylong GT GRMAX[0] THEN linelabel = ''
                ENDIF ELSE BEGIN
                   linelabel = strtrim(string(mylong,form=form),2);+'!Xo'
                   IF N_ELEMENTS(GRMAX) GT 0 THEN BEGIN
                     IF mylong GT GRMAX[0] THEN linelabel = ''
                   ENDIF ; ELSE linelabel = strtrim(string(mylong,form=form),2)
                ENDELSE
                ;linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                IF N_ELEMENTS(GRMIN) GT 0 THEN BEGIN
                   IF mylat LT GRMIN[1] THEN linelabel = '' ELSE linelabel = strtrim(string(mylat,form=form),2);+'!Xo'
                   IF N_ELEMENTS(GRMAX) GT 0 THEN IF mylat GT GRMAX[1] THEN linelabel = ''
                ENDIF ELSE BEGIN
                   linelabel = '!X'+strtrim(string(mylat,form=form),2);+'!Xo'
                   IF N_ELEMENTS(GRMAX) GT 0 THEN BEGIN
                     IF mylat GT GRMAX[1] THEN linelabel = ''
                   ENDIF ; ELSE linelabel = strtrim(string(mylat,form=form),2)
                ENDELSE
                ;linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            vec2moll, vv, u, v
            ;oplot_sphere, -flipconv * u, v, _extra = oplot_kw, linelabel=linelabel, charsize=charsize
            device, /HELVETICA , FONT_size=8
            oplot_sphere,  u, v, _extra = oplot_kw, linelabel=linelabel, charsize=charsize ;, ;FLUSH
            

;;            oplot_sphere, flipconv * u, v, _extra = oplot_kw
        endfor
    endfor
end
4: begin  ; orthographic : deal with boundaries
    if keyword_set(half_sky) then begin
        nd = 1 ; number of half-sky disc
        c0 = 0
    endif else begin
        nd = 2
        c0 = 1
    endelse
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            for sign = 1,1-nd,-2 do begin ; either (1,-1) or (1)
                k = where(vv[*,0]*sign ge 0, nk)
                if (nk gt 0) then begin
                    u = vv[k,1]
                    v = vv[k,2]
                    ;oplot_sphere, flipconv*(u+c0)*sign, v, _extra = oplot_kw
                    oplot_sphere,  flipconv*(u+c0)*sign, v, _extra = oplot_kw, linelabel=linelabel, charsize=charsize
                endif ; nk>0
            endfor ; loop on sign
        endfor
    endfor
end
3: begin  ; cartesian : straightforward
    for jg=0,1 do begin
        for i = bounds[0,jg],bounds[1,jg] do begin
;             if (jg eq 0) then ang2vec, vector*!pi,        replicate(i*dlong*!DtoR,nv), vv ; meridians
;             if (jg eq 1) then ang2vec, replicate((90.-i*dlat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            if (jg eq 0) then begin
                mylong = i*dlong ; longitude in Deg
                linelabel = strtrim(string(mylong,form=form),2)
                ang2vec, vector*!pi,        replicate(mylong*!DtoR,nv), vv ; meridians
            endif
            if (jg eq 1) then begin
                mylat = i*dlat ; latitude in Deg
                linelabel = strtrim(string(mylat,form=form),2)
                ang2vec, replicate((90.-mylat)*!DtoR,nv), vector*2*!pi, vv ; parallels
            endif
            if (do_ig ) then vv = rotate_coord(vv,in=coordsys[0],out=coordsys[1])
            if (do_rot) then vv = vv # transpose(eul_mat)

            phi = atan(vv[*,1],vv[*,0]) ; in [0,2pi]
      theta = asin(vv[*,2])       ; in [0,pi]
            ; OPLOT,-flipconv*phi,theta, _extra = oplot_kw
;            oplot_sphere, flipconv * phi, theta, _extra = oplot_kw
;            oplot_sphere, flipconv* phi, theta, _extra = oplot_kw,
;            linelabel=linelabel, charsize=charsize
            good = where(abs(phi) lt !x.crange[1]*1.1 and abs(theta) lt !y.crange[1]*1.1 ,ng)
                                ; reorder points to have one continuous segment across the plot
            bad = where(good-shift(good,1) ne 1, nbad)
            if (nbad gt 1) then good = shift(good, bad[1])
            if (ng gt 1) then oplot_sphere, flipconv *phi[good], theta[good], _extra = oplot_kw, linelabel=linelabel,charsize=charsize
        endfor
    endfor
end
endcase


;-----------------------
;     FOR i=-nmerid,nmerid-1 DO begin  ; meridians
;         ang2vec, vector*!pi, replicate(i*dlong*!DtoR,n_elements(vector)), vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         k = where(vv(*,0) gt 0, nk)
;         if (nk gt 0) then begin
;             u = vv(k,1)/vv(k,0)
;             v = vv(k,2)/vv(k,0)
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif
;     endfor
;     FOR i=-nparal,nparal DO begin  ; parallels
;         ang2vec, replicate((90.-i*dlat)*!DtoR,n_elements(vector)), vector*2*!pi, vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         k = where(vv(*,0) gt 0, nk)
;         if (nk gt 0) then begin
;             u = vv(k,1)/vv(k,0)
;             v = vv(k,2)/vv(k,0)
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif
;     endfor
; endif

; if (do_moll) then begin
;     ; mollweide : deal with boundaries
;     FOR i=-nmerid,nmerid-1 DO begin  ; meridians
;         ang2vec, vector*!pi, replicate(i*dlong*!DtoR,nv), vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         vec2moll, vv, u, v
;         bad = where(abs(u-shift(u,1)) gt .1, nbad)
;         if (nbad eq 0) then begin
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif else begin
;             bad = [0,bad,n_elements(u)-1]
;             for j=0,nbad do begin
;                 if (bad[j+1] gt bad[j]) then $
;                   oplot, flipconv * u[bad[j]:bad[j+1]-1], v[bad[j]:bad[j+1]-1], color=!p.color
;             endfor
;         endelse
;     endfor
;     FOR i=-nparal,nparal DO begin  ; parallels
;         ang2vec, replicate((90.-i*dlat)*!DtoR,nv), vector*2*!pi, vv
;         if (do_rot) then vv = vv # transpose(eul_mat)

;         vec2moll, vv, u, v
;         bad = where(abs(u-shift(u,1)) gt .1, nbad)
;         if (nbad eq 0) then begin
;             OPLOT, flipconv * u, v, COLOR = !P.COLOR
;         endif else begin
;             bad = [0,bad,n_elements(u)-1]
;             for j=0,nbad do begin
;                 if (bad[j+1] gt bad[j]) then $
;                   oplot, flipconv * u[bad[j]:bad[j+1]-1], v[bad[j]:bad[j+1]-1], color=!p.color
;             endfor
;         endelse
;     endfor
; endif

return
end
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;

; -----------------------------------------------------------------------------
pro LS_data2gnom, data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
               color, Tmax, Tmin, color_bar, dx, planvec, vector_scale, $
               PXSIZE=pxsize, PYSIZE=pysize, ROT=rot_ang, LOG=log, HIST_EQUAL=hist_equal, $
               MAX=max_set, MIN=min_set, $
               RESO_ARCMIN=reso_arcmin, FITS = fits, $
               FLIP=flip, DATA_plot = data_plot, $
               POLARIZATION=polarization, SILENT=silent, PIXEL_LIST=pixel_list, ASINH=asinh, $
               TRUECOLORS=truecolors, DATA_TC=data_tc, MAP_OUT=map_out

;+
;==============================================================================================
;     DATA2GNOM
;
;     turns a Healpix or Quad-cube map into in Gnomonic rectangular map
;
;     DATA2GNOM,  data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat,
;          color, Tmax, Tmin, color_bar, dx, planvec, vector_scale,
;          pxsize=, pysize=, rot=, log=, hist_equal=, max=, min=,
;          reso_arcmin=, fits=, flip=, data_plot=, polarization=, silent=,
;          pixel_list=, TRUECOLORS=, DATA_TC=, MAP_OUT=
;
; IN :
;      data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat
; OUT :
;      color, Tmax, Tmin, color_bar, dx, planvec, vector_scale
; KEYWORDS
;      Pxsize, Pysize, Rot, Log, Hist_equal, Max, Min, Reso_arcmin,
;      Fits, flip, data_plot, polarization, pixel_list, asinh, map_out
;
;  called by gnomview
;
;  HISTORY; Feb 2005: added small_file to avoid pol direction variation within pixels
; Sep 2007: added /silent
; April 2008: added pixel_list=
; July 2008: added asinh
; May 2009: can deal with maps without any valid pixel
;==============================================================================================
;-

do_true = keyword_set(truecolors) 
truetype = do_true ? truecolors : 0
proj_small = 'gnomic'
du_dv = 1.    ; aspect ratio
fudge = 1.00  ; 
if keyword_set(flip) then flipconv=1 else flipconv = -1  ; longitude increase leftward by default (astro convention)
if undefined(polarization) then polarization=0
do_polamplitude = (polarization[0] eq 1)
do_poldirection = (polarization[0] eq 2)
do_polvector    = (polarization[0] eq 3)

!P.BACKGROUND = 1               ; white background
!P.COLOR = 0                    ; black foreground

mode_col = keyword_set(hist_equal)
mode_col = mode_col + 2*keyword_set(log) + 4*keyword_set(asinh)

obs_npix = n_elements(data)
npix_full = (pix_type eq 'Q') ? 6*(4L)^(pix_param-1) : nside2npix(pix_param)
bad_data= !healpix.bad_value

if (do_poldirection or do_polvector) then begin
    ; compute new position of pixelisation North Pole in the plot coordinates
    north_pole = [0.,0.,1.]
    if (do_conv) then north_pole = SKYCONV(north_pole, inco= coord_in, outco=coord_out)
    if (do_rot) then north_pole = north_pole # transpose(eul_mat)
endif
; -------------------------------------------------------------
; create the rectangular window
; -------------------------------------------------------------
if defined(pxsize) then xsize = pxsize*1L else xsize = 500L
if defined(pysize) then ysize = pysize*1L else ysize = xsize
if defined(reso_arcmin) then resgrid = reso_arcmin/60. else resgrid = 1.5/60.
dx      = resgrid * !DtoR
zsize = (do_true) ? 3 : 1
N_uv = xsize*ysize
indlist = (n_elements(pixel_list) eq n_elements(data[*,0]))
small_file = ((!pi*4./dx^2 GT npix_full && do_poldirection))


if (~keyword_set(silent)) then begin
    print,'Input map  :  ',3600.*6.d0/sqrt(!dpi*npix_full),' arcmin / pixel ',form='(a,f8.3,a)'
    print,'gnomonic map :',resgrid*60.,' arcmin / pixel ',xsize,'*',ysize,form='(a,f8.3,a,i4,a,i4)'
endif

if (small_file) then begin
    ; file smaller than final map, make costly operation on the file
    ; initial data is destroyed and replaced by color
    if (do_poldirection or do_polvector) then begin
        phi = 0.
        if (do_rot or do_conv) then begin
            ; position of each map pixel after rotation and coordinate changes
            if (indlist) then begin
                id_pix = pixel_list
            endif else begin
                id_pix = lindgen(npix_full)
            endelse
            case pix_type of
                'R' : PIX2VEC_RING, pix_param, id_pix, vector ; Healpix ring
                'N' : PIX2VEC_NEST, pix_param, id_pix, vector; Healpix nest
                'Q' : vector = PIX2UV(pix_param, id_pix) ; QuadCube (COBE cgis software)
                else : print,'error on pix_type'
            endcase
            id_pix = 0
            if (do_conv) then vector = SKYCONV(vector, inco= coord_in, outco=coord_out)
            if (do_rot) then vector = vector # transpose(eul_mat)
            ; compute rotation of local coordinates around each vector
            tmp_sin = north_pole[1] * vector[*,0] - north_pole[0] * vector[*,1]
            tmp_cos = north_pole[2] - vector[*,2] * (north_pole[0:2] ## vector)
            if (flipconv lt 0) then tmp_cos = flipconv * tmp_cos
            phi = ATAN(tmp_sin, tmp_cos) ; angle in radians
            tmp_sin = 0. & tmp_cos = 0 & vector = 0.
        endif
        data_plot = data
        if (do_poldirection) then begin
            data = (data - phi + 4*!PI) MOD (2*!PI) ; angle
            min_set = 0. & max_set = 2*!pi
        endif
        if (do_polvector) then begin
            pol_data[*,1] = (pol_data[*,1] - phi + 4*!PI) MOD (2*!PI) ; angle is rotated
        endif
    endif else begin ; temperature only or polarisation amplitude only
        data_plot = data
    endelse
    ; color observed pixels
    if (do_true) then begin
        if (truetype eq 2) then begin
            for i=0,2 do begin
                find_min_max_valid, data_tc[*,i], mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
                data_tc[0,i] = COLOR_MAP(data_tc[*,i], mindata, maxdata, Obs, $
                                    color_bar = color_bar, mode=mode_col, silent=silent )
            endfor
        endif else begin
            find_min_max_valid, data_tc, mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
            data_tc = COLOR_MAP(data_tc, mindata, maxdata, Obs, $
                                color_bar = color_bar, mode=mode_col, $
                                minset = min_set, maxset = max_set, silent=silent )
        endelse
    endif else begin
        find_min_max_valid, data, mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
        data    = COLOR_MAP(data, mindata, maxdata, Obs, $
                         color_bar = color_bar, mode=mode_col, $
                         minset = min_set, maxset = max_set, silent=silent )
    endelse
    if (do_polvector) then begin ; rescale polarisation vector in each valid pixel
        pol_data[0,0] = vector_map(pol_data[*,0], Obs, vector_scale = vector_scale)
    endif
    if defined(Obs) then Obs = 0
    Tmin = mindata & Tmax = maxdata
    color = MAKE_ARRAY(/BYTE, xsize, ysize, zsize, Value = !P.BACKGROUND) ; white
    grid = FLTARR(xsize, ysize)
endif else begin ; large
    grid = FLTARR(xsize, ysize, zsize)
endelse
if do_polvector then planvec = MAKE_ARRAY(/FLOAT,xsize,ysize, 2, Value = bad_data) 
; -------------------------------------------------------------
; makes the projection around the chosen contact point
; -------------------------------------------------------------
; position on the planar grid  (1,u,v)
x0 = +1.
xll= 0 & xur =  xsize-1
yll= 0 & yur =  ysize-1
xc = 0.5*(xll+xur)  ; & deltax = (xur - xc)
yc = 0.5*(yll+yur)  ; & deltay = (yur - yc)

xsize=floor(xsize)
ysize=floor(ysize)

yband = LONG(5.e5 / FLOAT(xsize))
for ystart = 0, ysize - 1, yband do begin 
    yend   = (ystart + yband - 1) < (ysize - 1)
    nband = floor(yend - ystart + 1l)
    npb = xsize * nband
    u = flipconv*(FINDGEN(xsize) - xc)# REPLICATE(dx,nband)   ; minus sign = astro convention
    v =           REPLICATE(dx,xsize) # (FINDGEN(nband) + ystart - yc)


    x = replicate(x0, npb)
    help, u
    help, npb
    help, nband
    help, xsize

    npb=n_elements(u)

     vector = [[x],[reform(u,npb,/over)],[reform(v,npb,/over)]] ; non normalised vector
     ; --------------------------------
     ; deal with polarisation direction
     ; --------------------------------
     if (do_poldirection or do_polvector) then begin
         phi = 0.
         if (do_rot or do_conv) then begin
             vector = vector / (sqrt(total(vector^2, 2))#replicate(1,3)) ; normalize vector
             ; compute rotation of local coordinates around each vector
             tmp_sin = north_pole[1] * vector[*,0] - north_pole[0] * vector[*,1]
             tmp_cos = north_pole[2] - vector[*,2] * (north_pole[0:2] ## vector)
             if (flipconv lt 0) then tmp_cos = flipconv * tmp_cos
             phi = ATAN(tmp_sin, tmp_cos) ; angle in radians
             tmp_sin = 0. & tmp_cos = 0
         endif
     endif
     ; ---------
     ; rotation
     ; ---------
     if (do_rot) then vector = vector # eul_mat
     if (do_conv) then vector = SKYCONV(vector, inco = coord_out, outco =  coord_in)
           ; we go from the final Gnomonic map (system coord_out) to
           ; the original one (system coord_in)
     ; -------------------------------------------------------------
     ; converts the position on the sphere into pixel number
     ; and project the corresponding data value on the map
     ; -------------------------------------------------------------
     case pix_type of
         'R' : VEC2PIX_RING, pix_param, vector, id_pix ; Healpix ring
         'N' : VEC2PIX_NEST, pix_param, vector, id_pix ; Healpix nest
         'Q' : id_pix = UV2PIX(vector, pix_param) ; QuadCube (COBE cgis software)
         else : print,'error on pix_type'
     endcase
     if (small_file) then begin ; (data and data_pol are already rescaled and color coded)
         if (do_true) then begin
             for i=0,zsize-1 do color[ystart*xsize+i*n_uv] = data[id_pix,i]
             grid[ystart*xsize]  = data_plot[id_pix] ; unaltered data
         endif else begin
             color[ystart*xsize] = data[id_pix]
             grid[ystart*xsize]  = data_plot[id_pix] ; unaltered data
             if (do_polvector) then begin
                 planvec[ystart*xsize]       = pol_data[id_pix,0] ; amplitude
                 planvec[ystart*xsize+n_uv]  = pol_data[id_pix,1] ; direction
             endif
         endelse
     endif else begin            ; (large file : do the projection first)
         if (do_true) then begin
             for i=0,zsize-1 do grid[ystart*xsize+i*n_uv] = data_tc[id_pix,i]
         endif else begin
             if (do_poldirection) then begin
                 grid[ystart*xsize] = (data[id_pix] - phi + 4*!PI) MOD (2*!PI) ; in 0,2pi
             endif else if (do_polvector) then begin
                 grid[ystart*xsize]         = data[id_pix]
                 planvec[ystart*xsize]      = pol_data[id_pix,0]
                 planvec[ystart*xsize+n_uv] = (pol_data[id_pix,1] - phi + 4*!PI) MOD (2*!PI) ; angle
             endif else begin
 ;;;            grid[ystart*xsize] = data[id_pix]
                 grid[ystart*xsize] = sample_sparse_array(data, id_pix, in_pix=pixel_list, default= !healpix.bad_value)
             endelse
         endelse
     endelse
 endfor
 u = 0 & v = 0 & x = 0 & vector = 0

 ; -------------------------------------------------------------
 ; Test for unobserved pixels
 ; -------------------------------------------------------------
 if (small_file) then begin
     data = 0 & pol_data = 0
 endif else begin
     data_plot = temporary(data)
     pol_data = 0
     find_min_max_valid, grid, mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
 endelse

 ;-----------------------------------
 ; export in FITS and as an array the original gnomic map before alteration
 ;-----------------------------------

 ; grid -> IDL array
 if arg_present(map_out) then map_out = proj2map_out(grid, bad_data=bad_data)

 ; grid -> FITS file
 if keyword_set(fits) then begin 
     proj2fits, grid, fits, $
                projection = 'GNOM', flip=flip, $
                rot = rot_ang, coord=coord_out, reso = resgrid*60., unit = sunits, min=mindata, max = maxdata
 endif

 ; -------------------------------------------------------------
 ; set min and max and computes the color scaling
 ; -------------------------------------------------------------
 if (small_file) then begin

 endif else begin
     if (do_poldirection) then begin
         min_set = 0.
         max_set = 2*!pi
     endif
     if (truetype eq 2) then begin
         ; truecolors=2 map each field to its color independently
         color = bytarr(xsize,ysize,zsize)
         for i=0,zsize-1 do begin
             find_min_max_valid, grid[*,*,i], mindata, maxdata, valid=Obs, bad_data = 0.9 * bad_data
             color[0,0,i] = COLOR_MAP(grid[*,*,i], mindata, maxdata, Obs, $
                           color_bar = color_bar, mode=mode_col, silent=silent)
         endfor
     endif else begin
         ; same for truecolors=1 and false colors:
         color = COLOR_MAP(grid, mindata, maxdata, Obs, $
                           color_bar = color_bar, mode=mode_col, $
                           minset = min_set, maxset = max_set, silent=silent)
     endelse

     if (do_polvector) then begin ; rescale polarisation vector in each valid pixel
         planvec[*,*,0] = vector_map(planvec[*,*,0], Obs, vector_scale = vector_scale)
     endif
     Obs = 0
     grid = 0
     Tmin = mindata & Tmax = maxdata
 endelse

 return
 end
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ; -----------------------------------------------------------------------------
 ;
 ;  Copyright (C) 1997-2012  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
 ;
 ;
 ;
 ;
 ;
 ;  This file is part of HEALPix.
 ;
 ;  HEALPix is free software; you can redistribute it and/or modify
 ;  it under the terms of the GNU General Public License as published by
 ;  the Free Software Foundation; either version 2 of the License, or
 ;  (at your option) any later version.
 ;
 ;  HEALPix is distributed in the hope that it will be useful,
 ;  but WITHOUT ANY WARRANTY; without even the implied warranty of
 ;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ;  GNU General Public License for more details.
 ;
 ;  You should have received a copy of the GNU General Public License
 ;  along with HEALPix; if not, write to the Free Software
 ;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 ;
 ;  For more information about HEALPix see http://healpix.jpl.nasa.gov
 ;
 ; -----------------------------------------------------------------------------
 pro LS_proj2out, planmap, Tmax, Tmin, color_bar, dx, title_display, sunits,$
               coord_out, do_rot, eul_mat, planvec, vector_scale, $
               CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, $
               HXSIZE = hxsize, NOBAR = nobar, NOLABELS = nolabels, NOPOSITION = noposition, $
               PREVIEW = preview, $
               PS = ps, PXSIZE=pxsize, PYSIZE=pysize, ROT=rot_ang, SUBTITLE = subtitle, $
               TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
               POLARIZATION=polarization, $
               PNG = png, OUTLINE = outline, $
               PROJECTION=projection, MOLLWEIDE=mollweide, GNOMIC=gnomic, CARTESIAN=cartesian, $
               ORTHOGRAPHIC=orthographic, FLIP=flip, HALF_SKY=half_sky,COORD_IN=coord_in, $
               IGRATICULE = igraticule, HBOUND = hbound, DIAMONDS = diamonds, WINDOW = window_user, $
               TRANSPARENT = transparent, EXECUTE=execute, SILENT=silent, GLSIZE=glsize, IGLSIZE=iglsize, $
               SHADEMAP=SHADEMAP, RETAIN=retain, TRUECOLORS=truecolors, CHARTHICK=charthick, $
               STAGGER=stagger, AZEQ=azeq, JPEG=jpeg, keep_file_open=keep_file_open, bwhite=bwhite,$
               CTDIR=CTDIR, CTFILE=CTFILE, GRMIN=GRMIN, GRMAX=GRMAX, GRLS=GRLS, IGRMIN=IGRMIN, IGRMAX=IGRMAX, IGRLS=IGRLS, CBLBL=CBLBL

 ;

 ;===============================================================================
 ;+
 ;  PROJ2OUT
 ;  ouputs on X-screen or PS, GIF, PNG or JPEG file a gnomonic,
 ;  mollweide, cartesian, orthographic or azimuth equidistant map
 ;
 ;  IN:
 ;    planmap, Tmax, Tmin, color_bar, dx, title_display, sunits, 
 ;    coord_out, do_rot, eul_mat, planvec, vector_scale
 ;
 ;  KEYWORDS:
 ;     CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, HXSIZE = hxsize, $
 ;              NOBAR = nobar, NOLABELS = nolabels, NOPOSITION = noposition, $
 ;              PREVIEW = preview, PS = ps, $
 ;              PXSIZE=pxsize, PYSIZE=pysize, ROT = rot_ang, SUBTITLE = subtitle, $
 ;              TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
 ;              POLARIZATION=polarization,$
 ;              PNG=png, OUTLINE = outline,$
 ;              PROJECTION=projection, MOLLWEIDE=mollweide, $
 ;              GNOMIC=gnomic, CARTESIAN=cartesian,
 ;              ORTHOGRAPHIC=orthographic, $
 ;              FLIP=flip, HALF_SKY=half_sky,COORD_IN=coord_in, IGRATICULE=,
 ;              HBOUND=, DIAMONDS =, WINDOW =, TRANSPARENT=, EXECUTE=, SILENT=
 ;              GLSIZE=, IGLSIZE=, SHADEMAP=, STAGGER=, AZEQ=, JPEG=
 ;
 ;   for more information, see Gnomview.pro Mollview.pro
 ;
 ;   March 1999, EH
 ;   Nov 2000, EH, plot polarisation
 ;   May 2002, EH, merge gnom2out and moll2out
 ;   Jun 2002, EH, added the cartesian projection facility (hacked from
 ;       G. Giardino's pol2out)
 ;   Aug 2002, EH, added the orthographic projection facility
 ;   Jul 2002, EH, changed vector field loop index to LONG
 ;   Jan 2007, EH, added window keyword
 ;   Sep 2007, EH, the /CROP-ped image now include graticules, ...,
 ;   added /TRANSPARENT, EXECUTE=, /SILENT
 ;   May 2009, EH, added SHADEMAP (shade for orthographic PNG output)
 ;                 a single call to tvrd()
 ;                 uses Z buffer when window<0
 ;                 introduce RETAIN
 ;   Sep 2009, EH, TRANSPARENT can now be in {0,1,2,3}
 ;   Nov 2009, EH, retrofitting for GDL. 
 ;                 Everything works as in IDL except PS outputs and
 ;                 transparent pixels in PNG
 ;   Mar 2010, EH, corrected bug with use_z_buffer
 ;   Apr 2010, EH, accepts array of OUTLINE;
 ;                  supports CHARTHICK
 ;   Jan 2012, EH, turns off GRAT, IGRAT, HBOUND, OUTLINE when STAGGER is set
 ;                 added support of AZEQ and JPEG
 ;
 ;
 ; 2 problems with write_png,...,/transparent in GDL:
 ;  - it is currently (v0.9.2) not supported
 ;  - it expects a pixel mask (ie an array of the same size as the image) with
 ;    values in [0,255] setting the opacity of each pixel
 ;   while the IDL version expects a color mask (ie, a vector of size 255) with
 ;   values in [0,255] setting the opacity of each color
 ;   
 ;   Jan 2013, L. Spencer: Added CTDIR and CTFILE to allow a custom colortable file to be used (i.e. the Planck `parchment1' colourtable).
 ;                         Added GRMIN and GRMAX keywords [and IGRMIN and IGRMAX] to limit labelling the graticules on the edges of a mollview map.
 ;                         Added GRLS and IGRLS keyowrds to over-ride the linestyle of the graticules and igraticules.  
 ;-
 ;===============================================================================
 ;
 ;WIN = STRCOMP(!VERSION.OS, 'Win',3)
 ;IF WIN EQ 1 THEN ds = '\' ELSE ds = '/'
 IF N_ELEMENTS(CTFILE) EQ 0 THEN CTFILE = 'colors1.tbl'
 IF N_ELEMENTS(CTDIR) GT 0 THEN ColTabFl = CTDIR+CTFILE ELSE ColTabFl = !DIR+'//resource//colors//'+CTFILE
 ;
 IF N_ELEMENTS(GRLS) EQ 0 THEN GRLS = 0
 IF N_ELEMENTS(IGRLS) EQ 0 THEN BEGIN
   IF KEYWORD_SET(GRATICULE) THEN IGRLS = 1 ELSE IGRLS = 0
 ENDIF
 ;
 identify_projection, projtype, projection=projection, mollweide=mollweide, gnomic=gnomic, cartesian=cartesian, orthographic=orthographic,  diamonds = diamonds , azeq=azeq
 do_gnom = 0
 do_moll = 0
 do_cart = 0
 do_orth = 0
 do_azeq = 0
 do_fullsky = 0 ; dummy, only matters for orthview
 do_gif = keyword_set(gif)
 do_png = keyword_set(png)
 do_ps  = keyword_set(ps)
 do_jpeg  = keyword_set(jpeg)
 do_image = (do_gif || do_png || do_jpeg)
 do_true = keyword_set(truecolors)
 do_crop = keyword_set(crop)

 if undefined(polarization) then polarization=0
 do_polamplitude = (polarization[0] eq 1)
 do_poldirection = (polarization[0] eq 2)
 do_polvector    = (polarization[0] eq 3)
 ;-------------------------------------------------
 in_gdl = is_gdl()
 in_idl = ~in_gdl
 ;if (do_ps) then 
 test_preview
 @idl_default_previewer ; defines the paper size
 if (do_ps and undefined(papersize)) then papersize = 'a4'

 xsize = (size(planmap))(1)
 ysize = (size(planmap))(2)

 If (projtype eq 2) then begin
 ;  ---- Gnomonic specific definitions for the plot ----
     routine    = 'GNOMVIEW'
     proj_small = 'gnomic'
     proj_big   = 'Gnomic'
     do_gnom    = 1

 ;     du_dv = 1.                  ; aspect ratio
     du_dv = xsize/float(ysize)                  ; aspect ratio
     fudge = 1.00                ; 
     xc = (xsize-1)/2. & delta_x = (xsize-1 - xc)
     yc = (ysize-1)/2. & delta_y = (ysize-1 - yc)
 ; u and v range of the map
     umin = - dx * xc * fudge & umax = dx * xc * fudge
     vmin = - dx * yc * fudge & vmax = dx * yc * fudge
 ; position of the rectangle in the final window
     w_xll = 0.00 & w_xur = 1.00 & w_dx = w_xur - w_xll
     w_yll = 0.10 & w_yur = 0.90 & w_dy = w_yur - w_yll
     w_dx_dy = w_dx / w_dy       ; 1.4
 ; color bar, position, dimension
     cbar_dx = 1./3.
     cbar_dy = 1./70.
     cbar_xll = (1. - cbar_dx)/2.
     cbar_xur = (1. + cbar_dx)/2.
     cbar_yur = w_yll - cbar_dy
     cbar_yll = cbar_yur - cbar_dy
 ; polarisation color ring, position, dimension
     cring_dx = 1./15.
     cring_dy = 1./15.
     cring_xll = .025
 ; cring_xur = .025 + cring_dx
 ; cring_yur = w_yll
 ; cring_yll = cring_yur - cring_dy
     cring_yll = .025
 ; location of astro. coordinate
     x_aspos = 0.5
     y_aspos = 0.04
 ; location of pol vector scale
     vscal_x = 0.05
     vscal_y = 0.01
 ; location of title and subtitle
     x_title = 0.5 & y_title = 0.95
     x_subtl = 0.5 & y_subtl = 0.915
     if (do_ps) then begin
 ; default X dimension of hardcopy (cm)
         hxsize_def = 15.
 ; offset along the long axis of the page
         yoffset = (papersize eq 'a4') ? 2 : 1
 ;yoffset = 2  ; Europe (A4)
 ;yoffset = 1                 ; US (letter)
     endif
 endif

 if (projtype eq 1) then begin
 ;  ---- Mollweide specific definitions for the plot ----
     routine    = 'MOLLVIEW'
     proj_big   = 'Mollweide'
     proj_small = 'mollweide'
     do_moll    = 1

     du_dv = 2.                  ; aspect ratio
     fudge = 1.02                ; spare some space around the Mollweide egg
     xc = 0.5*(xsize-1) & delta_x = (xsize-1 - xc)
     yc = 0.5*(ysize-1) & delta_y = (ysize-1 - yc)
 ; x and y range of egg
     umin = - du_dv * fudge & umax = du_dv * fudge
     vmin = - fudge         & vmax =         fudge
 ; position of the egg in the final window
     w_xll = 0.0 & w_xur = 1.0 & w_dx = w_xur - w_xll
     ;w_yll = 0.1 & w_yur = 0.9 & w_dy = w_yur - w_yll
     w_yll = 0.16 & w_yur = 0.96 & w_dy = w_yur - w_yll
     w_dx_dy = w_dx / w_dy       ; 1./.8
 ; color bar, position, dimension
     ;cbar_dx = 1./3.
     ;cbar_dy = 1./70.
     cbar_dx = 2./3.
     cbar_dy = 1./32.
     cbar_xll = (1. - cbar_dx)/2.
     cbar_xur = (1. + cbar_dx)/2.
     cbar_yur = w_yll - cbar_dy
     cbar_yll = cbar_yur - cbar_dy
 ; polarisation color ring, position, dimension
     cring_dx = 1./10.
     cring_dy = 1./10.
     cring_xll = .025
 ; cring_xur = .025 + cring_dx
 ; cring_yur = w_yll
 ; cring_yll = cring_yur - cring_dy
     cring_yll = .025
 ; location of pol vector scale
     vscal_x = 0.05
     vscal_y = 0.02
 ; location of title and subtitle
     x_title = 0.5 & y_title = 0.95
     x_subtl = 0.5 & y_subtl = 0.905
     if (do_ps) then begin
 ; default X dimension of hardcopy (cm)
         hxsize_def = 26.
 ; offset along the long axis of the page
         yoffset = (papersize eq 'a4') ? 2 : 1
     ;yoffset = 2  ; Europe (A4)
     ;yoffset = 1                 ; US (letter)
     endif
 endif

 if (projtype eq 5) then begin
 ;  ---- Diamonds specific definitions for the plot ----
     routine    = 'DIAMONDS'
     proj_big   = 'Diamonds'
     proj_small = 'diamonds'
     do_moll    = 1

     du_dv = 2.                  ; aspect ratio
     fudge = 1.00                ; spare some space around the 12-Diamonds
     xc = 0.5*(xsize-1) & delta_x = (xsize-1 - xc)
     yc = 0.5*(ysize-1) & delta_y = (ysize-1 - yc)
 ; x and y range of egg
     umin = - du_dv * fudge & umax = du_dv * fudge
     vmin = - fudge         & vmax =         fudge
 ; position of the egg in the final window
     w_xll = 0.0 & w_xur = 1.0 & w_dx = w_xur - w_xll
     w_yll = 0.1 & w_yur = 0.9 & w_dy = w_yur - w_yll
     w_dx_dy = w_dx / w_dy       ; 1./.8
 ; color bar, position, dimension
     cbar_dx = 1./3.
     cbar_dy = 1./70.
     cbar_xll = (1. - cbar_dx)/2.
     cbar_xur = (1. + cbar_dx)/2.
     cbar_yur = w_yll - cbar_dy
     cbar_yll = cbar_yur - cbar_dy
 ; polarisation color ring, position, dimension
     cring_dx = 1./10.
     cring_dy = 1./10.
     cring_xll = .025
 ; cring_xur = .025 + cring_dx
 ; cring_yur = w_yll
 ; cring_yll = cring_yur - cring_dy
     cring_yll = .025
 ; location of pol vector scale
     vscal_x = 0.05
     vscal_y = 0.02
 ; location of title and subtitle
     x_title = 0.5 & y_title = 0.95
     x_subtl = 0.5 & y_subtl = 0.905
     if (do_ps) then begin
 ; default X dimension of hardcopy (cm)
         hxsize_def = 26.
 ; offset along the long axis of the page
         yoffset = (papersize eq 'a4') ? 2 : 1
     ;yoffset = 2  ; Europe (A4)
     ;yoffset = 1                 ; US (letter)
     endif
 endif

 if (projtype eq 4) then begin
 ;  ---- Orthographic specific definitions for the plot ----
     routine    = 'ORTHVIEW'
     proj_big   = 'Orthographic'
     proj_small = 'orthographic'
     do_orth    = 1

     if keyword_set(half_sky) then do_fullsky = 0 else do_fullsky = 1
     if (do_fullsky) then du_dv = 2. else du_dv = 1. ; aspect ratio

     fudge = 1.02                ; spare some space around the Orthographic disc
     xc = 0.5*(xsize-1) & delta_x = (xsize-1 - xc)
     yc = 0.5*(ysize-1) & delta_y = (ysize-1 - yc)
 ; x and y range of disc
     umin = - du_dv * fudge & umax = du_dv * fudge
     vmin = - fudge         & vmax =         fudge
 ; position of the disc in the final window
     w_xll = 0.0 & w_xur = 1.0 & w_dx = w_xur - w_xll
     w_yll = 0.1 & w_yur = 0.9 & w_dy = w_yur - w_yll
     w_dx_dy = w_dx / w_dy       ; 1./.8
 ; color bar, position, dimension
     cbar_dx = 1./3.
     cbar_dy = 1./70.
     cbar_xll = (1. - cbar_dx)/2.
     cbar_xur = (1. + cbar_dx)/2.
     cbar_yur = w_yll - cbar_dy
     cbar_yll = cbar_yur - cbar_dy
 ; polarisation color ring, position, dimension
     cring_dx = 1./10.
     cring_dy = 1./10.
     cring_xll = .025
 ; cring_xur = .025 + cring_dx
 ; cring_yur = w_yll
 ; cring_yll = cring_yur - cring_dy
     cring_yll = .025
 ; location of pol vector scale
     vscal_x = 0.05
     vscal_y = 0.02
 ; location of title and subtitle
     x_title = 0.5 & y_title = 0.95
     x_subtl = 0.5 & y_subtl = 0.905
     if (do_ps) then begin
 ; default X dimension of hardcopy (cm)
         hxsize_def = 26.
 ; offset along the long axis of the page
         yoffset = (papersize eq 'a4') ? 2 : 1
     ;yoffset = 2  ; Europe (A4)
     ;yoffset = 1                 ; US (letter)
     endif
 endif

 if (projtype eq 3) then begin
     ;------------ cartesion projection ----------------
     routine    = 'CARTVIEW'
     proj_small = 'cartesian'
     proj_big   = 'Cartesian'
     do_cart    = 1

 ;     du_dv = 1.                  ; aspect ratio
     du_dv = xsize/float(ysize)                  ; aspect ratio
     fudge = 1.00                ; 
     xc = (xsize-1)/2. & delta_x = (xsize-1 - xc)
     yc = (ysize-1)/2. & delta_y = (ysize-1 - yc)
 ; u and v range of the map
     umin = - dx * xc * fudge & umax = dx * xc * fudge
     vmin = - dx * yc * fudge & vmax = dx * yc * fudge
 ; position of the rectangle in the final window
     w_xll = 0.00 & w_xur = 1.00 & w_dx = w_xur - w_xll
     w_yll = 0.10 & w_yur = 0.90 & w_dy = w_yur - w_yll
     w_dx_dy = w_dx / w_dy       ; 1.4
 ; color bar, position, dimension
     cbar_dx = 1./3.
     cbar_dy = 1./70.
     cbar_xll = (1. - cbar_dx)/2.
     cbar_xur = (1. + cbar_dx)/2.
     cbar_yur = w_yll - cbar_dy
     cbar_yll = cbar_yur - cbar_dy
 ; polarisation color ring, position, dimension
     cring_dx = 1./15.
     cring_dy = 1./15.
     cring_xll = .025
     cring_yll = .025
 ; location of astro. coordinate
     x_aspos = 0.5
     y_aspos = 0.04
 ; pol vector scale
     vscal_x = 0.05
     vscal_y = 0.01
 ; location of title and subtitle
     x_title = 0.5 & y_title = 0.95
     x_subtl = 0.5 & y_subtl = 0.915
     if (do_ps) then begin
 ; default X dimension of hardcopy (cm)
         hxsize_def = 15.
 ; offset along the long axis of the page
         yoffset = (papersize eq 'a4') ? 2 : 1
     ;yoffset = 2  ; Europe (A4)
     ;yoffset = 1                 ; US (letter)
     endif
 endif

 if (projtype eq 6) then begin
     ;------------ azimuthal equidistant projection ----------------
     routine    = 'AZEQVIEW'
     proj_small = 'azimequid'
     proj_big   = 'AzimEquidistant'
     do_cart    = 1

 ;     du_dv = 1.                  ; aspect ratio
     du_dv = xsize/float(ysize)                  ; aspect ratio
     fudge = 1.00                ; 
     xc = (xsize-1)/2. & delta_x = (xsize-1 - xc)
     yc = (ysize-1)/2. & delta_y = (ysize-1 - yc)
 ; u and v range of the map
     umin = - dx * xc * fudge & umax = dx * xc * fudge
     vmin = - dx * yc * fudge & vmax = dx * yc * fudge
 ; position of the rectangle in the final window
     w_xll = 0.00 & w_xur = 1.00 & w_dx = w_xur - w_xll
     w_yll = 0.10 & w_yur = 0.90 & w_dy = w_yur - w_yll
     w_dx_dy = w_dx / w_dy       ; 1.4
 ; color bar, position, dimension
     cbar_dx = 1./3.
     cbar_dy = 1./70.
     cbar_xll = (1. - cbar_dx)/2.
     cbar_xur = (1. + cbar_dx)/2.
     cbar_yur = w_yll - cbar_dy
     cbar_yll = cbar_yur - cbar_dy
 ; polarisation color ring, position, dimension
     cring_dx = 1./15.
     cring_dy = 1./15.
     cring_xll = .025
     cring_yll = .025
 ; location of astro. coordinate
     x_aspos = 0.5
     y_aspos = 0.04
 ; pol vector scale
     vscal_x = 0.05
     vscal_y = 0.01
 ; location of title and subtitle
     x_title = 0.5 & y_title = 0.95
     x_subtl = 0.5 & y_subtl = 0.915
     if (do_ps) then begin
 ; default X dimension of hardcopy (cm)
         hxsize_def = 15.
 ; offset along the long axis of the page
         yoffset = (papersize eq 'a4') ? 2 : 1
     ;yoffset = 2  ; Europe (A4)
     ;yoffset = 1                 ; US (letter)
     endif
 endif
 ;====================================================

 do_shade = (do_orth && defined(shademap))
 ; set color table and character size
 ct          = defined(colt)     ? colt     : 33
 charsfactor = defined(charsize) ? charsize : 1.0
 mycharthick = defined(charthick)? charthick : 1.0
 be_verbose  = ~keyword_set(silent)

 ; alter the color table
 ; -----------------------
 if (be_verbose) then print,'... computing the color table ...'
 if (do_true) then begin
     loadct, 0, /silent, FILE=ColTabFl
     tvlct,red,green,blue,/get
 endif else begin
     if (do_poldirection) then begin
         LOADCT, 0, /SILENT, FILE=ColTabFl
         ncol = 256
         one = replicate(1.,ncol)
         tvlct,[0,0,0,findgen(ncol-3)]/(ncol-3)*720,one,one,/hsv ; hue is in degrees
     endif else begin
         loadct, get_names = color_names, FILE=ColTabFl
         nmax_col = n_elements(color_names)-1
         if (abs(ct) le nmax_col) then begin
             LOADCT, abs(ct), /SILENT, FILE=ColTabFl
         endif else begin
             if (be_verbose) then print,'... color table '+strtrim(abs(ct),2)+' unknown, using current color table instead ...'
         endelse
     endelse
     tvlct,red,green,blue,/get
     if (ct lt 0) then begin
         red = reverse(red) & green = reverse(green) & blue = reverse(blue)
     endif
 endelse
 ; set up some specific definitions
 ; reserve first colors for Black, White and Neutral grey
 idx_black = 0B & idx_white = 1B   & idx_grey = 2B   & idx_bwg = [idx_black, idx_white, idx_grey]
 col_black = 0B & col_white = 255B & col_grey = 175B & col_bwg = [col_black, col_white, col_grey]
;;if keyword_set(bwhite) then print, '******** Using White back ground *******'
;;if not keyword_set(bwhite) then print, '******** Using Grey back ground *******'
if keyword_set(bwhite) then begin
   col_black = 0B & col_white = 255B & col_grey = 175B & col_bwg = [col_black, col_white, col_white]
endif
 red  [idx_bwg] = col_bwg
 green[idx_bwg] = col_bwg
 blue [idx_bwg] = col_bwg
 TVLCT,red,green,blue

 ; ---------------------
 ; open the device
 ; ---------------------
 old_device=!d.name
 my_background = !p.background
 my_color = !p.color
 if (be_verbose) then print,'... here it is.'
 titlewindow = proj_big+' projection : ' + title_display
 back      = REPLICATE(BYTE(!P.BACKGROUND),xsize,(ysize*cbar_dy*w_dx_dy)>1)
 use_z_buffer = 0 ; set it to 0 (for postscript) 2010-03-18
 if (do_ps) then begin
     ; 2009-11-04: 'ps' in GDL does not support: COLOR, BITS, XSIZE, ...
     if DEFINED(hxsize) then hxsize = (hxsize > 3) < 200 else hxsize = hxsize_def
     if ((size(ps))(1) ne 7) then file_ps = 'plot_'+proj_small+'.ps' else file_ps = ps
     SET_plot,'ps'
     do_portrait = 0
     do_landscape = 0
     DEVICE, FILE=file_ps, /COLOR, BITS = 8 ; opens the file that will receive the PostScript plot
     if (do_gnom || (do_orth && ~do_fullsky)) then begin
         do_portrait = 1
         DEVICE, /PORTRAIT,  XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=4, YOFFSET=2
     endif
     if (do_moll || (do_orth && do_fullsky)) then begin
         ;do_landscape = 1
         ;DEVICE, /LANDSCAPE, XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=4, YOFFSET=hxsize+yoffset
         do_portrait = 1
         DEVICE, /PORTRAIT, XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=4, YOFFSET=hxsize+yoffset, /ENCAPSULATED, /HELVETICA, FONT_size=8
     endif
     if (do_cart || do_azeq) then begin
         do_landscape = 1
 ;         DEVICE, /LANDSCAPE, XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=4, YOFFSET=hxsize+yoffset
         DEVICE, /LANDSCAPE, XSIZE=hxsize, YSIZE=hxsize/du_dv*w_dx_dy, XOFFSET=0, YOFFSET=hxsize+yoffset
     endif
     TVLCT,red,green,blue
     thick_dev = 1. ; device dependent thickness factor
 endif else begin ; X, png, gif or jpeg output
     idl_window = defined(window_user) ? window_user : 32 ; idl_window = 32 or window_user
     free_window    =  (idl_window gt 31) ; random  window if idl_window > 31
     virtual_window =  (idl_window lt 0)  ; virtual window if idl_window < 0
     reuse_window   =  (~free_window && ~virtual_window && !d.window eq idl_window && !d.x_size eq long(xsize) && !d.y_size eq long(ysize*w_dx_dy))
     use_z_buffer   = (virtual_window && do_image)
     window_retain  = defined(retain) ? retain : 2
     if (use_z_buffer) then begin
         character_size = [!d.x_ch_size,!d.y_ch_size]
         set_plot,'z'
         pixel_depth = (do_true) ? 24 : 8
         if (in_gdl) then begin
 ; unresolved GDL0.9.2 bug: set_character_size ignored
             device,set_resolution= [xsize, ysize*w_dx_dy], set_character_size=character_size,z_buffering=0
         endif else begin
             device,set_resolution= [xsize, ysize*w_dx_dy], set_character_size=character_size,z_buffering=0, set_pixel_depth=pixel_depth
         endelse
     endif
     ;;if (!D.NAME eq 'X') then  DEVICE, PSEUDO = 8 ; for Windows compatibility ;
     ;;commented out 2009-10-28
     ;to_patch = ((!d.n_colors GT 256) && do_image  && ~do_crop)
     ;to_patch = ((!d.n_colors GT 256) && do_image && in_idl)
     n_colors = !d.n_colors
     if (in_gdl && (!d.name eq 'X' || !d.name eq 'WIN')) then begin ; work-around for GDL0.9.2 bug (!d.n_colors gets correct only after first call to WINDOW)
         device,get_visual_depth=gvd
         n_colors = 2L^gvd
     endif
     to_patch = (n_colors GT 256 && do_image)
     if (in_gdl) then begin
         if (use_z_buffer) then device else device,decomposed=0 ; in GDL0.9.2, decomposed is only available (but ignored) when device='X', or unvalid when device='Z'
         if (to_patch) then loadct,0,/silent, FILE=ColTabFl ; emulate decomposed
     endif else begin
         device, decomposed = use_z_buffer || to_patch
     endelse
     if (reuse_window) then begin
         wset, idl_window
     endif else begin
         if (~use_z_buffer) then begin
             WINDOW, idl_window>0, FREE=free_window, PIXMAP=virtual_window, $
                     XSIZE = xsize, YSIZE = ysize*w_dx_dy, TITLE = titlewindow, $
                     XPOS=(in_gdl && undefined(xpos)) ? 0 : xpos, $
                     YPOS=(in_gdl && undefined(ypos)) ? 0 : ypos, $ ; bug correction 2009-12-05
                     RETAIN=window_retain
         endif
         if (~virtual_window && (!d.x_size lt long(xsize) || !d.y_size lt long(ysize*w_dx_dy))) then begin
             message_patch,level=-1,/info,'==========================================================='
             message_patch,level=-1,/info,'WARNING: Because of screen and window manager limitations,'
             message_patch,level=-1,/info,'         the actual window is not as large as expected !'
             message_patch,level=-1,/info,strtrim(!d.x_size,2)+'*'+  strtrim(!d.y_size,2)+'    <    '+  strtrim(long(xsize),2)+'*'+strtrim(long(ysize*w_dx_dy),2)
             message_patch,level=-1,/info,'         The result is unpredictable.'            
             message_patch,level=-1,/info,' If you are only interested in GIF/PNG/JPEG output, you can use a virtual window (WINDOW<0) instead'            
             message_patch,level=-1,/info,'==========================================================='
         endif
     endelse
     if (in_idl) then TVLCT,red,green,blue
     thick_dev = 1. ; device dependent thickness factor
 endelse
 !p.background = my_background
 !p.color = my_color
 ; -------------------------------------------------------------
 ; make the plot
 ; -------------------------------------------------------------
 myplot={urange:[umin,umax],vrange:[vmin,vmax],position:[w_xll,w_yll,w_xur,w_yur],xstyle:5,ystyle:5}
 plot, /nodata, myplot.urange, myplot.vrange, pos=myplot.position, XSTYLE=myplot.xstyle, YSTYLE=myplot.ystyle
 ; Set the background when doing a plot in Z buffer
 if (use_z_buffer) then begin
     l64ysize = long64(ysize*w_dx_dy)
     if ~do_true then begin 
         tv, replicate(!p.background, xsize, l64ysize)
     endif else begin
         back = [red[!p.background],green[!p.background],blue[!p.background]] ## replicate(1b, xsize*l64ysize)
         tv, reform(back, xsize, l64ysize, 3, /overwrite),true=3
         back=0
     endelse
 endif
 ; ---------- projection independent ------------------
 ; map itself
 if (do_shade && ~do_image) then begin
     ; shaded for X or PS
     image = planmap
     image3d  =   make_array(/uint, xsize, ysize, 3)
     image3d[*,*,0] = uint( (256. * red  [image] * shademap) < 65535.)
     image3d[*,*,1] = uint( (256. * green[image] * shademap) < 65535.)
     image3d[*,*,2] = uint( (256. * blue [image] * shademap) < 65535.)
     if (do_ps) then loadct,0,/silent, FILE=ColTabFl ; must be in grey-scale for TrueColor PS output
     TV, bytscl(image3d),w_xll,w_yll,/normal,xsize=1.,true=3
     if (do_ps) then tvlct,red,green,blue ; revert to custom color table
     image3d = 0
 ; endif else if (do_true  && ~do_image) then begin
 endif else if (do_true) then begin
                                 ; truecolors for X or PS, red, green, blue are
                                 ; only useful for the {0,1,2}={Black, white, grey}
     image3d = make_array(/byte, xsize, ysize, 3)
     image3d[*,*,0] = red  [planmap[*,*,0]]
     image3d[*,*,1] = green[planmap[*,*,1]]
     image3d[*,*,2] = blue [planmap[*,*,2]]
 ;;;    if (do_ps) then loadct,0,/silent; must be in grey-scale for TrueColor PS output
     tv, image3d, w_xll,w_yll,/normal,xsize=1.,true=3
 ;;;    if (do_ps) then tvlct,red,green,blue ; revert to custom color table
     image3d = 0
 endif else begin
     TV, planmap,w_xll,w_yll,/normal,xsize=1.
 endelse
 hpxv11 = 0

 ; the polarisation color ring
 if (do_poldirection) then begin
     npring = xsize*cring_dx
     one = replicate(1.,npring)
     yy  = one # (findgen(npring) - npring/2)
     xx  = transpose(yy)
     tt  = (xx^2 + yy^2) / float(npring)^2
     mask = byte(tt lt .25 and tt gt 0.05)
     if (hpxv11) then begin
         ; angles are measured from horizontal
         angle = atan(yy,xx)  * !radeg + 180. ; angle in deg in [0,360]
     endif else begin
         ; angles are measured from vertical
         angle = atan(xx,-yy)  * !radeg + 180. ; angle in deg in [0,360]
     endelse
     color_ring = (bytscl(angle,TOP=252) + 3B) * mask + byte((1-mask)*!P.BACKGROUND); in [3,255] in ring and !p.background outside ring
     TV,color_ring,cring_xll,cring_yll,/normal,xsize = npring/float(xsize)
 endif

 ; polarisation vector field
 pgparam=[1.,1.]; [scale, grid_step] of grid of headless vector
 if (do_polvector) then begin
     pgparam = ([polarization[*], 1., 1.])[1:2] ; 2nd and 3rd elements of polarization (default:[1,1])
     pgparam = pgparam*(pgparam gt 0) + (pgparam le 0) ; replace non-positive values by 1.
     dyg = 10.
     pol_rescale = float(dyg)/ysize * pgparam[0]
     dyg *= pgparam[1] & dxg = dyg
     xg = (lindgen(xsize/dxg)+.5) #  replicate(dxg, ysize/dyg)
     yg = (lindgen(ysize/dyg)+.5) ## replicate(dyg, xsize/dxg)
     u = umin + xg * (umax-umin) / xsize
     v = vmin + yg * (vmax-vmin) / ysize
     for i=0L, n_elements(xg)-1 do begin
         norm = planvec[xg[i],yg[i],0] * pol_rescale * (vmax-vmin)
         angle = planvec[xg[i],yg[i],1]
         if (hpxv11) then begin
             ; angles are measured from horizontal
             if (norm gt 0) then plots, u[i]+norm*cos(angle)*[-.5,.5], v[i]+norm*sin(angle)*[-.5,.5]
         endif else begin
             ; angles are measured from vertical
             if (norm gt 0) then plots, u[i]-norm*sin(angle)*[-.5,.5], v[i]+norm*cos(angle)*[-.5,.5]
         endelse
     endfor
     xg = 0 & yg = 0 & u = 0 & v = 0
 endif

 ;  the color bar
 if (~(keyword_set(nobar) || do_poldirection || do_true)) then begin
     color_bar_out = BYTE(CONGRID(color_bar,xsize*cbar_dx)) # REPLICATE(1.,(ysize*cbar_dy*w_dx_dy)>1)
     back(xsize*cbar_xll,0) = color_bar_out
     TV, back,0,cbar_yll,/normal,xsize = 1.
 endif

 ;  the color bar labels
 if (~(keyword_set(nobar) || keyword_set(nolabels) || do_true || do_poldirection)) then begin
     format = '(g10.2)'
     if ((Tmax - Tmin) ge 50 and MAX(ABS([Tmax,Tmin])) le 1.e5) then format='(i8)'
     if ((Tmax - Tmin) ge 5  and MAX(ABS([Tmax,Tmin])) le 1.e2) then format='(f6.1)'
     strmin = STRING(Tmin,format=format)
     strmax = STRING(Tmax,format=format)
     ;XYOUTS, cbar_xll, cbar_yll,'!6'+STRTRIM(strmin,2)+' ',$
     ;        ALIGN=1.,/normal, chars=1.3*charsfactor, charthick=mycharthick
     ;XYOUTS, cbar_xur, cbar_yll,'!6 '+STRTRIM(strmax,2)+' '+sunits,$
     ;        ALIGN=0.,/normal, chars=1.3*charsfactor, charthick=mycharthick
     ;  LSedit to get colorbar label to work with latex/ psfrag
     IF KEYWORD_SET(CBLBL) THEN BEGIN
       XYOUTS, cbar_xll, cbar_yll,'!X'+STRTRIM(strmin,2)+' ',$
             ALIGN=1.,/normal, chars=1.3*charsfactor, charthick=mycharthick
       XYOUTS, cbar_xur, cbar_yll,' '+STRTRIM(strmax,2)+' '+sunits,$
             ALIGN=0.,/normal, chars=1.3*charsfactor, charthick=mycharthick
       XYOUTS, Cbar_xll + cbar_dx/2d, cbar_yll - cbar_dy/2d - DOUBLE(!D.Y_CH_SIZE)/DOUBLE(!D.Y_SIZE), CBLBL, $
             ALIGN=0.5,/NORMAL, chars=1.3*charsfactor, charthick=mycharthick 
     ENDIF ELSE BEGIN
       XYOUTS, cbar_xll, cbar_yll,'!X'+STRTRIM(strmin,2)+' ',$
             ALIGN=1.,/normal, chars=1.3*charsfactor, charthick=mycharthick
       XYOUTS, cbar_xur, cbar_yll,'!X '+STRTRIM(strmax,2)+' '+sunits,$
             ALIGN=0.,/normal, chars=1.3*charsfactor, charthick=mycharthick
       ;XYOUTS, Cbar_xll + cbar_dx/2d, cbar_yll - cbar_dy*2d, sunits, $
       ;      ALIGN=0.5,/NORMAL, chars=1.3*charsfactor, charthick=mycharthick 
     ENDELSE    
 endif

 ; the polarisation vector scale
 if (~keyword_set(nobar)  && do_polvector) then begin
     vp_plot = 5*pol_rescale[0] /pgparam[0]; length of ruler on plot
     vp_phys = 5*vector_scale[0]/pgparam[0] ; 'physical' length of ruler
     plots, vscal_x*[1,1], vscal_y+[0, vp_plot]*w_dy, /normal
     format = '(g10.2)'
     if (vp_phys lt 1.e3 && vp_phys ge 10)    then format = '(f5.1)'
     if (vp_phys lt 10   && vp_phys gt 1.e-1) then format = '(f4.2)'
     xyouts, vscal_x, vscal_y + .5*(vp_plot)*w_dy, '!6  '+strtrim(string(vp_phys,form=format),2)+' '+sunits,ALIGN=0.,/normal, chars=1.1*charsfactor, charthick=mycharthick
 endif

 ;  the title
 if (~ keyword_set(titleplot)) then title= '!6'+title_display else title='!6'+titleplot
 XYOUTS, x_title, y_title ,title, align=0.5, /normal, chars=1.6*charsfactor, charthick=mycharthick

 ;  the subtitle
 if (keyword_set(subtitle)) then begin
     XYOUTS, x_subtl, y_subtl ,'!6 '+subtitle, align=0.5, /normal, chars=1.6*.75*charsfactor, charthick=mycharthick
 endif

 ; ---------- projection dependent ------------------

 if (do_gnom) then begin
 ;  astronomical position of central point
     if (not keyword_set(noposition)) then begin
         if (undefined(rot_ang)) then rot_ang = [0.,0.,0.] else rot_ang = ([rot_ang,0,0])(0:2)
         rot_0 = STRTRIM(STRING(rot_ang(0),form='(f6.1)'),2)
         rot_1 = STRTRIM(STRING(rot_ang(1),form='(f6.1)'),2)
         XYOUTS,x_aspos,y_aspos,'('+rot_0+', '+rot_1+') '+decode_coord(coord_out),/normal,align=0.5
     endif

 ; ; cross in the middle of plot
 ; plots,0,0,ps=1,syms=5,thick=4
 ; phi = findgen(31)/30.*2*!pi
 ; x_circle = cos(phi)
 ; y_circle = sin(phi)
 ; radius = tan(1.*!dtor/2.) ; radius = fwhm/2
 ; xyouts,0.7*umax,-0.8*vmax,'100 GHz'
 ; oplot,0.92*umax+radius*x_circle,-0.8*vmax+radius*y_circle,thick=3
 ; radius = tan(1./1.5*!dtor/2.)
 ; xyouts,0.7*umax,-0.9*vmax,'150 GHz'
 ; oplot,0.92*umax+radius*x_circle,-0.9*vmax+radius*y_circle,thick=3

 endif

 ; do not plot graticules, outlines or pixel boundaries in stagger mode (orthview)
 skip_oplots = do_orth && keyword_set(stagger) && $
   ( keyword_set(graticule) || keyword_set(igraticule) || keyword_set(hbound) || keyword_set(outline))

 if (skip_oplots) then begin
     message_patch,/info,level=-1,'*Warning*: GRAT, IGRAT, HBOUND and OUTLINE keywords are ignored in STAGGER mode'
 endif else begin
     grattwice=0
 ;  the graticule in output astrophysical coordinates
     if (KEYWORD_SET(graticule)) then begin
         grattwice =1
         glabelsize = charsfactor * (keyword_set(glsize) ? glsize : 0 )

         ;stop
         LS_oplot_graticule, graticule, eul_mat, projection=proj_small, flip = flip, thick = 1.*thick_dev, color = !p.color, half_sky=half_sky, linestyle=GRLS, charsize=glabelsize, reso_rad=dx, GRMIN=GRMIN, GRMAX=GRMAX
     endif 

 ;  the graticule in input coordinates
     if (KEYWORD_SET(igraticule)) then begin
         lines_ig = 2*grattwice  ; either 0 or 2
         iglabelsize = charsfactor * (keyword_set(iglsize) ? iglsize : 0 )
         LS_oplot_graticule, igraticule, eul_mat, projection=proj_small, flip = flip, thick = 1.*thick_dev, color = !p.color, half_sky=half_sky, linestyle=IGRLS, coordsys=[coord_in,coord_out], charsize=iglabelsize, reso_rad=dx, GRMIN=IGRMIN, GRMAX=IGRMAX
     endif 

 ; outlines on the map
     if (keyword_set(outline)) then begin
         for iol=0, n_elements(outline)-1 do begin
             outline_coord2uv, outline[iol], coord_out, eul_mat, projection=proj_small, flip = flip, /show, thick = 1.*thick_dev, half_sky=half_sky
         endfor
     endif

 ; overplot pixel boundaries
     if keyword_set(hbound) then begin
         nhbd = n_elements(hbound)
         if (nhbd gt 3) then message_patch,/info,level=-1,'Hbound must have 3 elements at most'
         lnst = [0,2,1]          ; solid (smallest Nside), dashes (median Nside), dots (largest Nside)
         for i=0, (nhbd<3)-1 do begin
             if (hbound[i] gt 0) then oplot_healpix_bounds, hbound[i], eul_mat, projection=proj_small, flip = flip, thick = 1.*thick_dev, color = !p.color, half_sky=half_sky, linestyle=lnst[i], coordsys=[coord_in,coord_out]
         endfor
     endif
 endelse

 ; overplot user defined commands
 if keyword_set(execute) then begin
     junk=execute(execute)
     ; reset the plotting area for cursor to work properly
     plot, /nodata, myplot.urange, myplot.vrange, pos=myplot.position, XSTYLE=myplot.xstyle, YSTYLE=myplot.ystyle,/noerase
 endif

 ; -----------------------------------------------
 ;       output the PS/GIF/PNG/JPG
 ; -----------------------------------------------

 ;  gif/png/jpg output
 if do_image then begin
     jquality = 100 ; JPEG quality in [0,100]
     valid_transparent = 0
     if (keyword_set(transparent)) then begin
         itr = nint(transparent)
         if (itr lt 0 or itr gt 3) then begin
             message,/info,'keyword TRANSPARENT must be in {0,1,2,3}'
             message,/info,'current value '+string(transparent)+' will be ignored.'
         endif else valid_transparent = 1
     endif

     if (do_gif)  then file_image = (datatype(gif)  ne 'STR') ? 'plot_'+proj_small+'.gif'  : gif
     if (do_png)  then file_image = (datatype(png)  ne 'STR') ? 'plot_'+proj_small+'.png'  : png
     if (do_jpeg) then file_image = (datatype(jpeg) ne 'STR') ? 'plot_'+proj_small+'.jpeg' : jpeg

     image = (do_true) ? tvrd(true=3) : tvrd() ; a single call to tvrd()
     if (do_shade) then begin
         image3d  =   make_array(/uint, 3,!d.x_size,!d.y_size)
         allshade =   make_array(/float,  !d.x_size,!d.y_size,value=1.0)
         allshade[w_xll*!d.x_size,w_yll*!d.y_size] = shademap
         shademap = 0
         image3d[0,*,*] = uint( (256. * red  [image] * allshade) < 65535.)
         image3d[1,*,*] = uint( (256. * green[image] * allshade) < 65535.)
         image3d[2,*,*] = uint( (256. * blue [image] * allshade) < 65535.)
 ;         if (in_gdl) then image3d = bytscl(image3d) ; GDL's write_png won't deal correctly with 16bit integers
         image3d = bytscl(image3d) ; use 8 bit integers only
         allshade = 0
     endif
     if (do_true) then begin
         dim3d = valid_transparent ? 4 : 3
         image3d  =   make_array(/byte, dim3d, !d.x_size, !d.y_size)
         for i=0,2 do image3d[i,*,*] = image[*,*,i]
         if (valid_transparent) then begin
             image3d[3,*,*] = 255B
             if (itr   and 1) then begin ; turn grey  into transparent
                 pix_tr = where( total(abs(image3d[0:2,*,*]-col_grey ),1) eq 0, n_tr)
                 if (n_tr gt 0) then image3d[3 +4*pix_tr] = 0B
             endif
             if (itr/2 and 1) then begin ; turn white into transparent
                 pix_tr = where( total(abs(image3d[0:2,*,*]-col_white),1) eq 0, n_tr)
                 if (n_tr gt 0) then image3d[3 +4*pix_tr] = 0B
             endif
         endif
     endif
     ; deal with transparent colors for not TRUECOLORS, not SHADED images
     if ~(do_true || do_shade) then begin
         transp_colors = replicate(255B, 256) ; all colors are opaque
         if (valid_transparent) then begin
                                 ; transparent = {1,3} -> grey pixels  are transparent
                                 ; transparent = {2,3} -> white pixels are transparent
             if (itr   and 1) then transp_colors[idx_grey ] = 0B ; turn grey  into transparent
             if (itr/2 and 1) then transp_colors[idx_white] = 0B ; turn white into transparent
         endif
     endif
     if do_crop then begin
         y_crop_low = round(w_yll * n_elements(image[0,*])) & y_crop_hi  = y_crop_low + ysize - 1
         cropped = image[*,y_crop_low:y_crop_hi]
         if do_gif then write_gif,file_image,cropped,red,green,blue
         if do_png then begin
             if (do_shade || do_true) then begin
                 write_png, file_image, image3d[*,*,y_crop_low:y_crop_hi]
             endif else begin
                 if (keyword_set(transparent)) then begin
                     mytransp = (in_idl) ? transp_colors  :  0 ; transp_colors[cropped]
                     write_png,file_image,cropped,red,green,blue, transparent=mytransp
                 endif else begin
                     write_png,file_image,cropped,red,green,blue
                 endelse
             endelse
         endif
         if do_jpeg then begin
             if (do_shade || do_true) then begin
                 write_jpg_custom, file_image, image3d[*,*,y_crop_low:y_crop_hi], true=1, quality=jquality
             endif else begin
                 write_jpg_custom, file_image, cropped, red, green, blue,                 quality=jquality
             endelse
         endif
     endif else begin ; uncropped
         if do_gif then write_gif,file_image, image,red,green,blue
         if do_png then begin
             if (do_shade || do_true) then begin
                 write_png, file_image, image3d
             endif else begin
                 if (keyword_set(transparent)) then begin
                     mytransp = (in_idl) ? transp_colors  : 0 ; transp_colors[image]
                     write_png,file_image, image,red,green,blue, transparent=mytransp
                 endif else begin
                     write_png,file_image, image,red,green,blue
                 endelse
             endelse
         endif
         if do_jpeg then begin
             if (do_shade || do_true) then begin
                 write_jpg_custom, file_image, image3d,             true=1, quality=jquality
             endif else begin
                 write_jpg_custom, file_image, image,red,green,blue,        quality=jquality
             endelse
         endif
     endelse
     if (to_patch && ~use_z_buffer) then begin 
         if (in_gdl) then begin
 ; unresolved GDL0.9.2 bug: if a window is already open for a given color table
 ; (selected with loadct) subsequent tvlct are ignored for that window. Only a
 ; new loadct will do the job.
             device, decomposed=0
             tvlct,red,green,blue ; revert to custom color table
         endif else begin
             device,decomposed=0     ; put back colors on X window and redo color image
         endelse
         if (do_shade || do_true) then begin
             tv, bytscl(image3d),0,0,/normal,xsize=1.,true=1
         endif else begin
             tv, image
         endelse
     endif
     image = 0
     if (be_verbose) then print,'IMAGE file is in '+file_image
     if (keyword_set(preview)) then begin
         test_preview, found_preview ;, /crash
         if (found_preview gt 0) then begin
             resolve_routine,'preview_file',/either ; ,compile_full_file=in_idl
             if do_gif then preview_file, file_image, /gif
             if do_png then preview_file, file_image, /png
             if do_jpeg then preview_file, file_image, /jpeg
         endif
     endif
 endif

 if not keyword_set(keep_file_open) then begin
 if (do_ps) then begin
     device,/close
     set_plot,old_device
     if (be_verbose) then print,'PS file is in '+file_ps
     if (keyword_set(preview)) then begin
         test_preview, found_preview ;, /crash
         if (found_preview gt 0) then begin
             resolve_routine,'preview_file',/compile_full_file,/either
             preview_file, file_ps, /ps, landscape=do_landscape
         endif
     endif
 endif else if (use_z_buffer) then begin
     device,/close ;,decomp=~to_patch
     set_plot,old_device
 endif
 endif

 return
 end

 ; -----------------------------------------------------------------------------
 ;
 ;  Copyright (C) 1997-2010  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
 ;
 ;
 ;
 ;
 ;
 ;  This file is part of HEALPix.
 ;
 ;  HEALPix is free software; you can redistribute it and/or modify
 ;  it under the terms of the GNU General Public License as published by
 ;  the Free Software Foundation; either version 2 of the License, or
 ;  (at your option) any later version.
 ;
 ;  HEALPix is distributed in the hope that it will be useful,
 ;  but WITHOUT ANY WARRANTY; without even the implied warranty of
 ;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ;  GNU General Public License for more details.
 ;
 ;  You should have received a copy of the GNU General Public License
 ;  along with HEALPix; if not, write to the Free Software
 ;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 ;
 ;  For more information about HEALPix see http://healpix.jpl.nasa.gov
 ;
 ; -----------------------------------------------------------------------------
 PRO ls_gnomview, file_in, select_in, $
               ASINH = asinh, $
               CHARSIZE = charsize, $
               CHARTHICK = charthick, $
               COLT = colt, $
               COORD = coord, $
               CROP = crop, $
               EXECUTE=execute, $
               FACTOR = factor, $
               FITS = fits, $
               FLIP = flip, $
               GIF = gif, $
               GLSIZE = glsize, $
               GRATICULE = graticule, $
               HBOUND = hbound, $
               HELP = help, $
               HIST_EQUAL = hist_equal, $
               HXSIZE = hxsize, $
               IGLSIZE = iglsize, $
               IGRATICULE = igraticule, $
               LOG = log, $
               MAP_OUT = map_out, $
               MAX = max_set, $
               MIN = min_set, $
               NESTED = nested_online, $
               NOBAR = nobar, $
               NOLABELS = nolabels, $
               NOPOSITION = noposition, $
               OFFSET = offset, $
               ONLINE = online, $
               OUTLINE = outline, $
               PNG = png, $
               POLARIZATION = polarization, $
               PREVIEW = preview, $
               PS = ps, $
               PXSIZE = pxsize, $
               PYSIZE = pysize, $
               QUADCUBE = quadcube, $
               RESO_ARCMIN = reso_arcmin, $
               RETAIN = retain, $
               ROT = rot, $
               SAVE = save, $
               SILENT = silent, $
               SUBTITLE = subtitle, $
               TITLEPLOT = titleplot, $
               TRANSPARENT = transparent, $
               TRUECOLORS = truecolors, $
               UNITS = units, $
               WINDOW = window, $
               XPOS = xpos, $
               YPOS = ypos, $
               vector_scale = vector_scale,$
               CTDIR=CTDIR, $
               keep_file_open=keep_file_open, $
               CTFILE=CTFILE, GRMIN=GRMIN, GRMAX=GRMAX, GRLS=GRLS, IGRMIN=IGRMIN, IGRMAX=IGRMAX, IGRLS=IGRLS, CBLBL=CBLBL
 ;+
 ; for extended description see mollview or the paper documentation
 ;-

 defsysv, '!healpix', exists = exists
 if (exists ne 1) then init_healpix

 @viewcom ; define common
 data_plot = 0 ; empty common array
 ; record original color table and PLOTS settings
 record_original_settings, original_settings

 loadsky                         ; cgis package routine, define rotation matrices
 projection = 'GNOMIC'
 routine = 'gnomview'

 uroutine = strupcase(routine)
 if keyword_set(help) then begin
     doc_library,'mollview'
     return
 endif

 if keyword_set(gif) then begin
     message_gif, code=routine, error=error_gif
     if (error_gif) then return
 endif

 if (n_params() lt 1 or n_params() gt 2) then begin
     PRINT, 'Wrong number of arguments in '+uroutine
     print,'Syntax : '
     print, uroutine+', File, [Select, ]'
     print,'              [ASINH=, CHARSIZE=, COLT=, COORD=, CROP=, '
     print,'              EXECUTE=, FACTOR=, FITS=, FLIP=, GIF=, GLSIZE=, GRATICULE=, '
     print,'              HBOUND=, HELP=, '
     print,'              HIST_EQUAL=, HXSIZE=, '
     print,'              IGLSIZE=, IGRATICULE=,'
     print,'              LOG=, '
     print,'              MAP_OUT=, MAX=, MIN=, '
     print,'              NESTED=, NOBAR=, NOLABELS=, NOPOSITION = '
     print,'              OFFSET=, ONLINE=, OUTLINE=,'
     print,'              PNG=,'
     print,'              POLARIZATION=, PREVIEW=, '
     print,'              PS=, PXSIZE=, PYSIZE=, QUADCUBE= ,'
     print,'              RESO_ARCMIN=, RETAIN =, ROT=, '
     print,'              SAVE=, SILENT=, SUBTITLE=, '
     print,'              TITLEPLOT=, TRANSPARENT=, TRUECOLORS= '
     print,'              UNITS=, WINDOW=, XPOS=, YPOS=]'
     print
     print,' Type '+uroutine+', /help '
     print,'   for an extended help'
     return
 endif

 IF (undefined(file_in)) then begin
     print,routine+': Undefined variable as 1st argument'
     return
 endif
 ; file_in1   = file_in
 ; if defined(select_in) then select_in1 = select_in else select_in1=1
 ; if defined(save)      then save1 = save           else save1=0
 ; if defined(online)    then online1 = online       else online1=0
 do_flip = keyword_set(flip)

 if (!D.n_colors lt 4) then begin
     print,' : Sorry ... not enough colors ('+strtrim(string(!d.n_colors),2)+') available'
     return
 endif

 polar_type = 0
 if keyword_set(polarization) then polar_type = polarization

 loaddata_healpix, $
   file_in, select_in,$
   data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
   SAVE=save,ONLINE=online,NESTED=nested_online,UNITS=units,COORD=coord,FLIP=flip, $
   ROT=rot,QUADCUBE=quadcube,LOG=log,ERROR=error, $
   POLARIZATION=polarization, FACTOR=factor, OFFSET=offset, SILENT=silent, COMPRESS=1, PIXEL_LIST=pixel_list, $
   TRUECOLORS=truecolors, DATA_TC=data_tc
 if error NE 0 then return

 polar_type = 0
 if keyword_set(polarization) then polar_type = polarization


 LS_data2gnom, $
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
  planmap, Tmax, Tmin, color_bar, dx, planvec, vector_scale, $
  PXSIZE=pxsize, PYSIZE=pysize, ROT=rot, LOG=log, HIST_EQUAL=hist_equal, $
  MAX=max_set, MIN=min_set, $
    FITS = fits, FLIP=flip, DATA_plot = data_plot, $
  POLARIZATION=polarization, SILENT=silent, PIXEL_LIST=pixel_list, ASINH=asinh, $
  TRUECOLORS=truecolors, DATA_TC=data_tc, MAP_OUT=map_out ,  RESO_ARCMIN = reso_arcmin


LS_proj2out, $
  planmap, Tmax, Tmin, color_bar, dx, title_display, $
  sunits, coord_out, do_rot, eul_mat, planvec, vector_scale, $
  CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, HXSIZE = hxsize, $
  NOBAR = nobar, NOLABELS = nolabels, NOPOSITION = noposition, PNG = png, PREVIEW = preview, PS = ps, $
  PXSIZE=pxsize, PYSIZE=pysize, ROT = rot, SUBTITLE = subtitle, $
  TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
  POLARIZATION=polarization, OUTLINE=outline, /GNOM, FLIP=flip, COORD_IN=coord_in, IGRATICULE=igraticule, $
  HBOUND = hbound, WINDOW = window, EXECUTE=execute, SILENT=silent, GLSIZE=glsize, $
  IGLSIZE=iglsize, RETAIN=retain, TRUECOLORS=truecolors, TRANSPARENT=transparent, CHARTHICK=charthick, keep_file_open=keep_file_open,$
  JPEG=jpeg, CTDIR=CTDIR, CTFILE=CTFILE, GRMIN=GRMIN, GRMAX=GRMAX, GRLS=GRLS, IGRMIN=IGRMIN, IGRMAX=IGRMAX, IGRLS=IGRLS, CBLBL=CBLBL

  ;; planmap, Tmax, Tmin, color_bar, 0., title_display, $
  ;; sunits, coord_out, do_rot, eul_mat, planvec, vector_scale, $
  ;; CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, $
  ;; HXSIZE=hxsize, NOBAR = nobar, NOLABELS = nolabels, PNG = png, PREVIEW = preview, PS=ps, PXSIZE=pxsize, $
  ;; SUBTITLE = subtitle, TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
  ;; POLARIZATION=polarization, OUTLINE=outline, /MOLL, FLIP=flip, COORD_IN=coord_in, IGRATICULE=igraticule, $
  ;; HBOUND = hbound, WINDOW = window, EXECUTE=execute, SILENT=silent, GLSIZE=glsize, $
  ;; IGLSIZE=iglsize, RETAIN=retain, TRUECOLORS=truecolors, TRANSPARENT=transparent, CHARTHICK=charthick, keep_file_open=keep_file_open,$
  ;; JPEG=jpeg, CTDIR=CTDIR, CTFILE=CTFILE, GRMIN=GRMIN, GRMAX=GRMAX, GRLS=GRLS, IGRMIN=IGRMIN, IGRMAX=IGRMAX, IGRLS=IGRLS, CBLBL=CBLBL

;; proj2out, $


w_num = !d.window
; restore original color table and PLOTS settings
record_original_settings, original_settings, /restore

RETURN
END
; -----------------------------------------------------------------------------
;
;  Copyright (C) 1997-2012  Krzysztof M. Gorski, Eric Hivon, Anthony J. Banday
;
;
;
;
;
;  This file is part of HEALPix.
;
;  HEALPix is free software; you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation; either version 2 of the License, or
;  (at your option) any later version.
;
;  HEALPix is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with HEALPix; if not, write to the Free Software
;  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
;
;  For more information about HEALPix see http://healpix.jpl.nasa.gov
;
; -----------------------------------------------------------------------------
pro ls_mollview, file_in, select_in, $
              ASINH = asinh, $
              CHARSIZE = charsize, $
              CHARTHICK = charthick, $
              COLT = colt, $
              COORD = coord, $
              CROP = crop, $
              EXECUTE = execute, $
              FACTOR = factor, $
              FITS = fits, $
              FLIP = flip, $
              GAL_CUT = gal_cut, $
              GIF = gif, $
              GLSIZE = glsize, $
              GRATICULE = graticule, $
              HELP = help, $
              HBOUND = hbound, $
              HIST_EQUAL = hist_equal, $
              HXSIZE = hxsize, $
              IGLSIZE = iglsize, $
              IGRATICULE=igraticule, $
              JPEG=jpeg, $
              LOG = log, $
              MAP_OUT = map_out, $
              MAX = max_set, $
              MIN = min_set, $
              NESTED = nested_online, $
              NOBAR = nobar, $
              NOLABELS = nolabels, $
              NO_DIPOLE = no_dipole, $
              NO_MONOPOLE = no_monopole, $
              OFFSET = offset, $
              ONLINE = online, $
              OUTLINE = outline, $
              PNG = png, $
              POLARIZATION = polarization, $
              PREVIEW = preview, $
              PS = ps, $
              PXSIZE = pxsize, $
              QUADCUBE = quadcube, $
              RETAIN = retain, $
              ROT = rot, $
              SAVE = save, $
              SILENT = silent, $
              SUBTITLE = subtitle, $
              TITLEPLOT = titleplot, $
              TRANSPARENT = transparent, $
              TRUECOLORS = truecolors, $
              UNITS = units, $
              WINDOW = window, $
              XPOS = xpos, $
              YPOS = ypos, $
              CTDIR=CTDIR, $
              bwhite=bwhite,$
              keep_file_open=keep_file_open, $
              CTFILE=CTFILE, GRMIN=GRMIN, GRMAX=GRMAX, GRLS=GRLS, IGRMIN=IGRMIN, IGRMAX=IGRMAX, IGRLS=IGRLS, CBLBL=CBLBL

;+
; NAME:
; MOLLVIEW, GNOMVIEW, CARTVIEW, ORTHVIEW
;
; PURPOSE:
;   tools to view a Mollweide/gnomic/cartesian/orthographic projection of maps binned
; in Healpix or COBE Quad-Cube pixelisation
;
; CALLING SEQUENCE:
;   xxxxVIEW, File, [Select, ] $
;                       [ASINH=, CHARSIZE=, COLT=, COORD=, CROP=, $
;                       EXECUTE=execute, $
;                       FACTOR=, FITS=, FLIP=, $
;                       GAL_CUT=, GIF=, GLSIZE=, GRATICULE=, $
;                       HALF_SKY =, HBOUND =, HELP =, HIST_EQUAL=, HXSIZE=, $
;                       IGLSIZE=, IGRATICULE=, $
;                       JPEG=, $
;                       LOG=, $
;                       MAP_OUT=, MAX=, MIN=, $ 
;                       NESTED=, NOBAR=, NOLABELS=, NOPOSITION =, $
;                       OFFSET =, ONLINE=, OUTLINE=, $
;                       PNG=, POLARIZATION=, PREVIEW=,$
;                       PS=, PXSIZE=, PYSIZE=, $
;                       QUADCUBE= , $
;                       NO_DIPOLE=, NO_MONOPOLE=, $
;                       RESO_ARCMIN= , ROT=, $
;                       SAVE=, SHADED=, SILENT=, STAGGER=, SUBTITLE=, $
;                       TITLEPLOT=, TRANSPARENT=, TRUECOLORS= $
;                       UNITS=, WINDOW=, XPOS=, YPOS=]
;                        
;  all the arguments and parameters are identical for all the
;  routines, excepted stated otherwise.
;
;
; INPUTS:
;   File = 
;          by default,           name of a FITS file containing 
;               the healpix map in an extension or in the image field
;          if Online is set :    name of a variable containing
;               the healpix map
;          if Save is set   :    name of an IDL saveset file containing
;               the healpix map stored under the variable  data
;
; OPTIONAL INPUTS:
;       Select =  if the file read is a multi column BIN table, Select indicates
;                 which column is to be plotted (the default is to plot the
;                 first column containing some signal, as opposed to pixel index)
;               can be either a name : value given in TTYPEi of the FITS file
;                        NOT case sensitive and
;                        can be truncated, 
;                        (only letters, digits and underscore are valid)
;               or an integer        : number i of the column
;                            containing the data, starting with 1
;               (see the Examples below)
;
; OPTIONAL INPUT KEYWORDS:
;
;       ASINH: if set, the color table is altered in to emulate the effect of replacing
;            the data by sinh^{-1}(data) in order to enhance the low contrast regions.
;            Can be used in conjonction with FACTOR and OFFSET, but can not be
;            used with /LOG nor /HIST_EQUAL
;
;       CHARSIZE : overall multiplicative factor applied to the size of all
;               characters appearing on the plot
;                default = 1.0
;
;       CHARTHICK : character thickness (in TITLE, SUBTITLE and color bar labeling).  
;               Other characters thickness (such as graticule labels), can be 
;               controlled with !P.CHARTHICK.
;                default = 1
;
;   COLT : color table index:
;              -Indexes [0,40] are reserved for standard IDL color tables, while
;               [41,255] are used for user defined color tables read from disc (created and
;               written to disc with MODIFYCT), if any.
;              -If the index does not match any existing table, or if it is
;              above 255, the current
;               table (modifiable with TVLCT, XLOADCT, XPALETTE, ... 
;               or eg, J.Davenport's cubehelix.pro implementation of D. Green cubehelix
;               color scheme) is used instead.
;              -If not set, the color table will be 33 (Blue-Red).
;              -If colt<0, the IDL color table ABS(colt) is used, but the scale is
;              reversed (ie a red to blue scale becomes a blue to red scale).
;              (note: -0.1 can be used as negative 0)
;
;       COORD : vector with 1 or 2 elements describing the coordinate system of the map 
;                either 'C' or 'Q' : Celestial2000 = eQuatorial,
;                       'E'        : Ecliptic,
;                       'G'        : Galactic 
;               if coord = ['x','y'] the map is rotated from system 'x' to system 'y'
;               if coord = ['y'] the map is rotated to coordinate system 'y' (with the
;               original system assumed to be Galactic unless indicated otherwise in the file)
;                  see also : Rot
;
;       CROP : if set the image file (gif, png) only contains the mollweide map and
;               not the title, color bar, ...
;               (see also : GIF, PNG)
;
;       EXECUTE: character string containing an IDL command to be executed in
;                the plotting window
;
;       FACTOR : multiplicative factor to be applied to the data (default = 1.0)
;               the data plotted is of the form FACTOR*(data + OFFSET)
;               see also : OFFSET, LOG
;
;       FITS : string containing the name of an output fits file with
;       the projected map in the primary image
;       if set to 0 or not set : no .FITS done
;       if set to 1            : output the plot in plot_XXX.fits
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;       if set to a file name  : output the plot in that file 
;    * For compatibility with standard FITS viewers (including STIFF), 
;        unobserved pixels, and pixels outside the sphere, take the value {\tt
;        NaN} (ie {\tt !values.f\_nan} in IDL).
;          * The resulting FITS file can be read in IDL with eg. map=readfits(filename). 
;          * In the case of orthographic projection, HALF_SKY must be set.
;
;       FLIP : if set, the longitude increases to the right, whereas by
;               default (astro convention) it increases towards the left
;
;       GAL_CUT: (positive float) specifies the symmetric galactic cut in degree
;             outside of which the the monopole and/or dipole fitting is done
;             (see also: NO_DIPOLE, NO_MONOPOLE)
;             ** mollview and orthview only **
;
; GIF : string containing the name of a .GIF output
;       if set to 0 or not set : no .GIF done
;       if set to 1            : output the plot in plot_XXX.gif
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;       if set to a file name  : output the plot in that file 
;             (see also : CROP, JPEG, PNG, PS and PREVIEW)
;
;       GLSIZE : character size of the graticule labels in units of CHARSIZE
;             default: 0 (ie, no labeling of graticules)
;             (see also: CHARSIZE, GRATICULE)
;
;   GRATICULE : if set, puts a graticule with delta_long = delta_lat = default
;         if graticule is set to a scalar > gmin delta_long = delta_lat = graticule
;         if set to [x,y] with x,y > gmin then delta_long = x and delta_let = y
;         ** cartview : default =  5, gmin =  0 **
;         ** gnomview : default =  5, gmin =  0 **
;         ** mollview : default = 45, gmin = 10 **
;         ** orthview : default = 45, gmin = 10 **
;
;       HALF_SKY: if set, only shows only one half of the sky 
;          (centered on (0,0) or on the location parametrized by Rot) instead of the full sky
;             ** orthview only **
;        
;       HBOUND: scalar or vector of up to 3 elements.
;          For Hbound[i]>0, overplot the boundaries of Healpix pixels
;           for the resolution parameter Nside=hbound[i].
;           The first Nside will be plotted with solid lines, 
;           the second one (if any) with dashes, 
;           the third one (if any) with dots. Obviously, better results are
;           obtained for Hbounds elements in growing order.
;           Since 0-valued boundaries are not plotted, but used for linestyle
;           assignment, providing Hbound=[0,4] (or [0,0,4]) will
;           plot Nside=4 boundaries with dashes (resp. dots), while Hbound=4 would plot the same
;           boundaries with solid lines.
;
;       HELP : if set, the routine header is printed (by doc_library)
;             and nothing else is done
;
; HIST_EQUAL : if not set, uses linear color mapping and 
;                         puts the level 0 in the middle
;                         of the color scale (ie, green for Blue-Red)
;       unless MIN and MAX are not symmetric
;               if set,     uses a histogram equalized color mapping
;     (useful for non gaussian data field)
;                     (see also : LOG)
;
;   HXSIZE: horizontal dimension (in cm) of the Hardcopy plot : Only for postscript printout
;       ** mollview : default = 26 cm ~ 10 in **
;               ** mollview : default = 15 cm         **
;       (useful for large color printer)
;               (see also : PXSIZE)
;
;       IGLSIZE : character size of the input coordinates graticule labels in units of CHARSIZE
;             default: 0 (ie, no labeling of graticules)
;             (see also: CHARSIZE, IGRATICULE)
;
;       IGRATICULE: if set, puts a graticule in the input coordinates
;          if both graticule and igraticule are set, these ones will
;          be represented by dashes
;
; JPEG : string containing the name of a (lossless) .JPEG output
;       if set to 0 or not set : no .JPEG done
;       if set to 1            : output the plot in plot_XXX.jpeg
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;       if set to a file name  : output the plot in that file 
;             (see also : CROP, GIF, PNG, PS and PREVIEW)
;
;   LOG: display the log of map (see also : HIST)
;         only applies to pixel with signal > 0.
;         (see OFFSET to offset signal)
;
;       MAP_OUT : name of the IDL variable that will contain
;         an un-altered projected map.
;         Unobserved pixels, and pixels outside the sphere take 
;       value !healpix.bad_value=-1.6375e30
;
;   MAX : max value plotted, 
;   every data point larger than MAX takes the same color as MAX
;
;   MIN : min value plotted, 
;   every data point smaller than MIN takes the same color as MIN
;
; NESTED: specify that the online file is ordered in the nested scheme
;
;   NOBAR : if set, no color bar is present
;
; NOLABELS : if set, no color bar label (min and max) is present
;
; NOPOSITION : if set, the astronomical location of the map
;         central point is not indicated
;               ** gnomview only **
;
;       NO_DIPOLE: if set to 1 (and GAL_CUT is not set) 
;                the best fit monopole *and* dipole over all valid pixels are removed
;                * if GAL_CUT is set to b>0, the best monopole and dipole fit is done on all valid
;                pixels with |galactic latitude|>b (in deg) and is removed from all
;                pixels
;             can not be used together with NO_MONOPOLE 
;             (see: GAL_CUT, NO_MONOPOLE)
;               ** mollview and orthview only **
;
;       NO_MONOPOLE: if set to 1 (and GAL_CUT is not set) 
;                the best fit monopole over all valid pixels is removed
;                * if GAL_CUT is set to b>0, the best monopole fit is done on all valid
;                pixels with |galactic latitude|>b (in deg) and is removed from all
;                pixels
;             can not be used together with NO_DIPOLE 
;             (see: GAL_CUT, NO_DIPOLE)
;               ** mollview and orthview only **
;
;       OFFSET: additive offset to apply to data (default = 0)
;               the data plotted is of the form FACTOR*(data + OFFSET)
;               can be used together with LOG
;               see also : FACTOR, LOG
;               Note : does NOT apply to polarization direction or amplitude
;               when POLARIZATION=3. Will apply to polarization amplitude when POLARIZATION=1.
;
;   ONLINE: if set, you can put directly A HEALPIX VECTOR on File (and
;       without header): useful when the vector is already
;       available on line, and avoid to have to write it on disk
;       just to be read by mollview
;   N.B. : the content of file_in is NOT altered in the
;   process
;               **  can not be used with /SAVE  **    *** OBSOLETE ***
;
;       OUTLINE : single structure, or set of structures, 
;                 each containing the coordinates of one outline to be overplotted.
;           Each structure should contain the following fields : 
;           - 'COORD' coordinate system (either, 'C', 'G', or 'E') of the contour
;           - 'RA'  or longitude coordinates (array)
;           - 'DEC' or lattitude coordinates (array of same size)
;           - 'LINE[STYLE]' : +2 : black dashes
;                           +1 : black dots
;                            0 : black solid [default]
;                           -1 : black dots on white background
;                           -2 : black dashes on white background
;           - 'PSY[M]' symbol used to represent vertices of outline
;                    (same meaning as standard PSYM in IDL,
;                     if 9<=abs(psym)<=46, D. Fanning's SYMCAT symbols 
;                     definition will be used, for example psym=9 is an open circle)
;                    if <=0, the vertices are represented with the chosen symbols, and
;                        connected, by arcs of geodesics.
;                    if >0, only the vertices are shown
;                    (default = 0)
;           - 'SYM[SIZE]' symbol size (same meaning as SYMSIZE in IDL)
;          Outline can be either a single structure, or an array of structures,
;          or a structure of structures
;
; PNG : string containing the name of a .PNG output
;       if set to 0 or not set : no .PNG done
;       if set to 1            : output the plot in plot_XXX.png
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;       if set to a file name  : output the plot in that file 
;             (see also : CROP, GIF, JPEG, PNG, PS and PREVIEW)
;
;       POLARIZATION: 
;         if set to 0, no polarization information is plotted.
;
;         otherwise, and if the input data contains polarisation information
;             (ie, Stokes parameter Q and U for each pixel)
;
;         if set to 1 
;             the AMPLITUDE P = sqrt( U^2 + Q^2 ) of the polarisation is plotted
;
;         if set to 2 
;             the ANGLE phi = 0.5*ATAN(U/Q) of the polarisation is plotted
;             Note: the angles are color coded with a fixed color table (independant of Colt)
;
;         if set to 3 or [3, scale_factor, step_factor]
;             -the temperature is color coded (with a color table defined by Colt)
;             -and the polarisation is overplot as a headless vector
;             Polarization can be a 3-element vector (the first element being 3).
;             The second element controls the average length of vectors
;             (default=1), while the third one controls the distance between
;             vectors (default=1). Non positive values are replaced by 1.
;
; PREVIEW : if set, there is a 'ghostview' preview of the postscript file (see : PS)
;                    or a 'xv' preview of the gif or png file (see: CROP, GIF,
;                    JPEG, PNG and PS)
;
; PS :  if set to 0 or not set : output to screen
;       if set to 1            : output the plot in plot_XXX.ps
;                with XXX = azimequid, cartesian, gnomic, mollweide or orthographic
;       if set to a file name  : output the plot in that file 
;               (see: CROP, GIF, JPEG, PNG and PREVIEW)
;
;   PXSIZE: number of horizontal screen_pixels / postscript_color_dots of the plot
;       ** mollview : default = 800, gnomview and cartview : default = 500 **
;       (useful for high definition color printer)
;
;   PYSIZE: number of vertical screen_pixels or postscript_color_dots of the plot
;       default = PXSIZE
;       (useful for high definition color printer)
;                ** gnomview only **
;
;       RETAIN: backing store for graphics windows in {0,1,2}. Default=2
;
;       RESO_ARCMIN: resolution of the gnomic map in arcmin
;       (default=1.5)
;                ** gnomview only **
;
;   ROT :   vector with 1, 2 or 3 elements specifing the rotation angles in DEGREE
;               to apply to the map in the 'output' coordinate system (see coord)
;             = ( lon0, [lat0, rat0]) 
;               lon0 : longitude of the point to be put at the center of the plot
;          the longitude increases Eastward, ie to the left of the plot 
;                      (unless flip is set)
;           =0 by default
;               lat0 : latitude of the point to be put at the center of the plot
;           =0 by default
;               rot0 : anti clockwise rotation to apply to the sky around the
;                     center (lon0, lat0) before projecting
;                     =0 by default
;
;   SAVE: if set, tells that the file is in IDL saveset format, 
;       the variable saved should be DATA 
;                 ** can not be used with /ONLINE **
;
;       SHADED: if set, the orthographic sphere is shaded, using a Phong model, to emulate 3D viewing.
;              The sphere is illuminated by isotropic ambiant light plus a single light source.
;                 ** Can NOT be used with GIF. **
;                   ** orthview only **
;
;       SILENT: if set, the code runs silently
;
;       STAGGER: scalar or 2 element vector.
;            - if stagger[0] is in ]0,2], 
;             3 copies of the same sphere centered at [-stagger[0], 0, stagger[0]]
;             (expressed in radius units) along the plot horizontal axis are
;             shown in ORTHOGRAPHIC projection
;             - stagger[1] (if defined), defines the angle of rotation (in degrees) applied
;               to the left and right partial spheres:
;             the lhs sphere is rotated downward by the angle provided, while the rhs one
;             is rotated upward. Rotations are swapped if FLIP is set.
;               ** orthview only **
;
;   SUBTITLE : String containing the subtitle to the plot (see TITLEPLOT)
;
;   TITLEPLOT : String containing the title of the plot, 
;         if not set the title will be File (see SUBTITLE)
;
;       TRANSPARENT: some pixels are transparent in the produced PNG file
;            if set to 1: bad pixels (usually grey) are transparent
;            if set to 2: white background pixels are transparent
;            if set to 3: all of the above
;            only valid with PNG
;
;       TRUECOLORS: if the input data is of the form [Npix,3] then the 3 fields
;            are respectively understood as {Red, Green, Blue} True Colors
;
;
; UNITS : character string containing the units, to be put on the right
;   side of the color bar (see : NOBAR)
;
;       WINDOW: IDL window index (integer)
;                 if WINDOW < 0: virtual window: no visible window opened. Can be
;               used with PNG or GIF. The Z buffer will be used instead of the 
;               X server, allowing much faster production of the image over a slow network
;                 if WINDOW in [0,31]: the specified IDL window with index WINDOW is used
;               (or reused)
;                 if WINDOW > 31: a free (=unused) window with a random index > 31 will be
;               created and used : default

; XPOS : The X position on the screen of the lower left corner
;         of the window, in device coordinate
;
; YPOS : The Y position on the screen of the lower left corner 
;               of the window, in device coordinate
;
; NOTES
;   this routine doesn't use the IDL map_set because it is precisely bugged for 
;   the mollweide projection (at least in IDL 4.0)
;
; SIDE EFFECTS
;   this routine uses ghostview when PS and PREVIEW are selected 
; or xv when GIF or PNG and PREVIEW are selected
;
; EXAMPLES
;       ;to plot the signal of the COBE-DMR 4 year map at 53 GHz
;       read_fits_sb, 'dmr_skymap_53a_4yr.fits', dmr53a, /merge  ; read it only one time
;       mollview, dmr53a, /online, 'Sig', /quad
;
;       ;to plot it in Galactic coordinate instead of Ecliptic
;       mollview, drm53a, /online, 'Sig', /quad, coord='g'
;
; COMMONS USED : view_data
;
; PROCEDURES USED: 
;       in the Healpix package :
;   index_word, read_fits_sb, vec2pix_ring, vec2pix_nest, euler_matrix
;         see  http://www.tac.dk/~healpix
;       it also requires the IDL astro library
;         http://idlastro.gsfc.nasa.gov/homepage.html
;       and the COBE analysis software
;         http://www.gsfc.nasa.gov/astro/cobe/cgis.html
;
; MODIFICATION HISTORY:
;   October 1997, Eric Hivon, TAC
;   Nov, 5, 1997,  correction of some bugs for postscript output
;   13-Nov-97, H. Dole, IAS: save and log keywords
;   4-Dec-97, H. Dole, IAS: online keyword
;   16-Dec-97, E.H, TAC: added pxsize, hxsize, subtitle, nobar
; 17-Dec-97, split the loop for projection, added nolabels
; March-98, E.H. added UNITS keyword
; April-May-98 E.H. added NESTED_ONLINE, XPOS, YPOS, NOPREVIEW
;       March-99     E.H. Caltech, improved the GIF output
;              modified to deal with structures
;              added Select, COORD, ROT, QUADCUBE  suppressed LON0
;       April-99     E.H. Caltech, improved graticule
;       Nov-99         added flip
;       Feb-00   added rmmonopole and dipole, changed common
;       March-00   changed to no_monopole and no_dipole, changed common
;       Sept-00    added polarisation plotting (Polarization)
;       June-02  : EH, Caltech. Hacked G. Giardino's polview into cartview
;       June-02    partial consolidation of gnomview/mollview/cartview
;       Jan-07    added WINDOW keyword
;       Jun-07:  edited doc header about default data to plot from cut sky file
;       Sep-07:  added /SILENT
;       Mar-08:  added GLSIZE and IGLSIZE
;       Apr-08:  can deal with cut sky data set without creating full sky map
;       Nov-08:  restore original color table and plot settings when exiting
;       May-09:  added /SHADED to orthview, implemented EXECUTE in orthview, fix
;              Min-Max for LOG, use Z buffer when window<0, added RETAIN keyword
;       Oct-09:  added /TRUECOLORS to all routines and MAP_OUT= to Gnomview
;       Apr-10:  accept array of structures in Outline; added MAP_OUT= to
;       Cartview and Mollview
;       Jan-12: added STAGGER to orthview; created azeqview; added JPEG to all
;       Jan-2013, L. Spencer: Added CTDIR, CTFILE keywords to point to separate color table.
;                             Added GRMIN, GRMAX keywords to limit the graticule labels on map edges.
;                             Added GRLS to dictate graticule linestyle.
;                             Also added IGRMIN, IGRMAX, and IGRLS to do the same for input graticule.
;-

defsysv, '!healpix', exists = exists
if (exists ne 1) then init_healpix

@viewcom ; define common
data_plot = 0 ; empty common array
; record original color table and PLOTS settings
record_original_settings, original_settings

loadsky                         ; cgis package routine, define rotation matrices
projection = 'MOLLWEIDE'
routine = 'mollview'

uroutine = strupcase(routine)
if keyword_set(help) then begin
    doc_library,'mollview'
    return
endif

if keyword_set(gif) then begin
    message_gif, code=routine, error=error_gif
    if (error_gif) then return
endif

if (n_params() lt 1 or n_params() gt 2) then begin
    PRINT, 'Wrong number of arguments in '+uroutine
    print,'Syntax : '
    print, uroutine+', File, [Select, ]'
    print,'              [ASINH=, CHARSIZE=, COLT=, COORD=, CROP=, '
    print,'              EXECUTE=, FACTOR=, FLIP=, GAL_CUT=, GIF=, GLSIZE=, GRATICULE=, '
    print,'              HBOUND=, HELP=, '
    print,'              HIST_EQUAL=, HXSIZE=,'
    print,'              IGLSIZE=, IGRATICULE=,'
    print,'              JPEG=,'
    print,'              LOG=, '
    print,'              MAP_OUT=, MAX=, MIN=, NESTED=, NOBAR=, NOLABELS=, '
    print,'              NO_DIPOLE, NO_MONOPLE, '
    print,'              OFFSET=, ONLINE=, OUTLINE=,'
    print,'              PNG=,'
    print,'              POLARIZATION=, PREVIEW=, '
    print,'              PS=, PXSIZE=, PYSIZE=, QUADCUBE= ,'
    print,'              RETAIN=, ROT=, SAVE=, SILENT=, '
    print,'              SUBTITLE=, TITLEPLOT=, TRANSPARENT=, TRUECOLORS=, '
    print,'              UNITS=, WINDOW=, XPOS=, YPOS=]'
    print
    print,' Type '+uroutine+', /help '
    print,'   for an extended help'
    return
endif

IF (undefined(file_in)) then begin
    print,routine+': Undefined variable as 1st argument'
    return
endif
do_flip = keyword_set(flip)

if (!D.n_colors lt 4) then begin
    print,' : Sorry ... not enough colors ('+strtrim(string(!d.n_colors),2)+') available'
    return
endif

if (keyword_set(no_monopole) and keyword_set(no_dipole)) then begin
    print,routine+': choose either NO_MONOPOLE or NO_DIPOLE'
    print,'    (removal of best fit monopole only or best fit monopole+dipole)'
    return
endif

polar_type = 0
if keyword_set(polarization) then polar_type = polarization

loaddata_healpix, $
  file_in, select_in,$
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, title_display, sunits, $
  SAVE=save, ONLINE=online, NESTED=nested_online, UNITS=units, COORD=coord, FLIP=flip, $
  ROT=rot, QUADCUBE=quadcube, LOG=log, ERROR=error, $
  POLARIZATION=polarization, FACTOR=factor, OFFSET=offset, SILENT=silent, COMPRESS=1, PIXEL_LIST=pixel_list, $
  TRUECOLORS=truecolors, DATA_TC=data_tc
if error NE 0 then return

data2moll, $
  data, pol_data, pix_type, pix_param, do_conv, do_rot, coord_in, coord_out, eul_mat, $
  planmap, Tmax, Tmin, color_bar, planvec, vector_scale, $
  PXSIZE=pxsize, LOG=log, HIST_EQUAL=hist_equal, MAX=max_set, MIN=min_set, FLIP=flip,  $
  NO_DIPOLE=no_dipole, NO_MONOPOLE=no_monopole, UNITS=sunits, DATA_plot = data_plot, GAL_CUT=gal_cut, $
  POLARIZATION=polarization, SILENT=silent, PIXEL_LIST=pixel_list, ASINH=asinh, $
  TRUECOLORS=truecolors, DATA_TC=data_tc, MAP_OUT = map_out, ROT=rot, FITS=fits

LS_proj2out, $
  planmap, Tmax, Tmin, color_bar, 0., title_display, $
  sunits, coord_out, do_rot, eul_mat, planvec, vector_scale, $
  CHARSIZE=charsize, COLT=colt, CROP=crop, GIF = gif, GRATICULE = graticule, bwhite=bwhite,$
  HXSIZE=hxsize, NOBAR = nobar, NOLABELS = nolabels, PNG = png, PREVIEW = preview, PS=ps, PXSIZE=pxsize, $
  SUBTITLE = subtitle, TITLEPLOT = titleplot, XPOS = xpos, YPOS = ypos, $
  POLARIZATION=polarization, OUTLINE=outline, /MOLL, FLIP=flip, COORD_IN=coord_in, IGRATICULE=igraticule, $
  HBOUND = hbound, WINDOW = window, EXECUTE=execute, SILENT=silent, GLSIZE=glsize, $
  IGLSIZE=iglsize, RETAIN=retain, TRUECOLORS=truecolors, TRANSPARENT=transparent, CHARTHICK=charthick, keep_file_open=keep_file_open,$
  JPEG=jpeg, CTDIR=CTDIR, CTFILE=CTFILE, GRMIN=GRMIN, GRMAX=GRMAX, GRLS=GRLS, IGRMIN=IGRMIN, IGRMAX=IGRMAX, IGRLS=IGRLS, CBLBL=CBLBL


w_num = !d.window
; restore original color table and PLOTS settings
record_original_settings, original_settings, /restore


return
end
