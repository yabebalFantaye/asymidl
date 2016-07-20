;+
;+
; NAME:
; COLOURS
;
; PURPOSE:
; This procedure defines the default color table for IDL. The first 10 colors
; follow the resistor color code:
; 0 = black; 1 = brown; 2 = red; 3 = orange; 4 = yellow
;  5 = green; 6 = blue; 7 = violet; 8 = grey; 9 = white
; The rest of the colors are defined in the color table defined in colorindex.pro.
;
; CATEGORY:
; Utilities
;
; CALLING SEQUENCE:
; COLOURS
;
; INPUTS:
; None.
;
; KEYWORD PARAMETERS:
; None.
;
; OUTPUTS:
; This procedure sets the IDL color table
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;   Written by: GRD 28.11.92
; Modified Oct 2005  BGG - now uses the programrootdir() function to find the color table file.
; July 2006 (BGG) - now uses colorindex.pro to generate the color table, instead of a text file.
; 
; Copyright Brad. G. Gom 2006
;-
PRO COLOURS
  COMPILE_OPT HIDDEN
  ;
  ;basic resistor colors
  ; RED =   [0, 220,  255,  255,  255,  0,    0,    255,  160,  255]
  ; GREEN = [0, 140,  0,    127,  255,  255,  0,    0,    160,  255]
  ; BLUE =  [0, 127,  0,    0,    0,    0,    255,  255,  160,  255]

  names=['black','brown','red','orange','yellow','green','blue','violet',$
        'grey','white','dgrey','cyan','dbrown','lblue','skyblue','navyblue',$
        'iceblue','forrestgreen','lgreen','olive','lyellow','dyellow','brickred',$
        'hotpink','lpink','lviolet','lpurple','dpurple','turquoise','khaki',$
        'dorange','neon']

  ;  To get a color
  ; index for a named color, use ind=colorindex('blue') or plot,data,color=colorindex('khaki')
  ; The resistor color codes still occupy the 0-9 indices.

  rgb=colorindex(names,/rgb)

  TVLCT,rgb[0,*], rgb[1,*], rgb[2,*]

  if (!d.name eq 'WIN') || (!d.name eq 'X') || (!d.name eq 'Z') then begin
    device,decomposed=0 ;This makes IDL interpret color values as indices into a color table.
    !P.COLOR=colorindex('white')  ;set default color to white (9)
    endif else begin
    ;the other devices do not support decomposed colors.
    ;!P.COLOR=colorindex('white',/decomposed) ;set default color to white in 24-bit color
    !p.color=0  ;postscript plots color 0 as black.
    endelse

  message,'Loaded color table.',/info
  RETURN
END
