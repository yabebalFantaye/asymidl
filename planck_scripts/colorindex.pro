; NAME:
; COLORINDEX
;
; PURPOSE:
; This function returns the color index or RGB triplet for the given color name
;
; CATEGORY:
; Utilities
;
; CALLING SEQUENCE:
; Result = COLORINDEX(Name)
; (i.e., NO KEYWORDS)
;
; INPUTS:
; Name: The name(s) of the color(s) to return. This name must exist in the array of colors.
;
; KEYWORD PARAMETERS:
; ALL:  Set this keyword to return the full table of values
; RGB:  Set this keyword to return an RGB triple instead of an index
; COUNT: set this keyword to a named variable to return the total number of entries in the table.
;
; OUTPUTS:
; This function returns an index into the color table or an RGB value.
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;   Written by: Brad Gom 2000
; July 5 2006 (BGG) - no longer uses external text file. Colors are defined in-line.
; Copyright Brad G. Gom 2006
;-
function colorindex,name,rgb=rgb,count=count,decomposed=decomposed,all=all
  COMPILE_OPT HIDDEN

  names=['black','brown','red','orange','yellow','green','blue','violet','grey',$
    'white','dgrey','cyan','dbrown','lblue','skyblue','navyblue','iceblue',$
    'forrestgreen','lgreen','olive','lyellow','dyellow','brickred','hotpink',$
    'lpink','lviolet','lpurple','dpurple','turquoise','khaki','dorange','neon',$
    'lime','tan','coral','teal','indigo','ivory','gold']

  rgb_table=[[0,0,0],$
    [220,140,127],$
    [255,0,0],$
    [255,127,0],$
    [255,255,0],$
    [0,255,0],$
    [0,0,255],$
    [255,0,255],$
    [198,198,198],$
    [255,255,255],$
    [100,100,100],$
    [0,255,255],$
    [102,51,51],$
    [0,0,200],$
    [0,204,255],$
    [0,51,153],$
    [153,255,255],$
    [0,102,51],$
    [51,204,102],$
    [153,153,51],$
    [255,255,204],$
    [255,204,0],$
    [204,51,0],$
    [255,51,153],$
    [255,204,204],$
    [255,153,255],$
    [204,102,255],$
    [51,0,102],$
    [102,255,204],$
    [153,153,102],$
    [255,102,51],$
    [255,0,102],$
    [50,205,50],$
    [210,180,140],$
    [255,127,80],$
    [0,128,128],$
    [75,0,130],$
    [255,255,240],$
    [255,215,0] ]

  count=n_elements(names)

  n=n_elements(name)

  if n eq 0 then ind=0 else ind=intarr(n)

  for i=0,n-1 do begin
    ind[i]=where(names eq name[i], count)
    if count eq 0 then begin
      message,'Color '+name[i]+' does not exist!',/cont
      ind[i]=0
      endif
    endfor

  if n eq 1 then ret=ind[0] else ret=ind

  if keyword_Set(all) then ind=indgen(count)  ;return all values

  if keyword_set(rgb) then $ ;return a RGB triple
    return,rgb_table[*,ind]
  if keyword_set(decomposed) then begin ;return a 24 bit color value
    rgb_table=long(rgb_table)
    return,rgb_table[0,ind]+ishft(rgb_table[1,ind],8)+ishft(rgb_table[2,ind],16)
    endif
  return,ret ;return a color index
end
