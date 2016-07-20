function smart_yra, ymin, ymax, ny=ny,ytickv=ytickv,ytickname=ytickname, chatty=chatty, format=format
;;get smart yra values


;;check ymin and ymax are not equal
if n_params() lt 2 then begin
   print, 'usage: smart_yra, ymin, ymax, ny=ny,ytickv=ytickv,ytickname=ytickname, chatty=chatty, format=format'
   return,1
endif

print, 'smart_yra: in ymin, ymax = ',ymin, ymax

if not keyword_set(ny) then ny=3
if ymin eq ymax then on_error,2

yra=[ymin,ymax]
dy = float(ymax-ymin)/float(ny)
dystr = strtrim(string(dy, format='(g0)'),2)

ndig_dy = strlen(dystr)
ndec_dy = strlen(strmid(dystr,strpos(dystr,'.')+1))

print, 'dy, ndig_dy, ndec_dy',dy, ndig_dy, ndec_dy

ysize = size(yra)
ytype= ysize[n_elements(ysize)-2]  ;;the type of ymax 


format = '(g0)'  ;;case for int [1,10000]

if dy lt 1 or ytype gt 3 then format = '(f'+strn(ndig_dy)+'.'+strn(ndec_dy)+')'  ;;case for either ymin, ymax or dy is float
;;if  abs(ymax) gt 1.1 and ytype gt 3 then format = '(f6.1)' ;;case [1.0 999.0]
;;if  abs(ymax) le 1.1 then format = '(f5.2)' ;;case [0.001 0.9]
if  abs(ymax) gt 1e3 and ytype gt 3  then  format='exponent' ;;case [1000.0, inf]
if  abs(ymax) gt 1e4  then  format='exponent' ;;case [10000.0,inf]
if  abs(ymax) lt 1e-3 then  format='exponent' ;;case [-inf, 0.001]

;;convert yra to the corresponding rounded number   

if format eq 'exponent' then begin
   yra = float(strtrim(string(yra, format='(g0)'),2))
endif else begin
   yra = float(strtrim(string(yra, format=format),2))
endelse





;;get tick vals and tick names
ytickv = fltarr(ny+1)
ytickname = strarr(ny+1)

for jj=0,ny do begin
   yval = float(yra[0] + jj*float(yra[1]-yra[0])/ny)
   ytickv[jj] = yval
   if format eq 'exponent' then begin
      ytickname[jj] = exponent(0,1,yval)
   endif else begin
      ytickname[jj] = string(yval,format=format)
   endelse
endfor

print, 'smart_yra: format',format
print, 'smart_yra: ytickname',ytickname
;;remove trailing spaces
ytickname = STRTRIM(ytickname,2)

if keyword_set(chatty) then begin
   print, 'smart_yra: ytickv=',ytickv
   print, 'smart_yra: ytickname: ',ytickname
endif

yra = [min(ytickv),max(ytickv)]


return, yra

end



;; if keyword_set(symmetry) then ymin=-ymax ;;set symmetric

;; coef_ymin=1
;; coef_ymax=1

;; ndigit=10 ;;two digits
;; if keyword_set(digit) then ndigit=10l^(digit-1)


;; ;;if value is [-1 1]
;; if (abs(ymin) lt ndigit) and (abs(ymax) lt ndigit)  then begin

;;    coef_ymin=1l
;;    ymin_new=ymin
;;    while abs(ymin_new) lt ndigit) and do begin
;;       ymin_new = ymin_new*10l
;;       coef_ymin=coef_ymin*10l
;;    endwhile

;;    if keyword_set(chatty) then  print, 'coef_ymin',coef_ymin,ymin_new
   
;;    coef_ymax=1
;;    ymax_new=ymax
;;    while abs(ymax_new) lt 1 do begin
;;       ymax_new = ymax_new*10l
;;       coef_ymax=coef_ymax*10l
;;    endwhile
   
;;    if keyword_set(chatty) then print, 'coef_ymax',coef_ymax,ymax_new

;; endif


;; ;;get the first 10th component
;; coef = min([coef_ymin,coef_ymax])


;; ;;exponent part when value is
;; dexpo=''
;; n10 = round(alog10(coef))

;; ;;format = '(f'+strn(n10+3)+'.'+strn(n10)+')' ;;0.002 will have format '(f5.3)'
