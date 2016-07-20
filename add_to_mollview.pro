pro add_to_mollview,xep=xep,range=range

;;;insert this in proj2out before running mollview:
;xyouts,1.41,-0.29,'*',col=1,charsize=5,charthick=40 ;;;position found by eye
;xyouts,1.65,-0.15,'*',col=0,charsize=5,charthick=40 ;;;position found by eye

if keyword_set(xep) then begin
    xyouts,-1.02,0.365,'X(NEP)',charsize=2,col=0,charthick=1 ;;;position found by eye
    xyouts,0.82,-0.445,'X(SEP)',charsize=2,col=0,charthick=1 ;;;position found by eye
endif


cbar_dx = 2./3.
cbar_dy = 1./32.
cbar_xll = (1. - cbar_dx)/2.
cbar_xur = (1. + cbar_dx)/2.
cbar_yur = w_yll - cbar_dy
cbar_yll = cbar_yur - cbar_dy

if keyword_set(range) then begin
    dx_cbar = round((range[1]-range[0])/5l)
    for iijj=1,5 do begin
        cgtext,cbar_xll+(cbar_dx/5.)*iijj,-1.1,'.',col=0,charsize=1,charthick=2 ;;-0.7->0.7
        cgtext,cbar_xll+(cbar_dx/5.)*iijj-0.05,-1.2,strtrim(string(iijj*dx_cbar),2),charsize=1
    endfor
endif

end
