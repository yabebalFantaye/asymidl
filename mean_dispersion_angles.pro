pro  mean_dispersion_angles,thx, phx,mm,significance,outstruct=outstruct,nrep=nrep,iout=iout

if not keyword_set(nrep) then nrep=1
if not keyword_set(iout) then iout=0

nsim = n_elements(thx[0,*])
nband = n_elements(thx[*,0])

significance=dblarr(nband-1,nband-1)
mm=dblarr(nsim,nband-1,nband-1)

cor=dblarr(nband,nband)
dcor=cor
allcor=dblarr(nband,nband,nsim)


for i=0,nband-1 do begin
   for j=i,nband-1 do begin
      for k=0l,nsim-1 do begin
         ANG2VEC, thx(i,k), phx(i,k), vec1
         ANG2VEC, thx(j,k), phx(j,k), vec2
         ;if (i ne j) then cor(i,j)=cor(i,j)+acos(FLOAT(total(vec1*vec2)))*180d0/!pi
         ;if (i ne j) then dcor(i,j)=dcor(i,j)+(acos(FLOAT(total(vec1*vec2)))*180d0/!pi)^2
         if (i ne j) then allcor(i,j,k)=acos(FLOAT(total(vec1*vec2)))*180d0/!pi
      endfor
   endfor
endfor
print, 'for loops for sep angle and disp done'

;cor=cor/double(nsim)
;dcor=dcor/double(nsim) ;; error^2 = sum(err^2/n^2)
;dcor=sqrt(dcor-cor^2)

;;set one to lower diagonal matrix and 0 to diagonal and upper
;;diagonal parts

count = 0l

count2=0
for lastb=14,1,-1 do begin
   


   firstb_last = 0   
   if (lastb eq 14) then firstb_last = lastb-1

   count1 = 0            
   for firstb=0,firstb_last do begin
      
      v=dblarr(nband,nband)            
      for i=firstb,lastb do begin
         for j=i+1,lastb do begin
            v(i,j)=1.
         endfor
      endfor
      
      for i=0l,nsim-1 do begin
         mm(i,count1,count2)=total((allcor(*,*,i))*v)/total(v)
      endfor
      
      mndat = mm[0,count1,count2]
      mnangle = mm[1:nsim-1,count1,count2]
      significance[count1,count2] = max([total(mnangle lt mndat),0]*100d0/((nsim-1)*1d0))
      
      count1 = count1+1
   endfor

   count2 = count2+1
endfor

if keyword_set(outstruct) then begin

   if size(outstruct,/type) ne 8 then begin
      outstruct={sig:significance,mangle:mm}
      outstruct = replicate(outstruct,nrep)
      iout=0
   endif else begin
      outstruct[iout].sig = significance
      outstruct[iout].mangle = mm
   endelse

endif

end



   ;; mn=dblarr(nsim)
   ;; cnt=0d0
   ;; imin=1
   ;; imax=nband-1
   ;; jmin=0
   ;; jmax=nband-1

   ;; for i=imin,imax do begin
   ;;    jmax=i-1
   ;;    for j=jmin,jmax do begin
   ;;       mn=mn+acos(cos(thx(i,*))*cos(thx(j,*))+sin(thx(i,*))*sin(thx(j,*))*cos(phx(i,*)-phx(j,*)))*180d0/!Pi
   ;;       cnt=cnt+1d0
   ;;    endfor
   ;; endfor
   ;; mn=mn/cnt



