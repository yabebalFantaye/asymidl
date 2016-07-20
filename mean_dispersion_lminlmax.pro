pro  mean_dispersion_lminlmax,thx, phx,mm,significance,outstruct=outstruct,nrep=nrep,iout=iout,thdp=thdp,phdp=phdp,zstart=zstart

if not keyword_set(zstart) then zstart=1
if not keyword_set(nrep) then nrep=1
if not keyword_set(iout) then iout=0

nsim = n_elements(thx[0,*])
nband = n_elements(thx[*,0])


significance=dblarr(nband-1,4,nsim-zstart)
mm=dblarr(nsim,nband-1, 4)



;;for mm and significance indices 0 and 1 are for lmax and lmin
;;changes; 2 and 3 are the same but for powel method
;; 

if not keyword_set(thdp) then thdp = 41d0*!dtor
if not keyword_set(phdp) then phdp = 264d0*!dtor


cor=dblarr(nband,nband)
dcor=cor
allcor=dblarr(nband,nband,nsim)

;;===========================
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

;;a function of lmax 
v=dblarr(nband,nband)            
count1 = 0            
for lastb=1,nband-1,1 do begin

   for i=0,lastb do begin
      for j=i+1,lastb do begin
         v(i,j)=1.
      endfor
   endfor
   
   for i=0l,nsim-1 do begin
      mm[i,count1,0]=total((allcor(*,*,i))*v)/total(v)
   endfor
   
   mnangle = mm[zstart:nsim-1,count1,0]   

   count2 = 0
   for isim=0,zstart-1 do begin
      mndat = mm[isim,count1,0]      
      significance[count1, 0,count2] = total(mnangle lt mndat) ;*100d0/((nsim-zstart)*1d0)
      count2 = count2+1
   endfor
   
   count1 = count1+1
endfor


;;a function of lmin
v=dblarr(nband,nband)            
count1 = 0            
for firstb=0,nband-2,1 do begin

   for i=firstb,nband-1 do begin
      for j=i+1,nband-1 do begin
         v(i,j)=1.
      endfor
   endfor
   
   for i=0l,nsim-1 do begin
      mm[i,count1,1]=total((allcor(*,*,i))*v)/total(v)
   endfor

   mnangle = mm[zstart:nsim-1,count1,1]   

   count2 = 0
   for isim=0,zstart-1 do begin
      mndat = mm[isim,count1,1]
      
      significance[count1, 1,count2] = total(mnangle lt mndat) ;*100d0/((nsim-zstart)*1d0)
      count2 = count2+1
   endfor

   count1 = count1+1
endfor


;;=================================

;;significance as a function of lmax powel statistics
count1 = 0
for imax=1,nband-1 do begin

;;sum over cos(theta) 
   mnpowel=dblarr(nsim)
   
   imin=1
   jmin=0

  ;;the first bin wrt to (phdp,thdp)
   mnpowel=(cos(thx(0,*))*cos(thdp)+sin(thx(0,*))*sin(thdp)*cos(phx(0,*)-phdp))
   

   cnt2=1d0
   for i=imin,imax do begin
      mnpowel=mnpowel+(cos(thx(i,*))*cos(thdp)+sin(thx(i,*))*sin(thdp)*cos(phx(i,*)-phdp))
      cnt2=cnt2+1d0
   endfor
   mnpowel=mnpowel/cnt2
   
   mnangle_powel = mnpowel[zstart:nsim-1]

   count2 = 0
   for isim=0,zstart-1 do begin
      mndat_powel = mnpowel[isim]
      significance[count1, 2,count2] = min([total(mnangle_powel lt mndat_powel),total(mnangle_powel gt mndat_powel)]) ; *100d0/((nsim-zstart)*1d0)
      count2 = count2+1
   endfor
   mm[*,count1, 2] = mnpowel

   count1 = count1+1
endfor


;;significance as a function of lmin
count1 = 0
for imin=0,nband-2 do begin

   ;;sum over cos(theta) 
   mnpowel=dblarr(nsim)
   
   imax=nband-1
   jmin=imin

  ;;the imin bin wrt to (phdp,thdp)
   mnpowel=(cos(thx(imin,*))*cos(thdp)+sin(thx(imin,*))*sin(thdp)*cos(phx(imin,*)-phdp))

   cnt=0d0
   cnt2=1d0
   for i=imin+1,imax do begin
      mnpowel=mnpowel+(cos(thx(i,*))*cos(thdp)+sin(thx(i,*))*sin(thdp)*cos(phx(i,*)-phdp))
      cnt2=cnt2+1d0
   endfor
   mnpowel=mnpowel/cnt2
   
   mnangle_powel = mnpowel[zstart:nsim-1]

   count2 = 0
   for isim=0,zstart-1 do begin
      mndat_powel = mnpowel[isim]
      significance[count1, 3,count2] = min([total(mnangle_powel lt mndat_powel),total(mnangle_powel gt mndat_powel)]) ; *100d0/((nsim-zstart)*1d0)
      count2 = count2+1
   endfor
   mm[*,count1, 3] = mnpowel

   count1 = count1+1
endfor

print, 'mean_dispersion_lminlmax done!'

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






