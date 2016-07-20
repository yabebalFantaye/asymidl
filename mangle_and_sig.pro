pro  mangle_and_sig,thx, phx,mm,significance,outstruct=outstruct,nrep=nrep,iout=iout,thdp=thdp,phdp=phdp,zstart=zstart,rstart=rstart,median=meddian,nnb=nnb

if not keyword_set(zstart) then zstart=1
if not keyword_set(nrep) then nrep=1
if not keyword_set(iout) then iout=0

nsim = n_elements(thx[0,*])
nband = n_elements(thx[*,0])

significance=dblarr(nband-1,4,1)
if zstart gt 1 then significance=dblarr(nband-1,4,nsim-zstart)
mm=dblarr(nsim,nband-1, 4)


;;for mm and significance indices 0 and 1 are for lmax and lmin
;;changes; 2 and 3 are the same but for powel method
;; 

;; if not keyword_set(thdp) then thdp = 41d0*!dtor
;; if not keyword_set(phdp) then phdp = 264d0*!dtor

if not keyword_set(nnb) then nnb=nband-1

mthph= mean_healpix_ang(thx[0:nnb,*],phx[0:nnb,*],median=median)
help, mthph
print, 'mean angle until nnb=',nnb
print,mthph[*,2]*!rtod

;mthph=dblarr(2,nsim)
;mthph[0,*]=replicate(89d0*!dtor,nsim)
;mthph[1,*]=replicate(225d0*!dtor,nsim)
if not keyword_set(thdp) then thdp = mthph[0,*];89d0*!dtor
if not keyword_set(phdp) then phdp = mthph[1,*] ;225d0*!dtor


;;significance as a function of lmax
count1 = 0
for imax=1,nband-1 do begin

   mn=dblarr(nsim)
;;sum over cos(theta) 
   mnpowel=dblarr(nsim)
   
   
   imin=1
   jmin=0

  ;;the first bin wrt to (phdp,thdp)
   mnpowel=reform(acos(cos(thx(0,*))*cos(thdp[0,*])+sin(thx(0,*))*sin(thdp[0,*])*cos(phx(0,*)-phdp[0,*])))
   
   cnt=0d0
   cnt2=1d0
   for i=imin,imax do begin
      jmax=i-1
      for j=jmin,jmax do begin
         mn=mn+reform(acos(cos(thx(i,*))*cos(thx(j,*))+sin(thx(i,*))*sin(thx(j,*))*cos(phx(i,*)-phx(j,*)))*180d0/!Pi)
         cnt=cnt+1d0
         if j eq jmin then begin
            mnpowel=mnpowel+reform(acos(cos(thx[i,*])*cos(thdp[0,*])+sin(thx[i,*])*sin(thdp[0,*])*cos(phx[i,*]-phdp[0,*]))*180d0/!Pi)
            cnt2=cnt2+1d0
         endif
      endfor
   endfor
   mn=mn/cnt
   mnpowel=mnpowel/cnt2


;;reshufle such that the rstart index becomes the zeroth index
   if keyword_set(rstart) then begin
      if rstart gt 0 then begin
         temp = mn[0:rstart-1]
         mn[0]=mn[rstart]
         mn[1:rstart] = temp
         
         temp = mnpowel[0:rstart-1]
         mnpowel[0]=mnpowel[rstart]
         mnpowel[1:rstart] = temp
      endif
   endif

   mnangle = mn[zstart:nsim-1]
   mnangle_powel = mnpowel[zstart:nsim-1]

   count2 = 0
   for isim=0,zstart-1 do begin
      mndat = mn[isim]
      mndat_powel = mnpowel[isim]

      significance[count1, 0,count2] = total(mnangle lt mndat) ;*100d0/((nsim-zstart)*1d0)
      significance[count1, 2,count2] = min([total(mnangle_powel lt mndat_powel),total(mnangle_powel gt mndat_powel)]) ; *100d0/((nsim-zstart)*1d0)
      count2 = count2+1
   endfor
   mm[*,count1, 0] = mn
   mm[*,count1, 2] = mnpowel


   count1 = count1+1
endfor


;;significance as a function of lmin
count1 = 0
for imin=0,nband-2 do begin

   mn=dblarr(nsim)
;;sum over cos(theta) 
   mnpowel=dblarr(nsim)
   
   imax=nband-1
   jmin=imin

  ;;the imin bin wrt to (phdp,thdp)
   mnpowel=reform(acos(cos(thx(imin,*))*cos(thdp)+sin(thx(imin,*))*sin(thdp)*cos(phx(imin,*)-phdp)))

   cnt=0d0
   cnt2=1d0
   for i=imin+1,imax do begin
      jmax=i-1
      for j=jmin,jmax do begin
         mn=mn+reform(acos(cos(thx(i,*))*cos(thx(j,*))+sin(thx(i,*))*sin(thx(j,*))*cos(phx(i,*)-phx(j,*)))*180d0/!Pi)
         cnt=cnt+1d0
         if j eq jmin then begin
            mnpowel=mnpowel+reform(acos(cos(thx(i,*))*cos(thdp)+sin(thx(i,*))*sin(thdp)*cos(phx(i,*)-phdp))*180d0/!Pi)
            cnt2=cnt2+1d0
         endif
      endfor
   endfor
   mn=mn/cnt
   mnpowel=mnpowel/cnt2
   

;;reshufle such that the rstart index becomes the zeroth index
   if keyword_set(rstart) then begin
      if rstart gt 0 then begin
         temp = mn[0:rstart-1]
         mn[0]=mn[rstart]
         mn[1:rstart] = temp
         
         temp = mnpowel[0:rstart-1]
         mnpowel[0]=mnpowel[rstart]
         mnpowel[1:rstart] = temp
      endif
   endif


   mnangle = mn[zstart:nsim-1]
   mnangle_powel = mnpowel[zstart:nsim-1]

   count2 = 0
   for isim=0,zstart-1 do begin
      mndat = mn[isim]
      mndat_powel = mnpowel[isim]

      significance[count1, 1,count2] = total(mnangle lt mndat) ;,0]*100d0/((nsim-zstart)*1d0))
      significance[count1, 3,count2] = min([total(mnangle_powel lt mndat_powel),total(mnangle_powel gt mndat_powel)]) ;*100d0/((nsim-zstart)*1d0))
      count2 = count2+1
   endfor
   mm[*,count1, 1] = mn
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






