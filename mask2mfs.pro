function mask2mfs,u,maskmf,ell=ell,mpow=mpow,bb=bb,jmap=jmap,jstart=jstart,nj=nj,lmax=lmax,glpow=glpow,clvec=clvec,fsky=fsky,mlen=mlen


  if n_params() lt 2 then begin
     print, 'maskmfs: u must be specified!'
     on_error,2
  endif

  if not keyword_set(maskmf) then begin
     print, 'maskmfs: maskmf must be specified!'
     on_error,2
  endif


  lkl=2.*(2.+1)/2. ;; not used
  if keyword_set(ell) then lkl=double(ell)*(double(ell)+1.)/2.

  lk_mfs, u, lkl, lk0, lk1, lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy,lkxz=lkxz
    
  if not keyword_set(ell) then begin 

     

     if not keyword_set(glpow) then glpow=print_set(2,'mask2mfs: setting glpow=2')
     if not keyword_set(lmax) then lmax=print_set(1500,'mask2mfs: setting lmax=1500')
     if not keyword_set(jstart) then jstart=print_set(1d0,'mask2mfs: setting jstart=1')
     if not keyword_set(bb) then bb=print_set(2d0,'mask2mfs: setting bb=2.')
     if not keyword_set(jmap) then jmap=print_set(2,'mask2mfs: setting jmap=2')
     
     ;; get standard needlet coefficients
     print, 'getting needlet gl: bb,jstart,nj,lmax',bb,jstart,nj,lmax
     get_needlet,bb,jstart,nj,lmax,gl
     
     
     lkl_num=0d0
     lkl_dnom=0d0
     for jj=0,lmax do begin
        ell2=double(jj)*(double(jj)+1.)/2.
        coef=(2.*double(jj)+1.)/(4.*!pi)
        
        lkl_num = lkl_num + ell2*clvec[jj]*coef*gl[jj,jmap-jstart]^(2*glpow)
        lkl_dnom = lkl_dnom + clvec[jj]*coef*gl[jj,jmap-jstart]^(2*glpow)
     endfor
     lkl = lkl_num/lkl_dnom  ;; multipole dependent factor

  endif

;;note that L0=genus, L1=length and L2=area
;; maskmf is ordered L2 L0 L1 - so we write equations as
;; area_mask = lkxx*L2 + 0*L0 + 0*L1
;; genus_mask=lkl*lkxy/(4.*!pi)*L2+lkxx/2.*L0+sqrt(lkl)*lkyy/(4.*!pi)*L1
;; length_mask = sqrt(lkl)*lkyy/4.L2 + 0*L0 + lkxx*L1



estlk=dblarr(3)
nu=n_elements(u)

rho2=lkxy/(4.*!pi)
rho1=4*lkxz/(2*!pi)
rho0 = lkxx/2.


igenus=1
ilength=2
iarea=0

;print, 'maskmf: ',reform(maskmf)
;print, 'nu to use: ',nu

;; print, 'lkl ', lkl

;; if  not (keyword_set(fsky) and n_elements(mlen) gt 0) then begin

;;    for ii=0,nu-1 do begin
;;       l0 = maskmf[ii,igenus]
;;       l1 = maskmf[ii,ilength]/2.
;;       l2 = maskmf[ii,iarea]*(4.*!pi)
      
;;       M2 = (l2/rho0[ii])
;;       if keyword_set(fsky) then M2=fsky*4.*!pi
      
;;       M1 = (l1-(!pi/2.)*sqrt(lkl)*rho1[ii]*M2)/(rho0[ii])
;;       if keyword_set(mlen)  then M1=mlen/2.
      
;;       M0 = (l0-(sqrt(lkl)*rho1[ii]*M1 + lkl*rho2[ii]*M2))/rho0[ii]
      
;;                                 ;print, 'maskmfs: l0,l1,l2 = ',[estmf0,estmf1,estmf2]
      
;;       estlk = estlk + [M2,M0,M1]/float(nu) ;;take the mean
;;    endfor
   
;;    l2l0l1=estlk

;; endif else begin
   

   


   ;; M2=fsky*4.*!pi
   ;; M1=mlen/2.


   ;; Bvec = reform((l0-(sqrt(lkl)*rho1*M1 + lkl*rho2*M2)))
   ;; Amat = transpose(rho0)
   ;; M0  = LA_LEAST_SQUARES(Amat, Bvec)
   ;; l2l0l1 = [M2,M0,M1]


   l1 = maskmf[*,ilength]/2.
   Amat = dblarr(2,n_elements(u))
   Bvec = l1
   Amat[0,*] = transpose(rho0) ;;M1 coef
   Amat[1,*]=(!pi/2.)*sqrt(lkl)*rho1  ;;M2 coef
   M1M2 = LA_LEAST_SQUARES(Amat, Bvec)


;;----------
   l0 = maskmf[*,igenus]
   Amat = dblarr(1,n_elements(u))
   Bvec = l0 - (sqrt(lkl)*rho1*M1M2[0] + lkl*rho2*M1M2[1])
   Amat[0,*] = transpose(rho0) ;;M0 coef
   M0 = LA_LEAST_SQUARES(Amat, Bvec)

   l2l0l1_lg = [M1M2[1],M0,M1M2[0]]
   print, 'l2l0l1 from genus+length fit = ',l2l0l1_lg

   l0 = maskmf[*,igenus]
   Amat = dblarr(3,n_elements(u))
   Bvec = l0 
   Amat[0,*] = transpose(rho0) ;;M0 coef
   Amat[1,*]=sqrt(lkl)*rho1 ;;M1 coef
   Amat[2,*]=lkl*rho2 ;;M2 coef
   M0M1M2 = LA_LEAST_SQUARES(Amat, Bvec)
   l2l0l1 = [M0M1M2[2],M0M1M2[0],M0M1M2[1]]
   print, 'l2l0l1 from genus only fit = ',l2l0l1

;endelse

return, l2l0l1_lg

end




;; Amat = dblarr(3,n_elements(u)*3)
;; Bvec = dblarr(n_elements(u)*3)


;;   area_vec = [lk2[ii], 0, 0]
;;   genus_vec = [lkl*lkxy[ii]/(4.*!pi), lkxx[ii]/2., sqrt(lkl)*lkxz[ii]]
;;   length_vec = [sqrt(lkl)*lkxz[ii],0,lkxx[ii]]

;;   a = [[area_vec],$
;;        [genus_vec],$
;;        [length_vec]]

;;   Amat[*,3*ii:(3*(ii+1)-1)]=a

  
;;   
;;   Bvec[3*ii:(3*(ii+1)-1)] = lkmask


;; endfor
;; l2l0l1  = LA_LEAST_SQUARES(Amat, Bvec)



;  print, 'theory lks: ',[[lk2],[lk0],[lk1]]
;  print, 'sim lks: ',Bvec
  
  
  ;; print, '------------'
  ;; print, 'AX=B: A=',a
  ;; print, 'B = ',lkmask
  ;; print, '------------'
  
;;  l2l0l1 = LA_LINEAR_EQUATION(Amat, Bvec)  ;;square matrix only

