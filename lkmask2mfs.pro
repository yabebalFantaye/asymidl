function lkmask2mfs,u,lkmask,ell=ell,mpow=mpow,bb=bb,jmap=jmap,nj=nj,lmax=lmax,glpow=glpow,clvec=clvec,fsky=fsky


  if n_params() lt 1 then begin
     print, 'maskmfs: u must be specified!'
     on_error,2
  endif

  if not keyword_set(lkmask) then begin
     print, 'maskmfs: maskmf must be specified!'
     on_error,2
  endif


  lkl=2.*(2.+1)/!pi ;; not used
  if keyword_set(ell) then lkl=double(ell)*(double(ell)+1.)/2.

     lk_mfs, u, lkl, lk0, lk1, lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy,lkxz=lkxz


  if not keyword_set(ell) then begin 



     if not keyword_set(glpow) then glpow=print_set(2,'mask2mfs: setting glpow=2')
     if not keyword_set(lmax) then lmax=print_set(1500,'mask2mfs: setting lmax=1500')
     if not keyword_set(bb) then bb=print_set(2d0,'mask2mfs: setting bb=2.')
     if not keyword_set(jmap) then jmap=print_set(2,'mask2mfs: setting jmap=2')
     
     ;; get standard needlet coefficients
     get_needlet,bb,1,nj,lmax,gl
     
     
     lkl_num=0d0
     lkl_dnom=0d0
     for jj=0,lmax do begin
        ell2=double(jj)*(double(jj)+1.)/2.
        coef=(2.*double(jj)+1.)/(4.*!pi)
        
        lkl_num = lkl_num + ell2*clvec[jj]*coef*gl[jj,jmap-1]^(2*glpow)
        lkl_dnom = lkl_dnom + clvec[jj]*coef*gl[jj,jmap-1]^(2*glpow)
     endfor
     lkl = lkl_num/lkl_dnom  ;; multipole dependent factor

  endif

  rho2=lkxy/(4.*!pi)
  rho1=4*lkxz/(2*!pi)
  rho0 = lkxx/2.
  
  M0 = lkmask[1]
  M1 = lkmask[2]
  M2 = lkmask[0]
  

  area_vec = lk2 ;;rho0*M2
  if keyword_set(fsky) then   area_vec = rho0*M2
  length_vec = rho0*M1 + (!pi/2.)*sqrt(lkl)*rho1*M2
  genus_vec = rho0*M0 + sqrt(lkl)*rho1*M1 + lkl*rho2*M2

  
  lk0 = genus_vec
  lk1 = length_vec
  lk2 = area_vec

  area=lk2/(4*!pi)
  length=2.*lk1           ;/(2*!pi)
  genus=lk0
  
  mfs=[[area], [genus], [length]]

return, mfs

end
