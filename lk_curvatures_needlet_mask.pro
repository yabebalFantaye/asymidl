pro lk_curvatures_needlet_mask,sig=sig,area=area,length=length,genus=genus,mftry=mftry,lmax=lmax,glpow=glpow,bb=bb,$
                          jvec=jvec,bound=bound,phiu=phiu,mpow=mpow,clvec=clvec,varmftry=varmftry,nvec=nvec,lkmask=lkmask

  if not keyword_set(jvec) then begin
     print, 'lk_curvatures_needlet_mask: jvec must be specified!'
     on_error,2
  endif

  if not keyword_set(lkmask) then begin
     print, 'lk_curvatures_needlet_mask: lkmask must be specified!'
     on_error,2
  endif

  if not keyword_set(mpow) then mpow=1
  if not keyword_set(glpow) then glpow=2
  if not keyword_set(lmax) then lmax=1500
  if not keyword_set(bb) then bb=2d0

  if not keyword_set(phiu) then begin
     mf_treshold,u,sig=sig,bound=bound,nvec=nvec
     phiu=u
  endif else u=phiu

  nell=n_elements(jvec)
  j0=jvec[1]
  nj=nell-1

  num=n_elements(u)
  area=dblarr(num,nell)
  length=dblarr(num,nell)
  genus=dblarr(num,nell)

  mftry=dblarr(num,4,nell)

  varmftry=dblarr(num,4,nell)

  ;; get standard needlet coefficients
  get_needlet,bb,j0,nj,lmax,gl
  

  for ii=0,nell-1 do begin



     if ii eq 0 then begin  ;; get full map case upto lmax

        print, 'computing theory MFs for full data upto lmax,mpow=',lmax,mpow

        ell2=2.*(2.+1.)/2. ;;not used here finally
        lk_mfs, u, ell2, lk0, lk1, lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy,lkxz=lkxz

        lkl_num=0d0
        lkl_dnom=0d0
        for jj=0,lmax do begin
           ell2=double(jj)*(double(jj)+1.)/2.
           coef=(2.*double(jj)+1.)/(4.*!pi)

           lkl_num = lkl_num + ell2*clvec[jj]*coef
           lkl_dnom = lkl_dnom + clvec[jj]*coef
        endfor

        lkl = lkl_num/lkl_dnom  ;; multipole dependent factor

        lk0 = lkxx + lkl*lkxy
        lk1 = sqrt(lkl)*lkyy

        gvar=lk1*0d0


     endif else begin

        print, 'computing theory MFs for jmap, glpow, mpow=',jvec[ii],glpow,mpow

        jmap = jvec[ii]


        ell2=2.*(2.+1.)/2. ;will not be used
        lk_mfs, u, ell2, lk0, lk1, lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy,lkxz=lkxz

        lkl_num=0d0
        lkl_dnom=0d0
        for jj=0,lmax do begin
           ell2=double(jj)*(double(jj)+1.)/2.
           coef=(2.*double(jj)+1.)/(4.*!pi)

           lkl_num = lkl_num + ell2*clvec[jj]*coef*gl[jj,jmap-j0]^(2*glpow)
           lkl_dnom = lkl_dnom + clvec[jj]*coef*gl[jj,jmap-j0]^(2*glpow)
        endfor

        lkl = lkl_num/lkl_dnom  ;; multipole dependent factor

        area_vec = lkmask[0]*lk2
        genus_vec = lkmask[0]*lkl*lkxy/(4.*!pi) + lkmask[1]*lkxx/2. + lkmask[2]*sqrt(lkl)*lkxz
        length_vec = lkmask[0]*sqrt(lkl)*lkxz + lkmask[2]*lkxx
        
        lk0 = genus_vec
        lk1 = length_vec
        lk2 = area_vec


        gvar=lk1*0d0
     endelse

     area[*,ii]=lk2/(4*!pi)
     length[*,ii]=2.*lk1        ;/(2*!pi)
     genus[*,ii]=lk0

     mftry[*,1,ii]=area[*,ii]
     mftry[*,2,ii]=genus[*,ii]
     mftry[*,3,ii]=length[*,ii]


     varmftry[*,2,ii]=reform(gvar)
  endfor

end
