pro lk_curvatures_needlet,sig=sig,area=area,length=length,genus=genus,mftry=mftry,lmax=lmax,glpow=glpow,bb=bb,$
                          jvec=jvec,bound=bound,phiu=phiu,mpow=mpow,clvec=clvec,Bll=Bll, varmftry=varmftry,nvec=nvec,$
                          lkmask=lkmask,fsky=fsky,unnorm=unnorm,glvec=gl

  if not keyword_set(jvec) then begin
     print, 'lk_curvatures_needlet: jvec must be specified!'
     on_error,2
  endif

  if not keyword_set(mpow) then mpow=1
  if not keyword_set(glpow) then glpow=2
  if not keyword_set(lmax) then lmax=1500
  if not keyword_set(bb) then bb=2d0
  if not keyword_set(Bll) then Bll=1d0+0d0*clvec

  mf_treshold,u,sig=sig,bound=bound,nvec=nvec

  phiu=u

  nell=n_elements(jvec)
  j0=jvec[1]
  nj=nell-1

  num=n_elements(u)
  nsteps=num

  area=dblarr(num,nell)
  length=dblarr(num,nell)
  genus=dblarr(num,nell)

  mftry=dblarr(num,4,nell)

  varmftry=dblarr(num,4,nell)

  ;; get standard needlet coefficients
  print, 'getting gl using glpow, bb,j0,nj,lmax: ',glpow, bb,j0,nj,lmax
  get_needlet,bb,1,nj+j0+1,lmax,gl,/modif

  gl = gl[*,j0-1:j0+nj-2]^double(glpow)  ;;no glpow raise after this  

  for ii=0,nell-1 do begin



     if ii eq 0 then begin  ;; get full map case upto lmax

        print, 'computing theory MFs for full data upto lmax,mpow=',lmax,mpow

        ell2=2.*(2.+1.)/2. ;;not used here finally
        lk_mfs, u, ell2, lk0, lk1, lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy

        lkl_num=0d0
        lkl_dnom=0d0
        for jj=0,lmax do begin
           ell2=double(jj)*(double(jj)+1.)/2.
           coef=(2.*double(jj)+1.)/(4.*!pi)

           lkl_num = lkl_num + ell2*clvec[jj]*Bll[jj]*coef
           lkl_dnom = lkl_dnom + clvec[jj]*Bll[jj]*coef
        endfor

        lkl = lkl_num/lkl_dnom  ;; multipole dependent factor

        lk0 = lkxx + lkl*lkxy
        lk1 = sqrt(lkl)*lkyy

        gvar=lk1*0d0


     endif else begin

        jmap = jvec[ii]



        ell2=2.*(2.+1.)/2. ;will not be used
        lk_mfs, u, ell2, lk0, lk1, lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy,lkxz=lkxz

        lkl_num=0d0
        lkl_dnom=0d0
        for jj=1,lmax do begin
           ell2=double(jj)*(double(jj)+1d0)/2d0
           coef=(2d0*double(jj)+1d0)/(4d0*!pi)

           lkl_num = lkl_num + ell2*clvec[jj]*coef*gl[jj,jmap-j0]^2d0
           lkl_dnom = lkl_dnom + clvec[jj]*coef*gl[jj,jmap-j0]^2d0
        endfor

        lkl = lkl_num/lkl_dnom  ;; multipole dependent factor
        print, 'computing theory MFs for lmax, jmap, glpow, mpow, jkl, lkl=',lmax, jmap,glpow,mpow,jmap-j0,lkl

        ;print, 'lkl = ',lkl,lkl_num, lkl_dnom

        lk0 = lkxx + lkl*lkxy
        lk1 = sqrt(lkl)*lkyy

        if keyword_set(lkmask) then begin
           rho2=lkxy/(4.*!pi)
           rho1=4*lkxz/(2*!pi)
           rho0 = lkxx/2.
           
           M0 = lkmask[1]
           M1 = lkmask[2]
           M2 = lkmask[0]

           area_vec = lk2 ;;rho0*M2
           if keyword_set(fsky) then  area_vec = rho0*M2 
           length_vec = rho0*M1 + (!pi/2.)*sqrt(lkl)*rho1*M2
           genus_vec = rho0*M0 + sqrt(lkl)*rho1*M1 + lkl*rho2*M2
           
           lk0 = genus_vec
           lk1 = length_vec
           lk2 = area_vec
        endif


        gvar=lk1*0d0
     endelse

     area[*,ii]=lk2/(4*!pi)
     length[*,ii]=2.*lk1        ;/(2*!pi)
     genus[*,ii]=lk0

     mftry[*,1,ii]=area[*,ii]
     mftry[*,2,ii]=genus[*,ii]
     mftry[*,3,ii]=length[*,ii]
     
     if not keyword_set(unnorm) then begin
        if mpow eq 2 then begin
           mftry[(nsteps/2+10):nsteps-1,1,ii]=area[(nsteps/2+10):nsteps-1,ii]/max(area[(nsteps/2+10):nsteps-1,ii])
           mftry[(nsteps/2+10):nsteps-1,2,ii]=genus[(nsteps/2+10):nsteps-1,ii]/max(genus[(nsteps/2+10):nsteps-1,ii])
           mftry[(nsteps/2+10):nsteps-1,3,ii]=length[(nsteps/2+10):nsteps-1,ii]/max(length[(nsteps/2+10):nsteps-1,ii])
        endif else begin
           mftry[*,1,ii]=area[*,ii]/max(area[*,ii])
           mftry[*,2,ii]=genus[*,ii]/max(genus[*,ii])
           mftry[*,3,ii]=length[*,ii]/max(length[*,ii])
        endelse
     endif
     

     varmftry[*,2,ii]=reform(gvar)
  endfor

end
