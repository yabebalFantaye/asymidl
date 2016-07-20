pro lk_curvatures,sig=sig,area=area,length=length,genus=genus,mftry=mftry,ell=ell,bound=bound,$
                  phiu=phiu,mpow=mpow,clvec=clvec,varmftry=varmftry,nvec=nvec,lkmask=lkmask,unnorm=unnorm

  if not keyword_set(ell) then ell=100
  if not keyword_set(clvec) then clvec=dblarr(ell+1)+1d0
  if not keyword_set(mpow) then mpow=1
  if not keyword_set(Bll) then Bll=1d0+0d0*clvec

if not keyword_set(phiu) then begin
   print, 'lk_curvatures: computing phiu'
   mf_treshold,u,sig=sig,bound=bound,nvec=nvec
   phiu=u
endif else u=phiu

  nell=n_elements(ell)
  num=n_elements(u)
  area=dblarr(num,nell)
  length=dblarr(num,nell)
  genus=dblarr(num,nell)

  mftry=dblarr(num,4,nell)

  varmftry=dblarr(num,4,nell)

  

  for ii=0,nell-1 do begin

;;     print, min(u), max(u)

     if ii eq 0 then begin  ;; get full map case upto lmax

        lmax = ell[ii]

        print, 'computing theory MFs for full data upto lmax,mpow=',lmax,mpow

        ell2=2.*(2.+1.)/2.
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

        print, 'computing theory MFs for Multipole, mpow=',ell[ii],mpow
        ell2=double(ell[ii])*(double(ell[ii])+1.)/2.

        ;print, 'lk_cu mpow=',mpow
        lk_mfs, u, ell2, lk0, lk1, lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy,lkxz=lkxz

        if keyword_set(lkmask) then begin
        ;print, 'lk_cu lkmask=',lkmask

           rho2=lkxy/(4.*!pi)
           rho1=4*lkxz/(2*!pi)
           rho0 = lkxx/2.
           
           M0 = lkmask[1]
           M1 = lkmask[2]
           M2 = lkmask[0]

           area_vec = rho0*M2
           if keyword_set(fsky) then area_vec = lk2 ;;rho0*M2
           length_vec = rho0*M1 + (!pi/2.)*sqrt(lkl)*rho1*M2
           genus_vec = rho0*M0 + sqrt(lkl)*rho1*M1 + lkl*rho2*M2
           
           lk0 = genus_vec
           lk1 = length_vec
           lk2 = area_vec
        endif

        lk_variance,phiu=u,ell=ell[ii],var=gvar,mpow=mpow,nvec=num

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
     ;; mftry[*,1,ii]=area[*,ii]
     ;; mftry[*,2,ii]=genus[*,ii]
     ;; mftry[*,3,ii]=length[*,ii]


     varmftry[*,2,ii]=reform(gvar)
  endfor

end
