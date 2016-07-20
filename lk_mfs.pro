pro lk_mfs,u,ell2,lk0,lk1,lk2,mpow=mpow,lkxx=lkxx,lkxy=lkxy,lkyy=lkyy,lkxz=lkxz

  if not keyword_set(mpow) then mpow=1

     if mpow eq 1 then begin

        phi=gaussint(u)

        lkxx= 2.*(1.-phi)
        lkxy=u*exp(-u^2/2.)*4*!pi/sqrt((2*!pi)^3)
        lkxz=exp(-u^2/2.)/4.
        lk0 = lkxx + ell2*lkxy

        lkyy = !pi*exp(-u^2/2.)
        lk1 = sqrt(ell2)*lkyy


        lk2 = 4*!pi*(1.-phi)   
        
     endif 

     if mpow eq 2 then begin

        newu=u ;+1
        sqrtu1= sqrt(newu)
        phi=gaussint(sqrtu1)

        lkxx = 4.*(1.-phi)
        lkxy = 4*sqrtu1*exp(-newu/2.)/sqrt(2*!pi)
        lkxz = exp(-newu/2.)
        lk0 = lkxx + ell2*lkxy

;print, 'phi, squ,xx, xy, genus: ',sqrtu1[100],phi[100],lkxx[100],lkxy[100],lk0[100]

        lkyy=2.*!pi*exp(-newu/2.)
        lk1 = sqrt(ell2)*lkyy

        lk2 = 4*!pi*2.*(1.-phi)

     endif 


     if mpow eq 3 then begin

        newu=(abs(u)/u)*(abs(u))^(1d0/3d0) ;^(3d0/2d0)
        ;newu=u^(3d0/2d0)
        phi=gaussint(newu)

        lkxx=2.*(1.-phi)
        lkxy=2*newu*exp(-newu^2./2.)/sqrt(2*!pi)
        lkxz=2*exp(-newu^2./2.)/4.
        lk0 = lkxx + ell2*lkxy


        lkyy=!pi*exp(-newu^2./2.)
        lk1 = sqrt(ell2)*lkyy

        lk2 = 4*!pi*(1.-phi)

     endif 


end
