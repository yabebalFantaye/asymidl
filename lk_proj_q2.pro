function lk_proj_q2,u,lvec,clhat
;; Compute Wigner chaos projecte Lipschitz-Killing curvatures


  ;;constants
  half=0.5
  over2pi=1d0/(2d0*!pi)

   ;;important terms
  bigPhi=gaussint(u)              ;;Gaussian CDF
  phiu=exp(-u^2/2d0)/sqrt(2*!pi)    ;;Gaussian density

  ;;Hermite polynomial n=-1
  H_n1u= sqrt(2.*!pi)*(1d0 - bigPhi)*exp(u^2/2.) 

  ;;Hermite polynomials n>=0
  H_0u = 1d0
  H_1u = u
  H_2u = u^2 - 1d0

  genus_upart=H_2u*H_1u*phiu

  ;; print, 'u',u[300:350]
  ;; print, 'u-dependent part: ',genus_upart[300:350]
  ;; print, 'phiu',phiu[300:350]

  nell=n_elements(lvec)
  nu=n_elements(u)
  proj=dblarr(nu,3,nell)

  for ii=0,nell-1 do begin
     ell=double(lvec[ii])
     cell=clhat[ii]

     ;;important terms  
     lambda_l=ell*(ell+1d0)/2d0 ;;all terms are lambda_l/2
     

     ;;Integral( H_2(f_l(x))dx )
     intH2flx=(cell-1d0) ;;*4d0*!pi 

     ;print, 'lk_proj_12: ell, cell, intH2, lambda_l',ell, cell,intH2flx,sqrt(lambda_l)
     ;;print, 'u[250:260]',u[250:260]

     ;;k=0; realised - expected Euler-Poincare characteristic: proj(L0)
     proj[*,0,ii] = lambda_l*genus_upart*intH2flx ;;/(2*!pi)
  
     ;;k=1; realised - expected (half) boundary length: proj(L1)
     proj[*,1,ii] = sqrt(lambda_l*!pi/8d0)*H_1u^2*phiu*intH2flx*(4*!pi) ;*half
     
     ;;k=2; realised - expected Area: proj(L2)  
     proj[*,2,ii] = half*H_0u*H_1u*phiu*intH2flx

  endfor
  
  return, proj

end
