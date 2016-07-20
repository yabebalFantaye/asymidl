pro lk_variance,phiu=phiu,ell=ell,var=var,_extra=extra
;;compute variance of the Euler-Poincare characteristics

  if not keyword_set(ell) then ell=1
  if not keyword_set(mpow) then mpow=1
  if not keyword_set(phiu) then phiu=2d0*3*(findgen(100)/(100-1)-0.5) ;;[-3,3]

  u=phiu


  vx=(u^2)*exp(-u^2)/(4*!pi)
  vy=(u-u^3)^2*exp(-u^2)/(8*!pi)
  var = double(ell)^3*vy



end
