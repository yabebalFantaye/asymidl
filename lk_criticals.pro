function lk_crit_integrand,u
     
  return, sqrt(3d0/(8*!pi))*(2*exp(-u^2)+u^2-1)*exp(-u^2/2d0)
  
end

function lk_extrema_integrand,u
     
  return, sqrt(3d0/(8*!pi))*(exp(-u^2)+u^2-1)*exp(-u^2/2d0) 
  
end

function lk_saddle_integrand,u
     
  return, sqrt(3d0/(8*!pi))*exp(-3*u^2/2d0)
  
end


function lk_criticals_int,umin, umax, name=name
  ;;
  ;;Get the PDF of expected value of
  ;;criticals, extrema, or saddle by
  ;;band integration
  ;;
  if not keyword_set(name) then name='critical'
  
  if name eq 'critical' then $
     func='lk_crit_integrand'
     
  if name eq 'extrema' then $
     func='lk_extrema_integrand'
     
  if name eq 'saddle' then $
     func='lk_saddle_integrand'

  return, QROMB(func, umin, umax,jmax=10,eps=1e-6)

end

function lk_criticals_center, umin, umax, name=name
  ;;
  ;;Get the PDF value of the expected 
  ;;criticals, extrema, or saddle by
  ;;evauating at u
  ;;
  if not keyword_set(name) then name='critical'

  u_cent = umin + (umax-umin)/2d0
  
  if name eq 'critical' then $
     res=lk_crit_integrand(u_cent)     
     
  if name eq 'extrema' then $
     res=lk_extrema_integrand(u_cent)          
     
  if name eq 'saddle' then $
     res=lk_saddle_integrand(u_cent)             


  return, res

end

function lk_criticals,u,name=name,prec=prec
  ;;computes the expectation of criticals,
  ;;extremas and saddles
  ;;
  ;; Inputs: excursion level u and type of extrema
  
  if not keyword_set(name) then name='critical'

  coef=2d0/sqrt(3d0)
  nsteps=n_elements(u)
  res=dblarr(nsteps)

  delta_u=(u[1]-u[0])  
  umin=u[0]-delta_u
  
  for i=0,nsteps-1 do begin
     if i ne 0 then umin=u[i-1]
     
     if delta_u gt 1e-1 or keyword_set(prec) then  $
        res[i]=lk_criticals_int(umin,u[i],name=name) $
     else $
        res[i]=lk_criticals_center(umin,u[i],name=name)        
  endfor
  
  return,coef*res 
end
