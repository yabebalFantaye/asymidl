function lk_varcrit_integrand,t
     
  return, exp(-3*t^2d0/2d0)*(2d0-6d0*t^2d0-exp(t^2d0)*(1d0-4d0*t^2d0+t^4d0))
  
end

function lk_varextrema_integrand,t

  return, exp(-3*t^2d0/2d0)*(1d0-3d0*t^2d0-exp(t^2d0)*(1d0-4d0*t^2d0+t^4d0))  
  
end

function lk_varsaddle_integrand,t
     
  return, exp(-3d0*t^2d0/2d0)*(1d0-3d0*t^2d0)
  
end

function lk_var_criticals_cdf, u, name=name

  if not keyword_set(name) then name='critical'

  ;;The limit of the intergration is very important
  ;;The shape of the variance entirely depends on that
  maxu=max(u) 

  if maxu gt 3.5 then begin
     res=1d0/(8*!pi)*exp(-3*u^2)*u^2*(2+exp(u^2)*(u^2-1))^2
     
     if name eq 'extrema' then $
        res=1d0/(8*!pi)*exp(-3*u^2)*u^2*(1+exp(u^2)*(u^2-1))^2     
     
     if name eq 'saddle' then $
        res=1d0/(8*!pi)*exp(-3*u^2)*u^2
  endif else begin
     
     coef=1d0/(8*!pi)
     nsteps=n_elements(u)
     res=dblarr(nsteps)
          
     for i=0,nsteps-1 do begin
        ;;print, 'qromb integration interval: ',name,u[i],maxu
        
        if name eq 'critical' then $
           res[i]=QROMB('lk_varcrit_integrand', u[i],maxu,jmax=10,eps=1e-6)
        
        if name eq 'extrema' then $
           res[i]=QROMB('lk_varextrema_integrand', u[i],maxu,jmax=10,eps=1e-6)
        
        if name eq 'saddle' then $
           res[i]=QROMB('lk_varsaddle_integrand', u[i],maxu,jmax=10,eps=1e-6)
     endfor
     res=coef*res^2
  endelse
  
  return, res

end  

;;-------------------------------

function lk_var_criticals,u,name=name,pdf=pdf
  ;;computes the variance of criticals,
  ;;extremas and saddles
  ;;
  ;; Inputs: excursion level u and type of extrema
  
  if not keyword_set(name) then name='critical'

  if not keyword_set(pdf) then begin

     res=lk_var_criticals_cdf(u, name=name)
     
  endif else begin
     
     coef=1d0/(8*!pi)
     nsteps=n_elements(u)
     res=dblarr(nsteps)
     minu0=min(u)-(u[1]-u[0])
     
     for i=0,nsteps-1 do begin
        if i eq 0 then minu=minu0 else minu=u[i-1]        
        ;;print, 'qromb integration interval: ',minu,u[i]
        
        if name eq 'critical' then $
           res[i]=QROMB('lk_varcrit_integrand', minu, u[i],jmax=10,eps=1e-6)
     
        if name eq 'extrema' then $
           res[i]=QROMB('lk_varextrema_integrand', minu, u[i],jmax=10,eps=1e-6)
     
        if name eq 'saddle' then $
           res[i]=QROMB('lk_varsaddle_integrand', minu, u[i],jmax=10,eps=1e-6)
     endfor
     res=coef*res^2
     
  endelse
  
  return,res
end
