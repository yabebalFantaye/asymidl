function get_cdf_mms,minmaxsad
  ;;
  ;;Return MinMaxSad as cumulative distribution [u,inf]
  ;;
  
  nsteps=n_elements(minmaxsad[*,0,0,0])  
  mms=minmaxsad
  ;;cumulative distribution [u, inif]
  for i=0,nsteps-1 do begin
     mms[i,*,*,*]=total(minmaxsad[i:nsteps-1,*,*,*], 1) ;/total(minmaxsad[*,*,*,*], 1)
  endfor
  return, mms

end

function get_pdf_mms,minmaxsad
  ;;
  ;;Return MinMaxSad as cumulative distribution [u,inf]
  ;;
  
  nsteps=n_elements(minmaxsad[*,0,0,0])  
  mms=minmaxsad
  ;;cumulative distribution [u, inif]
  for i=0,nsteps-1 do begin
     mms[i,*,*,*]=reform(minmaxsad[i,*,*,*])/total(minmaxsad, 1)
  endfor
  return, mms

end


function get_mms,minmaxsad,t=t,j=j
  ;;
  ;;Return Minima, Maxima or Saddle (mms)
  ;;
  
  if not keyword_set(t) then t=2 ;;default saddle
  
  if keyword_set(j) then begin
     mms = reform(minmaxsad[*,t,j,*])
  endif else begin
     mms = reform(minmaxsad[*,t,*,*])
  endelse
  
  return, mms
end

function get_mean_mms,minmaxsad,t=t,j=j
  ;;
  ;;Return the mean across simulations of
  ;;Minima, Maxima, or Saddle
  ;;
  
  if not keyword_set(t) then t=2 ;;default saddle
  
  if keyword_set(j) then begin
     mean_mms = mean(reform(minmaxsad[*,t,j,*]),dimension=2)
  endif else begin
     mean_mms = mean(reform(minmaxsad[*,t,*,*]),dimension=3)
  endelse
  
  return, mean_mms
end

function get_var_mms,mms,t=t,j=j
  ;;
  ;;Return the variance across simulation of
  ;;Minima, Maxima, or Saddle
  ;;
  
  if not keyword_set(t) then t=2 ;;default saddle
  
  if keyword_set(j) then begin
     var_mms = variance(reform(mms[*,t,j,*]),dimension=2)
  endif else begin
     var_mms = variance(reform(mms[*,t,*,*]),dimension=3)
  endelse
  
  return, var_mms
end


;;--------------------------------


function get_crit,minmaxsad,j=j
  
  if keyword_set(j) then begin
     crit = reform(minmaxsad[*,0,j,*]+minmaxsad[*,1,j,*]+minmaxsad[*,2,j,*])
  endif else begin
     crit = reform(minmaxsad[*,0,*,*]+minmaxsad[*,1,*,*]+minmaxsad[*,2,*,*])
  endelse
     
  return, crit
end

function get_mean_crit,minmaxsad,j=j
  mcrit = mean(get_crit(minmaxsad,j=j),dimension=2)
  return, mcrit
end

function get_var_crit,mms,j=j
  if not keyword_set(j) then j=1
  
  dcrit = variance(get_crit(mms,j=j),dimension=2)
  return, dcrit
end

;;--------------------------------

function get_extrema,minmaxsad,j=j
  if keyword_set(j) then begin  
     extrema = reform(minmaxsad[*,0,j,*] + minmaxsad[*,1,j,*])
  endif else begin
     extrema = reform(minmaxsad[*,0,*,*] + minmaxsad[*,1,*,*])     
  endelse
  
  return, extrema
end

function get_mean_extrema,minmaxsad,j=j
  if not keyword_set(j) then j=1
  mextrema = mean(get_extrema(minmaxsad,j=j),dimension=2)
  return, mextrema
end

function get_var_extrema,mms,j=j
  if not keyword_set(j) then j=1

  dextrema = variance(get_extrema(mms,j=j),dimension=2)
  return, dextrema
end


;;--------------------------------
;;   Wrapper functions to all the above
;;--------------------------------
function get_peak,minmaxsad,j=j,name=name
  ;;
  ;;returns PDF of local peaks
  ;;
  
  if not keyword_set(name) then name='critical'
  case name of
     'extrema':res=get_extrema(minmaxsad,j=j) 
     'saddle': res=get_mms(minmaxsad,t=2,j=j)
     'maxima': res=get_mms(minmaxsad,t=1,j=j)
     'minima': res=get_mms(minmaxsad,t=0,j=j)
     else:     res=get_crit(minmaxsad,j=j)
  endcase
  return, res
end

function get_mean_peak,minmaxsad,j=j,name=name
  ;;
  ;;returns mean PDF of local peaks
  ;;
  
  if not keyword_set(name) then name='critical'
  if not keyword_set(j) then j=1
  
  case name of
     'extrema':res=get_mean_extrema(minmaxsad,j=j) 
     'saddle': res=get_mean_mms(minmaxsad,t=2,j=j)
     'maxima': res=get_mean_mms(minmaxsad,t=1,j=j)
     'minima': res=get_mean_mms(minmaxsad,t=0,j=j)
     else:     res=get_mean_crit(minmaxsad,j=j)
  endcase
  return, res
end

function get_var_peak,minmaxsad,j=j,name=name,norm=norm
  ;;
  ;;returns variance of CDF of local peaks
  ;;
  
  if not keyword_set(name) then name='critical'
  if not keyword_set(j) then j=1
  mms=minmaxsad
  if keyword_set(norm) then mms=get_cdf_mms(minmaxsad) 
  
  case name of
     'extrema':res=get_var_extrema(mms,j=j) 
     'saddle': res=get_var_mms(mms,t=2,j=j)
     'maxima': res=get_var_mms(mms,t=1,j=j)
     'minima': res=get_var_mms(mms,t=0,j=j)
     else:     res=get_var_crit(mms,j=j)
  endcase
  return, res
end


;;--------------------------------
;;  Main start here
;;--------------------------------

pro load_peaks, fres,minmaxsad, nmms=nmms, peakmap=peakmap,$
                nsim_start=nsim_start,ipol=ipol

  ;;compile all functions with simple init
  ;;
  resolve_routine, 'get_peak', /is_function  
  resolve_routine, 'get_mean_peak', /is_function
  resolve_routine, 'get_var_peak', /is_function
  
  
  nsimb=1
  if keyword_set(nsim_start) then nsimb=nsim_start

  fend=''
  if keyword_set(ipol) then begin
     if ipol gt 1 then fend='_ipol'+strn(ipol)
  endif
  
  nsteps=n_elements(minmaxsad[*,0,0,0])
  nj=n_elements(minmaxsad[0,0,*,0])
  nsim=n_elements(minmaxsad[0,0,0,*])

  ;;define working arrays 
  nmms=lindgen(1,3,nj,nsim)
  res0_count=lindgen(nsteps+1,3,nj)
  
  for i=0,nsim-1 do begin
     isim=nsimb+i-1

     simno= strn(isim)          ;trim(string(i),2)
     
     if not keyword_set(peakmap) then begin
        
        file=fres+'.unf_minmaxsad_insim'+simno+fend
        runf,res0_count,file
        
     endif else begin
        
        for ij=0,nj-1 do begin           
                                ;file=fres+'.unf_ipix_val_map'+simno+fend+'_j'+str(ij)           
           
           ;; if phiu is given, vmaxima etc. are pdf of maxima etc.
           peaks_pix=read_mf_peaks(fres,i,ij,iminima=iminima,vminima=vminima,phiu=phiu,$
                                   imaxima=imaxima,vmaxima=vmaxima,rms=rms,avg=avg,$                             
                                   isaddle=isaddle,vsaddle=vsaddle,minmax=minmaxb, $
                                   /nopmap,/nested, /inring) ;; ,pixrad=pixrad)
           
           res0_count[1:nsteps,0,ij]=vminima
           res0_count[1:nsteps,1,ij]=vmaxima
           res0_count[1:nsteps,2,ij]=vsaddle                 
        endfor              
     endelse         
     
     ;;temp=res0_count[0,0,*,*]
     ;;res0_count[0,0,*,*]=temp-res0_count[1,0,*,*]+res0_count[2,0,*,*]
     res0_count[1,0,*,*]=res0_count[2,0,*,*]           
     nmms[*,*,*,i]=res0_count[0,*,*,*]
     minmaxsad[*,*,*,i]=res0_count[1:nsteps,*,*,*]
  endfor
  


end  
