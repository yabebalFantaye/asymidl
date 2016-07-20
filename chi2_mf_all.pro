pro chi2_mf_all,mfpol, chi2mat,out_struct=out_struct, nrep=nrep, iout=iout


  if n_params() lt 1 then begin
     print, 'usage: chi2_mf,mcres, out_struct, nrep=nrep, iout=iout'
     return
  endif

  if not keyword_set(nrep) then nrep=1
  if not keyword_set(iout) then iout=0

  nsteps=n_elements(mfpol[0].mcres[*,0,0,0])
  nmf=n_elements(mfpol[0].mcres[0,*,0,0])
  nj = n_elements(mfpol[0].mcres[0,0,*,0])
  nsim = n_elements(mfpol[0].mcres[0,0,0,*])

  nmethod=n_elements(mfpol)

  print, 'chi2_mf: nsim, nj,nsteps,nrep ',nsim, nj, nsteps,nrep

  ;; mres=dblarr(nsteps,nmf,nj)
  ;; dres=dblarr(nsteps,nmf,nj)
  ;; cov_mat=dblarr(nsteps,nsteps,nmf,nj)
  ;; xcor_mat=dblarr(nsteps,nsteps,nmf,nj)

  chi2mat=dblarr(nsim,nj,nmethod)
  
  newsteps=n_elements(mfpol[0].mcres[0:nsteps-1,1,0,0])

  ;;compute covmat of area, genus and length
  for imd=0,nmethod-1 do begin
     for k=0,nj-1 do begin          
     
        simres=dblarr(3*newsteps,nsim)
        for i=0,nsim-1 do begin 
           zzz = [mfpol[imd].mcres[0:nsteps-1,1,k,i],mfpol[imd].mcres[0:nsteps-1,2,k,i],mfpol[imd].mcres[0:nsteps-1,3,k,i]]
           simres[*,i]=reform(zzz)
        endfor

        simres=transpose(simres)

        mc_reduce, simres, mres_ik,dres_ik, covmat_ik, xcormat_ik

        cinv=invert(covmat_ik) 

        for isim=0,nsim-1 do begin
           diff_sim = reform(simres[isim,*] - reform(mres_ik))  
           chi2mat[isim,k,imd] = transpose(diff_sim) # (cinv # diff_sim)
        endfor
     
     endfor
  endfor

help, chi2mat

;;=================================================
;;=============== structure out put ================

  if (size(out_struct,/type) eq 8) and (nsim gt 1) then begin ;;structure is given as input
     print, 'allchi2 saving load_res parameters on already DEFINED structure: nsim>10. nrep, iout=', nrep, iout
     out_struct[iout].chi2mat=chi2mat
  endif


  if size(out_struct,/type) ne 8 and (nsim gt 1) then begin ;;structure is given as input
     print, 'allchi2 saving load_res DEFINING NEW structure: nsim>10. nrep=',nrep
     out_struct = {chi2mat:chi2mat}
     out_struct = replicate(out_struct,nrep)
  endif


end
