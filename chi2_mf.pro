pro chi2_mf,mcres, out_struct=out_struct, nrep=nrep, iout=iout


  if n_params() lt 1 then begin
     print, 'usage: chi2_mf,mcres, out_struct, nrep=nrep, iout=iout'
     return
  endif

  if not keyword_set(nrep) then nrep=1
  if not keyword_set(iout) then iout=0

  nsteps=n_elements(mcres[*,0,0,0])
  nmf=n_elements(mcres[0,*,0,0])
  nj = n_elements(mcres[0,0,*,0])
  nsim = n_elements(mcres[0,0,0,*])

  print, 'chi2_mf: nsim, nj,nsteps,nrep ',nsim, nj, nsteps,nrep


  mres=dblarr(nsteps,nmf,nj)
  dres=dblarr(nsteps,nmf,nj)

  cov_mat=dblarr(nsteps,nsteps,nmf,nj)
  xcor_mat=dblarr(nsteps,nsteps,nmf,nj)

  chi2mat=dblarr(nsim,nmf,nj)
  

  ;;compute covmat of area, genus and length
  for k=0,nj-1 do begin                         
     for i=1,nmf-1 do begin 
        simres = transpose(reform(mcres[*,i,k,1:nsim-1]))

        mc_reduce, simres, mres_ik,dres_ik, covmat_ik, xcormat_ik
     
        ;;help, mres_ik
        ;;help, covmat_ik

        mres[*,i,k]=reform(mres_ik)
        dres[*,i,k]=reform(dres_ik)
        cov_mat[*,*,i,k]=covmat_ik 
        xcor_mat[*,*,i,k]=xcormat_ik

        cinv=invert(covmat_ik) 

        for isim=0,nsim-1 do begin
           diff_sim = reform(mcres[*,i,k,isim]-mres[*,i,k])  

           ;; help, diff_sim # covmat_ik
           ;; help,  (diff_sim # covmat_ik) # diff_sim 

           chi2mat[isim,i,k] = transpose(diff_sim) # (cinv # diff_sim)
        endfor

     endfor
  endfor

help, chi2mat

;;=================================================
;;=============== structure out put ================

  if (size(out_struct,/type) eq 8) and (nsim gt 1) then begin ;;structure is given as input
     print, 'saving load_res parameters on already DEFINED structure: nsim>10. nrep, iout=', nrep, iout
     out_struct[iout].mres=mres
     out_struct[iout].dres=dres
     out_struct[iout].cov=cov_mat
     out_struct[iout].xcor=xcor_mat
     out_struct[iout].chi2=chi2mat
  endif


  if size(out_struct,/type) ne 8 and (nsim gt 1) then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim>10. nrep=',nrep
     out_struct = {mres:mres,dres:dres,cov:cov_mat,xcor:xcor_mat,chi2:chi2mat}
     out_struct = replicate(out_struct,nrep)
  endif


end
