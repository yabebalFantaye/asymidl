pro load_res_mf,fres,mcres,minmaxsad,out_struct=out_struct,dobin=dobin,nzero=nzero, nrep=nrep, iout=iout,$
                ires=ires,nsplit=nsplit,overwrite=overwrite,mres=mres,dres=dres,mcount=mcount,dcount=dcount,nsimb=nsimb,nsimf=nsimf,$
                nproc=nproc, nodir=nodir,fkey=fkey,restry=restry,mmstry=mmstry,feedback=feedback,ipol=ipol,mpow=mpow,unnorm=unnorm,$
                nopeaks=nopeaks, gnorm=gnorm, lnorm=lnorm,mfnorm=mfnorm

;;read spectra from cmb sim:

  if not keyword_set(ires) then ires=0
  if not keyword_set(nzero) then nzero=4
  if not keyword_set(nrep) then  nrep=1
  if not keyword_set(iout) then  iout=0
  if not keyword_set(fkey) then fkey=''
  if not keyword_set(feedback) then feedback=0
  if not keyword_set(mpow) then mpow=1

  if not keyword_set(gnorm) then gnorm=1d0
  if not keyword_set(lnorm) then lnorm=1d0


;print,'out_struct defined or not: type=',size(out_struct,/type)
;help, out_struct

  infores = size(mcres)
  print, 'infro mcres', infores

  nsteps=n_elements(mcres[*,0,0,0])
  nmf=n_elements(mcres[0,*,0,0])
  nj = n_elements(mcres[0,0,*,0])
  nsim = n_elements(mcres[0,0,0,*])
  nsim_orig=nsim

  print, 'nsim, nj,nsteps ',nsim, nj, nsteps


  if not keyword_set(nsimb) then nsimb=1
  if not keyword_set(nsimf) then nsimf=nsim
  nsim=nsimf-nsimb+1

  res00=dblarr(nsteps,nmf,nj)
  if nmf gt 4 then nmf=4

  res0=dblarr(nsteps,nmf,nj)
  mcres = dblarr(nsteps,nmf,nj,nsim)



  if not keyword_set(nopeaks) then begin    
     minmaxsad=dblarr(nsteps,3,nj,nsim)
     tminmaxsad=lindgen(1,3,nj,nsim)
     norm_count=dblarr(nsteps,3,nj)
     res0_count=lindgen(nsteps+1,3,nj)
  endif


  fend=''
  if keyword_set(ipol) then begin
     if ipol gt 1 then fend='_ipol'+strn(ipol)
  endif


  count=0
  for i=0,nsim-1l do begin
     isim=nsimb+i-1

     simno= strn(isim)          ;trim(string(i),2)
                                ;simno=str_replicate('0',nzero-strlen(simno))+simno
     file=fres+'.unf_insim'+simno+fend
                                ;print, file
     

     if file_test(file,/regular) then begin

        if feedback gt 1 then print, 'reading: ',file
        runf,res00,file          
        res0[*,0:nmf-1,*] = res00[*,0:nmf-1,*]

        ;;normalized MF
        if keyword_set(mfnorm) then begin
           sig = reform(res0[0,0,*]) ;;rms of map
           dersig = reform(res0[1,0,*]) ;;rms of der(map)
           lk_data_norm, sig, dersig, Ak=mfnorm ;;obtain A_k

           for ij=0,nj-1 do begin
              res0[*,1,ij]=mfnorm[ij,0]*res0[*,1,ij] ;;area
              res0[*,3,ij]=mfnorm[ij,1]*res0[*,3,ij] ;;length
              res0[*,2,ij]=mfnorm[ij,2]*res0[*,2,ij] ;;genus
           endfor
        endif

        ;; if not keyword_set(unnorm) then begin
        ;;    if mpow eq 2 then begin
        ;;       for jj=1,nmf-1 do begin                               
        ;;          for kk=0,nj-1 do begin              
        ;;             res0[(nsteps/2+10):nsteps-1,jj,kk]=double(res0[(nsteps/2+10):nsteps-1,jj,kk])/max(res0[(nsteps/2+10):nsteps-1,jj,kk]*1d0)
        ;;          endfor
        ;;       endfor
        ;;    endif else begin
        ;;       for jj=1,nmf-1 do begin                               
        ;;          for kk=0,nj-1 do begin              
        ;;             res0[*,jj,kk]=double(res0[*,jj,kk])/max(res0[*,jj,kk]*1d0)
        ;;          endfor
        ;;       endfor
        ;;    endelse
        ;; endif

        mcres[*,*,*,i]=res0

        if not keyword_set(nopeaks) then begin        
           file=fres+'.unf_minmaxsad_insim'+simno+fend
           if feedback gt 1 then print, 'minmax file: ',file
           runf,res0_count,file
           
           norm_count=double(res0_count[1:nsteps,*,*,*])
           
           tminmaxsad[*,*,*,i]=res0_count[0,*,*,*]
           minmaxsad[*,*,*,i]=norm_count
        endif

        count = count+1
        
     endif else print, 'skipping not found file: '+file
  endfor

  nsim=count

  print, 'load_res_mf: nsim_orig, nsim_read',nsim_orig,nsim

  mcres = mcres[*,*,*,0:nsim-1]


  mcres[*,2,*,*]=mcres[*,2,*,*]*gnorm        ;;length and genus normalized by mfnorm
  mcres[*,3,*,*]=mcres[*,3,*,*]*lnorm


  mcres_diff=mcres
  mcres_ratio=mcres
  res =reform(mcres[*,*,*,ires])
  help, res

  if not keyword_set(nopeaks) then begin        
     minmaxsad=minmaxsad[*,*,*,0:nsim-1]
     minmaxsad_diff = minmaxsad
     minmaxsad_ratio = minmaxsad
     res_count = reform(minmaxsad[*,*,*,ires])
     help, res_count
  endif




  
;;;;;this is set by hand here because the simulation number 2 for
;;;;;planck ffp6 is set to be noise only (no CMB)
;;if keyword_set(repires) then mcres[*,*,*,repires]=mcres[*,*,*,3]


  if nsim gt 1 then begin

     
     decompose,fres,d,path,name,ext,ver
     ext='.unf'
     dirout=path+'/post_nj'+strn(nj)+'_'+fkey+'/'

     print,'dirout = '+dirout
     spawn, 'mkdir -p '+dirout
     
     fmres = dirout+'mres_'+name+ext
     fdres = dirout+'dres_'+name+ext

     fmcres_diff = dirout+'mcres_diff_'+name+ext
     fmcres_ratio = dirout+'mcres_ratio_'+name+ext
     fdres_diff = dirout+'dres_diff_'+name+ext
     fdres_ratio = dirout+'dres_ratio_'+name+ext

     mres=dblarr(nsteps,nmf,nj)
     dres=mres
     dres_diff=mres
     dres_ratio=mres


     if not keyword_set(nopeaks) then begin        
        fminmaxsad_diff = dirout+'minmaxsad_diff_'+name+ext
        fminmaxsad_ratio = dirout+'minmaxsad_ratio_'+name+ext
        fdcount_diff = dirout+'dcount_diff_'+name+ext
        fdcount_ratio = dirout+'dcount_ratio_'+name+ext

        fmcount = dirout+'mcount_'+name+ext
        fdcount = dirout+'dcount_'+name+ext

        mtcount=dblarr(nj)
        dtcount=mtcount
        ttt = (total(minmaxsad[*,0,*,*]+minmaxsad[*,1,*,*]+minmaxsad[*,2,*,*],1))
        for k=0,nj-1 do begin              
           mtcount[k]=mean(reform(ttt[0,k,*]))
           dtcount[k]=stdev(reform(ttt[0,k,*]))
        endfor
        
        mcrit=dblarr(nsteps,nj)
        dcrit=mcrit
        
        mcount=dblarr(nsteps,3,nj)
        dcount=mcount
        dcount_diff=mcount
        dcount_ratio=mcount

     endif



     if  not file_test(fmres,/regular) or keyword_set(overwrite) then begin
        
        for j=0,nsteps-1 do begin
           
           if j mod 100 eq 0 then print,'mean var nsteps = ',j
           
           for k=0,nj-1 do begin              

              for i=0,nmf-1 do begin                               
                 mres[j,i,k]=mean(mcres[j,i,k,*])
                 dres[j,i,k]=stdev(mcres[j,i,k,*])

                 if keyword_set(restry) then begin
                    mcres_diff[j,i,k,*] = mcres[j,i,k,*]-restry[j,i,k]
                    dres_diff[j,i,k]=stdev(mcres_diff[j,i,k,*])
                    
                    mcres_ratio[j,i,k,*] = (mcres[j,i,k,*]-restry[j,i,k])/(restry[j,i,k]+1e-12)
                    dres_ratio[j,i,k]=stdev(mcres_ratio[j,i,k,*])
                 endif else begin
                    mcres_diff[j,i,k,*] = mcres[j,i,k,*]-mres[j,i,k]
                    dres_diff[j,i,k]=stdev(mcres_diff[j,i,k,*])
                    
                    mcres_ratio[j,i,k,*] = (mcres[j,i,k,*]-mres[j,i,k])/(mres[j,i,k]+1e-12)
                    dres_ratio[j,i,k]=stdev(mcres_ratio[j,i,k,*])
                 endelse
              endfor


              ;;peaks stat
              if not keyword_set(nopeaks) then begin        

                 ttt=reform(total(minmaxsad[j:nsteps-1,0,k,*]+minmaxsad[j:nsteps-1,1,k,*]+minmaxsad[j:nsteps-1,2,k,*],1))
                 mcrit[j,k]=mean(ttt)
                 dcrit[j,k]=stdev(ttt)
                 
                 for i=0,2 do begin                               
                    mcount[j,i,k]=mean(minmaxsad[j,i,k,*])
                    dcount[j,i,k]=stdev(minmaxsad[j,i,k,*])
                    
                    
                    if keyword_set(mmstry) then begin
                       minmaxsad_diff[j,i,k,*] = minmaxsad[j,i,k,*]-mmstry[j,i,k]
                       dcount_diff[j,i,k]=stdev(minmaxsad_diff[j,i,k,*])
                       
                       minmaxsad_ratio[j,i,k,*] = (minmaxsad[j,i,k,*]-mmstry[j,i,k])/(mmstry[j,i,k]+1e-12)
                       dcount_ratio[j,i,k]=stdev(minmaxsad_ratio[j,i,k,*])
                    endif else begin 
                       minmaxsad_diff[j,i,k,*] = minmaxsad[j,i,k,*]-mcount[j,i,k]
                       dcount_diff[j,i,k]=stdev(minmaxsad_diff[j,i,k,*])
                       
                       minmaxsad_ratio[j,i,k,*] = (minmaxsad[j,i,k,*]-mcount[j,i,k])/(mcount[j,i,k]+1e-12)
                       dcount_ratio[j,i,k]=stdev(minmaxsad_ratio[j,i,k,*])
                    endelse
                    
                 endfor

              endif ;;for no peaks


           endfor
        endfor
        
        ;; wunf, mres, fmres
        ;; wunf, dres, fdres

        ;; wunf, mcount, fmcount
        ;; wunf, dcount, fdcount


        ;; wunf, mcres_diff, fmcres_diff
        ;; wunf, mcres_ratio, fmcres_ratio
        ;; wunf, dres_diff, fdres_diff
        ;; wunf, dres_ratio, fdres_ratio

        ;; wunf, minmaxsad_diff, fminmaxsad_diff
        ;; wunf, minmaxsad_ratio, fminmaxsad_ratio
        ;; wunf, dcount_diff, fdcount_diff
        ;; wunf, dcount_ratio, fdcount_ratio
     endif else begin
        print, 'loading  ',fmres
        runf, mres, fmres
        print, 'loading  ',fdres
        runf, dres, fdres

        runf, mcres_diff, fmcres_diff
        runf, mcres_ratio, fmcres_ratio
        runf, dres_diff, fdres_diff
        runf, dres_ratio, fdres_ratio


        if not keyword_set(nopeaks) then begin        
           print, 'loading  ',fmcount
           runf, mcount, fmcount
           print, 'loading  ',fdcount
           runf, dcount, fdcount
           
           runf, minmaxsad_diff, fminmaxsad_diff
           runf, minmaxsad_ratio, fminmaxsad_ratio
           runf, dcount_diff, fdcount_diff
           runf, dcount_ratio, fdcount_ratio
        endif

     endelse

     print, 'sanity mres, dres: ',total(mres),total(dres)

  endif

  res_diff =reform(mcres_diff[*,*,*,ires])
  res_ratio =reform(mcres_ratio[*,*,*,ires])

  if not keyword_set(nopeaks) then begin        
     res_count_diff = reform(minmaxsad_diff[*,*,*,ires])
     res_count_ratio = reform(minmaxsad_ratio[*,*,*,ires])
  endif

  if (size(out_struct,/type) eq 8) and (nsim gt 1) then begin ;;structure is given as input
     print, 'saving load_res parameters on already DEFINED structure: nsim>10'
     out_struct[iout].res=res
     out_struct[iout].mcres=mcres
     out_struct[iout].mres=mres
     out_struct[iout].dres=dres

     out_struct[iout].mcres_diff=mcres_diff
     out_struct[iout].dres_ratio=dres_ratio
     out_struct[iout].dres_diff=dres_diff
     out_struct[iout].mcres_ratio=mcres_ratio


     if not keyword_set(nopeaks) then begin        
        ;;out_struct[iout].tminmaxsad=tminmaxsad
        out_struct[iout].res_count=res_count
        out_struct[iout].mcount=mcount
        out_struct[iout].dcount=dcount
        
        out_struct[iout].mtcount=mtcount
        out_struct[iout].dtcount=dtcount
        
        out_struct[iout].mcrit=mcrit
        out_struct[iout].dcrit=dcrit

        out_struct[iout].minmaxsad_diff=minmaxsad_diff
        out_struct[iout].minmaxsad_ratio=minmaxsad_ratio

        out_struct[iout].res_count_diff=res_count_diff
        out_struct[iout].res_count_ratio=res_count_ratio
        out_struct[iout].dcount_diff=dcount_diff
        out_struct[iout].dcount_ratio=dcount_ratio

     endif

  endif
  
  if size(out_struct,/type) eq 8 and nsim eq 1 then begin ;;structure is given as input

     print, 'saving load_res parameters on already DEFINED structure nsim<10'
     print, 'iout = ',iout

     out_struct[iout].res = reform(res)
     if keyword_set(mres) then out_struct[iout].mres=mres
     if keyword_set(dres) then out_struct[iout].dres=dres
     ;;out_struct[iout].mcres = reform(mcres)

     if not keyword_set(nopeaks) then begin        
        ;;out_struct[iout].minmaxsad=minmaxsad
        out_struct[iout].res_count=res_count
        if keyword_set(mres) then out_struct[iout].mcount=mcount
        if keyword_set(dres) then out_struct[iout].dcount=dcount
     endif

  endif
  
  
  if size(out_struct,/type) ne 8 and (nsim gt 1) then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim>10'

     if not keyword_set(nopeaks) then begin        
        out_struct = {mcres:mcres,minmaxsad:minmaxsad,tminmaxsad:tminmaxsad, res:res, $
                      res_count:res_count, mres:mres,dres:dres,mcount:mcount,dcount:dcount,$
                      mtcount:mtcount,dtcount:dtcount,mcrit:mcrit,dcrit:dcrit,$
                      fmres:fmres,fdres:fdres, fmcount:fmcount,fdcount:fdcount,$
                      mcres_diff:mcres_diff,minmaxsad_diff:minmaxsad_diff,res_count_diff:res_count_diff,dcount_diff:dcount_diff,dres_diff:dres_diff,$
                      mcres_ratio:mcres_ratio,minmaxsad_ratio:minmaxsad_ratio,res_count_ratio:res_count_ratio,dcount_ratio:dcount_ratio,dres_ratio:dres_ratio}
     endif else begin
        out_struct = {mcres:mcres,res:res, mres:mres,dres:dres,fmres:fmres,fdres:fdres,$
                      mcres_diff:mcres_diff,dres_diff:dres_diff,mcres_ratio:mcres_ratio,dres_ratio:dres_ratio}
     endelse

     out_struct = replicate(out_struct,nrep)
  endif 
  
end
