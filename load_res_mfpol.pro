pro load_res_mfpol,fres,mcres,minmaxsad,out_struct=out_struct,dobin=dobin,nzero=nzero, nrep=nrep, iout=iout,$
                      ires=ires,nsplit=nsplit,overwrite=overwrite,mres=mres,dres=dres,mcount=mcount,dcount=dcount,nsimb=nsimb,nsimf=nsimf,$
                      nproc=nproc, nodir=nodir,fkey=fkey,pfres=pfres,mpow=mpow,norm=norm,mfnorm=mfnorm

;;read spectra from cmb sim:

  if not keyword_set(ires) then ires=0
  if not keyword_set(nzero) then nzero=4
  if not keyword_set(nrep) then  nrep=1
  if not keyword_set(iout) then  iout=0
  if not keyword_set(fkey) then fkey=''
  if not keyword_set(mpow) then mpow=1


;print,'out_struct defined or not: type=',size(out_struct,/type)
;help, out_struct

  infores = size(mcres)
  print, 'infro mcres', infores

  nsteps=n_elements(mcres[*,0,0,0,0])
  nmf=n_elements(mcres[0,*,0,0,0])
  nj = n_elements(mcres[0,0,*,0,0])
  npol = n_elements(mcres[0,0,0,*,0])
  nsim = n_elements(mcres[0,0,0,0,*])

  nsim_orig=nsim

  if not keyword_set(nsimb) then nsimb=1
  if not keyword_set(nsimf) then nsimf=nsim
  nsim=nsimf-nsimb+1


  mcres = dblarr(nsteps,nmf,nj,npol,nsim)
  res0=dblarr(nsteps,nmf,nj,npol)
  
  minmaxsad=dblarr(nsteps,3,nj,npol,nsim)
  tminmaxsad=lindgen(1,3,nj,npol,nsim)
  norm_count=dblarr(nsteps,3,nj,npol)
  res0_count=lindgen(nsteps+1,3,nj,npol)

  print, 'nsim, nj,nsteps ',nsim, nj, nsteps

  count=0
  for i=0,nsim-1l do begin
     isim=nsimb+i-1

     simno= strn(isim)             ;trim(string(i),2)
                                ;simno=str_replicate('0',nzero-strlen(simno))+simno
     file=fres+'.unf_insim'+simno
     
     ;;print, 'reading ',file

     if file_test(file,/regular) then begin

        if not keyword_set(pfres) then begin

           file=fres+'.unf_insim'+simno
           runf,res0,file          
           
           file=fres+'.unf_minmaxsad_insim'+simno
           runf,res0_count,file
        endif else begin

           restp=res0[*,*,*,1]

           file=fres+'.unf_insim'+simno
           runf,restp,file
           res0[*,*,*,0]=restp

           ;;print,'temp: ',total(restp)
           pfile=pfres+'.unf_insim'+simno          
           runf,restp,pfile 
           res0[*,*,*,1]=restp
           ;;print,'pol: ',total(restp)

           restp=res0_count[*,*,*,0]
           file=fres+'.unf_minmaxsad_insim'+simno
           runf,restp,file
           res0_count[*,*,*,0]=restp

           pfile=pfres+'.unf_minmaxsad_insim'+simno
           runf,restp,pfile
           res0_count[*,*,*,1]=restp
     endelse

        ;;normalized MF
        if keyword_set(mfnorm) then begin
           for ij=0,nj-1 do begin
              res0[*,1,ij,*]=mfnorm[ij,0]*res0[*,1,ij,*]    ;;area
              res0[*,3,ij,*]=mfnorm[ij,1]*res0[*,3,ij,*]    ;;length
              res0[*,2,ij,*]=mfnorm[ij,2]*res0[*,2,ij,*]    ;;genus
           endfor
        endif
        


        ;; if keyword_set(norm) then begin
        ;;    if mpow eq 2 then begin
        ;;       for jj=1,nmf-1 do begin                               
        ;;          for kk=0,nj-1 do begin              
        ;;             res0[(nsteps/2+10):nsteps-1,jj,kk,0]=double(res0[(nsteps/2+10):nsteps-1,jj,kk,0])/max(res0[(nsteps/2+10):nsteps-1,jj,kk,0]*1d0)
        ;;             res0[(nsteps/2+10):nsteps-1,jj,kk,1]=double(res0[(nsteps/2+10):nsteps-1,jj,kk,1])/max(res0[(nsteps/2+10):nsteps-1,jj,kk,1]*1d0)
        ;;          endfor
        ;;       endfor
        ;;    endif else begin
        ;;       for jj=1,nmf-1 do begin                               
        ;;          for kk=0,nj-1 do begin              
        ;;             res0[*,jj,kk,0]=double(res0[*,jj,kk,0])/max(res0[*,jj,kk,0]*1d0)
        ;;             res0[*,jj,kk,1]=double(res0[*,jj,kk,1])/max(res0[*,jj,kk,1]*1d0)
        ;;          endfor
        ;;       endfor
        ;;    endelse
        ;; endif

        mcres[*,*,*,*,count]=res0

        
        norm_count=double(res0_count[1:nsteps,*,*,*,*])
        
        tminmaxsad[*,*,*,*,count]=res0_count[0,*,*,*,*]
        minmaxsad[*,*,*,*,count]=norm_count
        
        count=count+1
     endif

  endfor

nsim=count

print, 'load_res_mf: nsim_orig, nsim_read',nsim_orig,nsim

mcres = mcres[*,*,*,*,0:nsim-1]
minmaxsad=minmaxsad[*,*,*,*,0:nsim-1]

;;print,reform(res0_count[0,*,*,*])

  res =reform(mcres[*,*,*,*,ires])
  res_count = reform(minmaxsad[*,*,*,*,ires])

help, res
help, res_count
   
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

     fmcount = dirout+'mcount_'+name+ext
     fdcount = dirout+'dcount_'+name+ext



     mres=dblarr(nsteps,nmf,nj,npol)
     dres=mres

     mcount=dblarr(nsteps,3,nj,npol)
     dcount=mcount

     if  not file_test(fmres,/regular) or keyword_set(overwrite) then begin
        
        for j=0,nsteps-1 do begin
           
           if j mod 100 eq 0 then print,'mean var nsteps = ',j
              
           for k=0,nj-1 do begin              

              for z=0,npol-1 do begin

                 for i=0,nmf-1 do begin                               
                    mres(j,i,k,z)=mean(mcres(j,i,k,z,*))
                    dres(j,i,k,z)=stdev(mcres(j,i,k,z,*))
                 endfor
                 
                 for i=0,2 do begin                               
                    mcount(j,i,k,z)=mean(minmaxsad(j,i,k,z,*)*1.)
                    dcount(j,i,k,z)=stdev(minmaxsad(j,i,k,z,*)*1.)
                 endfor

              endfor

           endfor
        endfor
        
        wunf, mres, fmres
        wunf, dres, fdres

        wunf, mcount, fmcount
        wunf, dcount, fdcount
        
     endif else begin
        print, 'loading  ',fmres
        runf, mres, fmres
        print, 'loading  ',fdres
        runf, dres, fdres

        print, 'loading  ',fmcount
        runf, mcount, fmcount
        print, 'loading  ',fdcount
        runf, dcount, fdcount
        
     endelse

     print, 'sanity mres, dres: ',total(mres),total(dres)

  endif



  if (size(out_struct,/type) eq 8) and (nsim gt 1) then begin ;;structure is given as input
     print, 'saving load_res parameters on already DEFINED structure: nsim>10'
     out_struct[iout].res=res
     ;;out_struct[iout].mcres=mcres
     out_struct[iout].mres=mres
     out_struct[iout].dres=dres

     ;out_struct[iout].tminmaxsad=tminmaxsad
     out_struct[iout].res_count=res_count
     out_struct[iout].mcount=mcount
     out_struct[iout].dcount=dcount
     
  endif
  
  if size(out_struct,/type) eq 8 and nsim eq 1 then begin ;;structure is given as input

     print, 'saving load_res parameters on already DEFINED structure nsim<10'
     print, 'iout = ',iout

     out_struct[iout].res = reform(res)
     ;;out_struct[iout].mcres = reform(mcres)

     ;;out_struct[iout].minmaxsad=minmaxsad
     out_struct[iout].res_count=res_count
     
     if keyword_set(mres) then out_struct[iout].mres=mres
     if keyword_set(dres) then out_struct[iout].dres=dres
     if keyword_set(mres) then out_struct[iout].mcount=mcount
     if keyword_set(dres) then out_struct[iout].dcount=dcount


  endif
  
  
  if size(out_struct,/type) ne 8 and (nsim gt 1) then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim>10'
     out_struct = {mcres:mcres,minmaxsad:minmaxsad,tminmaxsad:tminmaxsad, res:res, res_count:res_count, mres:mres,dres:dres,mcount:mcount,dcount:dcount,fmres:fmres,fdres:fdres, fmcount:fmcount,fdcount:fdcount}
     out_struct = replicate(out_struct,nrep)
  endif 
  
end
