pro load_res_varpol,fres,mcres,out_struct=out_struct,dobin=dobin,bins=bins,nzero=nzero, nrep=nrep, iout=iout,nrad=nrad,$
                      ires=ires,nsplit=nsplit,overwrite=overwrite,mres=mres,dres=dres,nsimuse=nsimuse,$
                      minnlist=minnlist,ivar=ivar,nproc=nproc, nodir=nodir,mask=mask,unnorm=unnorm,repires=repires

;;read spectra from cmb sim:

  if not keyword_set(ires) then ires=0
  if not keyword_set(nzero) then nzero=4
  if not keyword_set(minnlist) then minnlist=0.1
  if not keyword_set(nrep) then  nrep=1
  if not keyword_set(iout) then  iout=0
  if not keyword_set(ivar) then  ivar=3
  if not keyword_set(unnorm) then  unnorm=1

funnorm = ''
if unnorm ne 1 then funnorm='_unnorm'+strn(unnorm)

;print,'out_struct defined or not: type=',size(out_struct,/type)
;help, out_struct

  infores = size(mcres)
  print, 'infro mcres', infores

  nspots=n_elements(mcres[*,0,0,0])
  nbins=n_elements(mcres[0,*,0,0])
  nrad = n_elements(mcres[0,0,*,0])
  nsim = n_elements(mcres[0,0,0,*])

  if not keyword_set(nproc) then  nproc = min([nsim,50])

                                ; if infores[0] gt 2 then nsim=infores[3]
  if not keyword_set(nsimuse) then nsimuse=nsim


  res0=dblarr(nspots,nbins,nrad)

  count=0
  for i=0l,nsim-1l do begin
     simno= strn(i)             ;trim(string(i),2)
                                ;simno=str_replicate('0',nzero-strlen(simno))+simno
     file=fres+'.unf_insim'+simno
     if file_test(file,/regular) then begin
        print, 'readin: '+file
        runf,res0,file
        
        
        mcres[*,*,*,count]=res0
        count=count+1
     endif else print, '****missing file skipped ***',file
     
  endfor

if count eq 0 then begin
   print, 'load_res_varpol: no file is read. Quitting ..'
   on_error, 2
endif

;;in case missing files use only those that are available
if nsim ne count then begin
   print
   print, 'nsim, nsim_read',nsim, nsim_count
   nsim=count
   mcres=mcres[*,*,*,0:nsim-1]
endif

;;;;;this is set by hand here because the simulation number 2 for
;;;;;planck ffp6 is set to be noise only (no CMB)


  nmcres = reform(mcres[*,ivar,*,*])
  res = reform(nmcres[*,*,ires])
  dmcres=res

  help, res
  help, mcres

nside_map = 512 ; npix2nside(nspots)

  if nsim gt 10 then begin

     
     decompose,fres,d,path,name,ext,ver
     ext='.unf'
     dirout=path+'/post_nrad'+strn(nrad)+'_fnlist0'+strn(round(100*minnlist))+'/'

     print,'dirout = '+dirout
     spawn, 'mkdir -p '+dirout
     
     fmres = dirout+'mres_'+name+ext
     fdres = dirout+'dres_'+name+ext
     fmask = dirout+'res_mask_'+name+ext

     if  keyword_set(mres) then mres_use=mres    
     if  keyword_set(dres) then dres_use=dres    

     mres=dblarr(nspots,nrad)
     dres=mres
     mask=fltarr(nspots,nrad)
     max_nlist=dblarr(nrad)
     for i=0,nrad-1 do begin              
        max_nlist[i] = max(mcres[*,0,i,*])
     endfor

     print, 'max_nlist',max_nlist     

     print, 'max nlist of every radius patch: ',round(max_nlist)

     if  not file_test(fmres,/regular) or keyword_set(overwrite) then begin

        for j=0,nspots-1 do begin

           if j mod 100 eq 0 then print,'mean var nspot = ',j
           
           for i=0,nrad-1 do begin              

               ;pix2vec_ring,nside_map,j,vec
               ;query_disc,nside_map,vec,mcres[j,1,i,0],listpix
               ;print, i, j, n_elements(listpix)

              if mcres[j,0,i,0] gt 50 then begin ; 
               ;if mcres[j,0,i,0] gt minnlist*n_elements(listpix) then begin
                 mask[j,i]=1.
                 mres(j,i)=mean(nmcres(j,i,0:nsimuse-1))
                 dres(j,i)=stdev(nmcres(j,i,0:nsimuse-1))
                 
              endif else  begin
                 mask[j,i]=0.
                 mres(j,i)=0.
                 dres(j,i)=1e10
                 nmcres[j,i,*] = 0.
              endelse

           endfor
        endfor

        wunf, mres, fmres
        wunf, dres, fdres
        wunf, mask, fmask
     endif else begin
        print, 'loading  ',fmres
        runf, mres, fmres
        print, 'loading  ',fdres
        runf, dres, fdres
        print, 'loading  ',fmask
        runf, mask, fmask
        
     endelse

     if  keyword_set(mres_use) and keyword_set(dres_use) then begin 
        print, 'mres, dres, mres_use, dres_use: ',total(mres),total(dres),total(mres_use),total(dres_use)
        print, 'diff_mres, diff_dres: ',total(mres-mres_use),total(dres-dres_use)
        mres=mres_use
     endif



     print, 'total_mres, total_dres, total_mcres', mean(mres), mean(dres), mean(mcres)
     
     theta=dblarr(nrad,nsim+1)
     phi=theta
     dipamp=theta

     dtheta = reform(theta[*,0])
     dphi = reform(phi[*,0])
     ddipamp = reform(dipamp[*,0])


     ftheta = dirout+'theta_'+name+funnorm+ext
     fphi = dirout+'phi_'+name+funnorm+ext
     fdipamp = dirout+'dipamp_'+name+funnorm+ext


     if not keyword_set(nodir) then  begin


        if not file_test(ftheta,/regular) or keyword_set(overwrite) then begin ; 
           get_variance_dipdir,res,dres,nmcres,mres,$
                               nproc=nproc,mask=mask,$
                               theta=theta,phi=phi,bres=bres,dipamp=dipamp,unnorm=unnorm
           
           wunf, theta, ftheta
           wunf, phi, fphi    
           wunf,dipamp,fdipamp
           
        endif else begin
           
           print, 'loading  theta ',ftheta
           runf, theta, ftheta
           print, 'loading  ',fphi
           runf, phi, fphi
           print, 'loading  ',fdipamp
           runf, dipamp, fdipamp        
        endelse
        
        dtheta = reform(theta[*,0])
        dphi = reform(phi[*,0])
        ddipamp = reform(dipamp[*,0])
        
        theta = theta[*,1:nsim]
        phi = phi[*,1:nsim]
        dipamp = dipamp[*,1:nsim]
        
     endif   
  endif

  if (size(out_struct,/type) eq 8) and (nsim gt 10) then begin ;;structure is given as input
     print, 'saving load_res parameters on already DEFINED structure: nsim>10'
     out_struct[iout].res=res
     out_struct[iout].mcres=mcres
     out_struct[iout].mres=mres
     out_struct[iout].dres=dres
     out_struct[iout].fmres=fmres
     out_struct[iout].fdres=fdres
     out_struct[iout].fmask=fmask
     out_struct[iout].mask=mask
     
     out_struct[iout].theta=theta
     out_struct[iout].ftheta=ftheta
     out_struct[iout].phi=phi
     out_struct[iout].fphi=fphi
     out_struct[iout].dipamp=dipamp
     out_struct[iout].fdipamp=fdipamp
     
     out_struct[iout].ddipamp=ddipamp

  endif
  
  if size(out_struct,/type) eq 8 and nsim eq 1 then begin ;;structure is given as input

     print, 'saving load_res parameters on already DEFINED structure nsim<10'
     print, 'iout = ',iout

     mask = out_struct[iout].mask

     out_struct[iout].res = reform(res)
     out_struct[iout].dmcres = reform(nmcres)     
     if keyword_set(mres) then out_struct[iout].mres=mres
     if keyword_set(dres) then out_struct[iout].dres=dres

     get_variance_dipdir,res,out_struct[iout].dres,nmcres,out_struct[iout].mres,$
                         nproc=2,mask=mask,$
                         theta=theta,phi=phi,bres=bres,dipamp=dipamp,unnorm=unnorm
     

     out_struct[iout].dtheta = reform(theta[*,0])
     out_struct[iout].dphi = reform(phi[*,0])            
     out_struct[iout].ddipamp = reform(dipamp[*,0])
  endif
  
  
  if size(out_struct,/type) ne 8 and (nsim gt 10) then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim>10'
     out_struct = {dmcres:dmcres, res:res, mcres:mcres, mres:mres,dres:dres,mask:mask,fmres:fmres,fdres:fdres, fmask:fmask,theta:theta,phi:phi,dipamp:dipamp,ftheta:ftheta,fphi:fphi,fdipamp:fdipamp,dtheta:dtheta,dphi:dphi,ddipamp:ddipamp}
     out_struct = replicate(out_struct,nrep)
  endif 
  
end
