pro load_res_varneedle,fres,mcres,out_struct=out_struct,dobin=dobin,bins=bins,nzero=nzero, nrep=nrep, iout=iout,nrad=nrad,$
                       ires=ires,nsplit=nsplit,overwrite=overwrite,mres=mres,dres=dres,nsimuse=nsimuse,sfactor=sfactor,dfactor=dfactor,$
                       minnlist=minnlist,ivar=ivar,nproc=nproc, nodir=nodir,mask=mask,unnorm=unnorm,repires=repires,nomcres=nomcres

;;read spectra from cmb sim:

  if not keyword_set(ires) then ires=0
  if not keyword_set(nzero) then nzero=4
  if not keyword_set(minnlist) then minnlist=0.1
  if not keyword_set(nrep) then  nrep=1
  if not keyword_set(iout) then  iout=0
  if not keyword_set(ivar) then  ivar=3
  if not keyword_set(unnorm) then  unnorm=1
  if not keyword_set(dfactor) then  dfactor=1.
  if not keyword_set(sfactor) then  sfactor=1.

  funnorm = ''
  if unnorm ne 1 then funnorm='_unnorm'+strn(unnorm)

;print,'out_struct defined or not: type=',size(out_struct,/type)
;help, out_struct

  infores = size(mcres)
  print, 'infro mcres', infores

  nspots=n_elements(mcres[*,0,0,0,0])
  nbins=n_elements(mcres[0,*,0,0,0])
  nrad = n_elements(mcres[0,0,*,0,0])
  nj = n_elements(mcres[0,0,0,*,0])
  nsim = n_elements(mcres[0,0,0,0,*])

  if not keyword_set(nproc) then  nproc = min([nsim+1,50])

                                ; if infores[0] gt 2 then nsim=infores[3]
  if not keyword_set(nsimuse) then nsimuse=nsim

;;  if keyword_set(nsplit) then begin

  res0=dblarr(nspots,nbins,nrad,nj)

  if not keyword_set(nomcres) then begin  
     for i=0l,nsim-1l do begin
        simno= strn(i)          ;trim(string(i),2)
                                ;simno=str_replicate('0',nzero-strlen(simno))+simno
        file=fres+'.unf_insim'+simno
        if i lt 10 then print, 'readin: '+file
        ;;if total(i eq [205,277, 326, 400, 783, 804, 856, 947, 956, 969, 988]) eq 0 then begin
        runf,res0,file
        ;;endif else print, 'not readying faulty file i',i
        if i eq ires then res0[*,ivar,*,*] = dfactor*res0[*,ivar,*,*]
        if i eq ires then print, 'dfactor, var(6deg): ',dfactor, reform(res0[10,ivar,*,*])
        if i ne ires then res0[*,ivar,*,*] = res0[*,ivar,*,*]/sfactor
        
        mcres[*,*,*,*,i]=res0
     endfor
  endif


;print,'data', mcres[0:20,ivar,3,0,ires]
;print,'sim', mcres[0:20,ivar,3,0,ires+1]
;print,'sim', mcres[0:20,ivar,3,0,nsim-10]


;;print, 'mean(data)=',mean(mcres[0:10,ivar,*,*,ires])
;;print, 'mean(sim)=',mean(mcres[0:10,ivar,*,*,ires+1])


;;poorman data/sim unit control 
  if nsim gt ires+1 then begin
     ;;if data values are very smaller or very high from sims stop

     msim=mean(mcres[0:10,ivar,nrad-1,*,ires+1])
     mdat=mean(mcres[0:10,ivar,nrad-1,*,ires])

     if mdat gt msim*1e6 then begin
        print, 'dfactor, sfactor',dfactor, sfactor
        print, 'mean(data)=',mdat
        print, 'mean(sim)=',msim
        message,fres+': may be data values are much higher than simulations. dfactor=1e-12'
     endif
     if mdat*1e6 lt msim then begin
        print, 'dfactor, sfactor',dfactor, sfactor
        print, 'mean(data)=',mdat
        print, 'mean(sim)=',msim
        message,fres+': may be data values are much smaller than simulations. Use dfactor=1e12'
     endif
  endif


;;;;;this is set by hand here because the simulation number 2 for
;;;;;planck ffp6 is set to be noise only (no CMB)
  if keyword_set(repires) then mcres[*,*,*,*,repires]=mcres[*,*,*,*,repires+1]

  nmcres = mcres[*,ivar,*,*,*]

  ;;make data the 0th row and don't use it in mean/var
  if ires ne 0 then begin
     zz=nmcres[*,0,*,*,0]
     nmcres[*,0,*,*,0] = nmcres[*,0,*,*,ires]
     nmcres[*,0,*,*,ires] = zz
     ires=0
  endif

  ;;weighted variance 
  ;;wmcres=mcres

  res = reform(nmcres[*,0,*,*,ires])
  dmcres=res

  help, res
  help, mcres

  nside_map = 512               ; npix2nside(nspots)

;;number of pixels per disk
  fpixdisk='/mn/owl1/d3/yabebalf/planck/common_files/n_pixels_nside2048_nspot'+strn(nspots)+'_nrad16.unf'
  npix_per_disk = lindgen(nspots, nrad)
  runf, npix_per_disk, fpixdisk


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

     mres=dblarr(nspots,nrad,nj)
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

              if mcres[j,0,i,0] gt minnlist*mcres[j,1,i,0]  then begin ;;npix_per_disk[j,i]

                 mask[j,i]=1.
                 for jned=0,nj-1 do begin
                    mres(j,i,jned)=mean(nmcres(j,0,i,jned,1:nsimuse-1))
                    dres(j,i,jned)=stdev(nmcres(j,0,i,jned,1:nsimuse-1))
                 endfor
              endif else  begin
                 mask[j,i]=0.
                 mres(j,i,*)=1e10
                 dres(j,i,*)=1e10
                 nmcres[j,0,i,*,*] = 1e10
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

     print, 'mcres,mres, dres: ',total(mcres),total(mres),total(dres)
     if  keyword_set(mres_use) and keyword_set(dres_use) then begin 
        print, 'mres, dres, mres_use, dres_use: ',total(mres),total(dres),total(mres_use),total(dres_use)
        print, 'diff_mres, diff_dres: ',total(mres-mres_use),total(dres-dres_use)
        mres=mres_use
     endif


     print, 'total_mres, total_dres, total_mcres', mean(mres), mean(dres), mean(mcres)
     
     theta=dblarr(nrad,nj,nsim)
     phi=theta
     dipamp=theta
     
     dtheta = dblarr(nrad,nj)
     dphi = dtheta
     ddipamp = dtheta

     xxtheta = dblarr(nrad,nsim+1)
     xxphi = xxtheta
     xxdipamp = xxtheta

     ftheta = dirout+'theta_j'+strn(0)+'_'+name+funnorm+ext
     fphi = dirout+'phi_j'+strn(0)+'_'+name+funnorm+ext
     fdipamp = dirout+'dipamp_j'+strn(0)+'_'+name+funnorm+ext


     if not keyword_set(nodir) then  begin
        
        for jned=0,nj-1 do begin

           ftheta = dirout+'theta_j'+strn(jned)+'_'+name+funnorm+ext
           fphi = dirout+'phi_j'+strn(jned)+'_'+name+funnorm+ext
           fdipamp = dirout+'dipamp_j'+strn(jned)+'_'+name+funnorm+ext
           

           xxres=res[*,*,jned]
           xxdres = dres[*,*,jned]
           xxmres = mres[*,*,jned]
           xxnmcres = reform(nmcres[*,0,*,jned,*])

           if not file_test(ftheta,/regular) or keyword_set(overwrite) then begin ; 
              get_variance_dipdir,xxres,xxdres,xxnmcres,xxmres,$
                                  nproc=nproc,mask=mask,$
                                  theta=xxtheta,phi=xxphi,bres=xxbres,dipamp=xxdipamp,unnorm=unnorm
              
              wunf, xxtheta, ftheta
              wunf, xxphi, fphi    
              wunf,xxdipamp,fdipamp
              
           endif else begin
              
              print, 'loading  theta ',ftheta
              runf,xxtheta, ftheta
              print, 'loading  ',fphi
              runf, xxphi, fphi
              print, 'loading  ',fdipamp
              runf, xxdipamp, fdipamp        
           endelse
           
           dtheta[*,jned] = reform(xxtheta[*,0])
           dphi[*,jned] = reform(xxphi[*,0])
           ddipamp[*,jned] = reform(xxdipamp[*,0])
           
           theta[*,jned,0:nsim-1] = xxtheta[*,1:nsim]
           phi[*,jned,0:nsim-1] = xxphi[*,1:nsim]
           dipamp[*,jned,0:nsim-1] = xxdipamp[*,1:nsim]
           
        endfor
        
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

     for jned=0,nj-1 do begin
        get_variance_dipdir,res[*,*,jned],out_struct[iout].dres[*,*,jned],res[*,*,jned],out_struct[iout].mres[*,*,jned],$
                            nproc=2,mask=mask,$
                            theta=xxtheta,phi=xxphi,bres=xxbres,dipamp=xxdipamp,unnorm=unnorm
        

        help, res
        help, xxtheta
        help, out_struct[0].dtheta
        print, jned
        out_struct[iout].dtheta[*,jned] = reform(xxtheta[*,0])
        out_struct[iout].dphi[*,jned] = reform(xxphi[*,0])            
        out_struct[iout].ddipamp[*,jned] = reform(xxdipamp[*,0])
     endfor
  endif
  
  
  if size(out_struct,/type) ne 8 and (nsim gt 10) then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim>10'
     out_struct = {dmcres:dmcres, res:res, mcres:mcres, mres:mres,dres:dres,mask:mask,fmres:fmres,fdres:fdres,$ 
                   fmask:fmask,theta:theta,phi:phi,dipamp:dipamp,ftheta:ftheta,fphi:fphi,fdipamp:fdipamp,$
                   dtheta:dtheta,dphi:dphi,ddipamp:ddipamp}
     out_struct = replicate(out_struct,nrep)
  endif 
  
end
