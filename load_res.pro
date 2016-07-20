pro load_res,fres,mcres,out_struct=out_struct,dobin=dobin,bins=bins,nzero=nzero, nrep=nrep, iout=iout,$
             ires=ires,nsplit=nsplit,nband=nband,nnbin=nnbin,overwrite=overwrite,mres=mres,dres=dres,nsimuse=nsimuse,nproc=nproc,unnorm=unnorm

;;read spectra from cmb sim:

  if not keyword_set(ires) then ires=0
  if not keyword_set(nzero) then nzero=4
  if not keyword_set(nband) then nband=15l
  if not keyword_set(nnbin) then  nnbin=6l
  if not keyword_set(nrep) then  nrep=1
  if not keyword_set(iout) then  iout=0
  if not keyword_set(unnorm) then unnorm=1

;print,'out_struct defined or not: type=',size(out_struct,/type)
;help, out_struct

  infores = size(mcres)
  nbins=infores[1]
  nspots=infores[2]

  if not keyword_set(bins) then  bins = indgen(nbins)*16l+2l

  nsim=1
  if infores[0] gt 2 then nsim=infores[3]
  if not keyword_set(nsimuse) then nsimuse=nsim

  if not keyword_set(nproc) then  nproc = min([nsim+1,50])

print, 'load_res: nbins, nspots, nsim, nsimuse: ',nbins, nspots, nsim,nsimuse

  if keyword_set(nsplit) then begin

     res0=dblarr(nbins,3072,nsplit)
     nblock=nsim/nsplit

     for i=0l,nblock-1l do begin
        simno=trim(string(i),2)
        simno=str_replicate('0',nzero-strlen(simno))+simno
        file=fres+simno+'.unf'
        print, 'readin: '+file
        runf,res0,file
        mcres(*,*,i*nsplit:(i+1l)*nsplit-1l)=res0
     endfor

  endif else runf,mcres,fres

  print, 'mean(mcres): ',mean(mcres)


  if nsimuse eq 20 then begin
     ;;TEMPORARYYYYYY set NAN sims
     ;; mcres[*,*,3]=mcres[*,*,24]
     ;; mcres[*,*,7]=mcres[*,*,23]
     
     for ii=0,nsimuse-1 do begin
        print, 'sim, mean(mcres): ',ii,mean(mcres[*,*,ii])
     endfor
     
  endif
  
  mcres = mcres[*,*,0:nsimuse-1]
  nsim=nsimuse
  
  print, 'load res: unnorm = ',unnorm

  fnorm=''
  if unnorm gt 1 then fnorm='_norm'+strn(unnorm)


  if nsim gt 10 then begin

     
     decompose,fres,d,path,name,ext,ver
     dirout=path+'/post_nband'+strn(nband)+'/'

     print,'dirout = '+dirout
     spawn, 'mkdir -p '+dirout
     
     fmres = dirout+'mres_'+name+ext
     fdres = dirout+'dres_'+name+ext

     if not keyword_set(mres) then begin     

        mres=dblarr(nbins,nspots)
        dres=mres
        bres=dblarr(nband,nspots)

        if  not file_test(fmres,/regular) or keyword_set(overwrite) then begin
           for i=0,nbins-1 do begin
              if i mod 20 eq 0 then print,i
              for j=0,nspots-1 do begin
                 mres(i,j)=mean(mcres(i,j,0:nsimuse-1))
                 dres(i,j)=stdev(mcres(i,j,0:nsimuse-1))
              endfor
           endfor
           wunf, mres, fmres
           wunf, dres, fdres
        endif else begin
           print, 'loading  ',fmres
           runf, mres, fmres
           print, 'loading  ',fdres
           runf, dres, fdres
           
        endelse
     endif else print, 'diff_mres, diff_dres: ',total(out_struct[0].mres-mres),total(out_struct[0].dres-dres)

     res = reform(mcres[*,*,ires])
     
     
     theta=dblarr(nband,nsim+1)
     phi=theta
     dipamp=theta
     ftheta = dirout+'theta_'+name+fnorm+ext
     fphi = dirout+'phi_'+name+fnorm+ext
     fdipamp = dirout+'dipamp_'+name+fnorm+ext
     
     if not file_test(ftheta,/regular) or keyword_set(overwrite) then begin ; 
        get_cross_dipdir,bins,res,dres,mcres,mres,$
                         nband=nband,nnbin=nnbin,nproc=nproc,$
                         theta=theta,phi=phi,bres=bres,dipamp=dipamp,unnorm=unnorm
        
        wunf, theta, ftheta
        wunf, phi, fphi    
        wunf,dipamp,fdipamp
        
     endif else begin
        
        print, 'loading  theta ',ftheta
        runf, theta, ftheta
        print, 'loading  ',fphi
        runf, phi, fphi
        print, 'loading ',fdipamp
       if file_test(fdipamp,/regular) then  runf,dipamp,fdipamp
        
     endelse
     
     dtheta = reform(theta[*,0])
     dphi = reform(phi[*,0])
     ddipamp = reform(dipamp[*,0])

     theta = theta[*,1:nsim]
     phi = phi[*,1:nsim]
     dipamp = dipamp[*,1:nsim]

  endif   
  

  if (size(out_struct,/type) eq 8) and (nsim gt 10) then begin ;;structure is given as input
     print, 'saving load_res parameters on already DEFINED structure: nsim>10'
     out_struct[iout].mres=mres
     out_struct[iout].dres=dres
     out_struct[iout].fmres=fmres
     out_struct[iout].fdres=fdres
     
     out_struct[iout].theta=theta
     out_struct[iout].ftheta=ftheta
     out_struct[iout].phi=phi
     out_struct[iout].fphi=fphi

     out_struct[iout].dipamp=dipamp     
     out_struct[iout].fdipamp=fdipamp
     out_struct[iout].ddipamp=ddipamp

     out_struct[iout].dtheta=dtheta
     out_struct[iout].dphi=dphi
  endif
  
  if size(out_struct,/type) eq 8 and nsim eq 1 then begin ;;structure is given as input

     print, 'saving load_res parameters on already DEFINED structure nsim<10'
     
     res = mcres[*,*,ires]
     mcres = res
     
     
     get_cross_dipdir,bins,res,out_struct[iout].dres,mcres,out_struct[iout].mres,$
                      nband=nband,nnbin=nbin2bin,nproc=2,$
                      theta=theta,phi=phi,bres=bres,dipamp=dipamp,unnorm=unnorm
     
     out_struct[iout].dtheta = reform(theta[*,0])
     out_struct[iout].dphi = reform(phi[*,0])            
     out_struct[iout].ddipamp = reform(dipamp[*,0])       
  endif
  
  
  if size(out_struct,/type) ne 8 and (nsim gt 10) then begin ;;structure is given as input
     print, 'saving load_res DEFINING NEW structure: nsim>10'
     out_struct = {mres:mres,dres:dres,fmres:fmres,fdres:fdres, theta:theta,phi:phi,dipamp:dipamp,ftheta:ftheta,fphi:fphi,fdipamp:fdipamp,dtheta:dtheta,dphi:dphi,ddipamp:ddipamp}
     out_struct = replicate(out_struct,nrep)
  endif 
  
end
