pro dipdir_lmin_lmax,res, mcres, mres,dres, dirroot=dirroot,method_mask=method_mask,binvec=binvec,$
                     nside=nside,lmax=lmax,overwrite=overwrite,aps=aps,$
                     mctheta=mctheta,mcphi=mcphi,mcdipamp=mcdipamp,nproc=nproc,dskip=dskip,noplot=noplot,unnorm=unnorm
  


  if not keyword_set(aps) then aps = 0 & daps=0
  if not keyword_set(dirroot) then dirroot = '../output/ddx9/'
  if not keyword_set(method_mask) then method_mask = ''
  if not keyword_set(nproc) then nproc=125
  if not keyword_set(dskip) then dskip=1

  nband=15
  dipole_nside=512l


  if not keyword_set(binvec) then  begin
     runf,binvec,'/mn/owl1/d3/yabebalf/planck/DX9/bins_150_lmax2400.unf',prec=1
  endif

if not keyword_set(unnorm) then unnorm=1


  if n_params() eq 1 then begin
     mcres = res
     dres=res
     mres=res

     unnorm=2
  endif
  
  bins=binvec
  get_bin_dipdir,val,binmin,binmax,show=show,bvec=bvec

  nsim = n_elements(mcres[0,0,*])
  nmask = 1

  nbin2bin = 6l

  mctheta = dblarr(nband,nsim+1,nband) 
  mcphi = mctheta
  mcdipamp = mctheta

  llmin=2

  allband = 'oneband/'+method_mask+'/lmin_func/'
  dirout = dirroot+allband
  spawn, 'mkdir -p '+dirout


  count3=0
  for lastb = nband-1,0,-1*dskip do begin


     nbj = bvec[lastb,1]      ;;the last dl=16 band to consider

     nband_lastb = lastb+1 ;;this is the number of lmin to consider 

     alltheta=dblarr(nband_lastb,nsim+1)
     allphi=alltheta
     dipamp=alltheta

     llmax = lastb*100+100

     print, 'llmin, llmax, Aps,nbj,nbin2bin,nband_lastb : ',llmin, llmax, aps,nbj,nbin2bin,nband_lastb
     

     falltheta = dirout+'alltheta'+'_lmax'+strn(llmax)+'.unf'
     fallphi = dirout+'allphi'+'_lmax'+strn(llmax)+'.unf'
     fdipamp = dirout+'dipamp'+'_lmax'+strn(llmax)+'.unf'


     if not file_test(falltheta,/regular) or keyword_set(overwrite) then begin ; 
        get_cross_lmin_dipdir,bins,res,dres,mcres,mres,lmax=lmax,nbj=nbj,$
                              nside=nside,nband=nband_lastb,aps=aps,nproc=nproc,nnbin=nbin2bin,$
                              theta=alltheta,phi=allphi,bres=bres,load=load,unnorm=unnorm,dipamp=dipamp
        
        wunf, alltheta, falltheta
        wunf, allphi, fallphi    
        wunf, dipamp, fdipamp 
        
        
     endif else begin
        print, 'loading  ',falltheta
        runf, alltheta, falltheta
        print, 'loading  ',fallphi
        runf, allphi, fallphi
        print, 'loading  ',fdipamp
        runf, dipamp, fdipamp

     endelse

;help, alltheta
;print, nband_lastb

     mctheta(0:nband_lastb-1,0:nsim,count3) = alltheta(0:nband_lastb-1,0:nsim)
     mcphi(0:nband_lastb-1,0:nsim,count3) = allphi(0:nband_lastb-1,0:nsim)
     mcdipamp(0:nband_lastb-1,0:nsim,count3) = dipamp(0:nband_lastb-1,0:nsim)

     count3 = count3+1
  endfor
  ncount3 = count3-1
  print,'ncount3 = ',ncount3

  falltheta = dirout+'all_case_theta'+'.unf'
  fallphi = dirout+'all_case_phi'+'.unf'
  falldipamp = dirout+'all_case_dipamp'+'.unf'
  wunf,mctheta,falltheta
  wunf,mcphi,fallphi
  wunf,mcdipamp,falldipamp


  figdir = dirout+'figure/'
  spawn, 'mkdir -p '+figdir


  val_max = indgen(15)*100+100

  if not keyword_set(noplot) then begin

     count3 = 0
     for lastb = nband-1,0,-1*dskip do begin

        nband_lastb = lastb+1 ;;this is the number of lmin to consider 

        theta = mctheta(0:nband_lastb-1,0,count3)
        phi = mcphi(0:nband_lastb-1,0,count3)

        dipamp = mcdipamp(0:nband_lastb-1,0,count3)

        llmax = val_max[lastb]
        val = [2,val_max[0:lastb-1]]
        map = make_disk_to_mollview(theta, phi,val,2d0,dipole_nside)
        ;;plot dipole directions 
        plot_mollview,map,2,max([2,lastb*100]),file_ps=figdir+'asymm_'+method_mask+'_lmin_func_lmax'+strn(llmax)+'.ps',$
                      /show_dipdir,title=' ',mtheta=theta,mphi=phi,midval=0.5,ctitle='subttl',tkey='lttl',ax=['subttl','lttl'],ay=['$\ell_{\rm min}$','$\ell_{\rm max} = '+strn(llmax)+'$'],val=val


        count3 = count3 + 1
     endfor


     ;; 

     count3 = 0
     for lastb = nband-1,0,-1*dskip do begin


        theta = mctheta(count3,0,0:(ncount3-count3))
        phi = mcphi(count3,0,0:(ncount3-count3))

        ;;print, 'lmin func theta: ',theta
        ;;print, 'lmin func phi: ',phi

        theta = reverse(theta)
        phi = reverse(phi)

        ;;print,'theta and phi ]', [[theta],[phi]]

        llmin = round(max([2,count3*100*dskip]))

        ellsize=0.4
        if count3 lt 2 then ellsize = 0.2

        val = reverse(val_max[nband-1:(count3*dskip):-dskip])
        print, val
        map = make_disk_to_mollview(theta, phi,val,2d0,dipole_nside)
        ;;plot dipole directions 
        plot_mollview,map,min(val),max(val),file_ps=figdir+'asymm_'+method_mask+'_lmax_func_lmin'+strn(llmin)+'.ps',$
                      /show_dipdir,title=' ',mtheta=theta,mphi=phi,midval=0.5,ctitle='subttl',tkey='lttl',ax=['subttl','lttl'],ay=['$\ell_{\rm max}$','$\ell_{\rm min} = '+strn(llmin)+'$'],ellsize=ellsize,val=val


        ;; if lastb lt 2 then begin
        ;;    plot_gnomview,map,val(0),1500,file_ps=figdir+'asymm_gnome_'+method_mask+'_lmax_func_lmin'+strn(llmin)+'_aps'+strn(daps)+'.ps',$
        ;;                  /show_dipdir,title=' ',mtheta=theta,mphi=phi,midval=0.5,ctitle='subttl',tkey='lttl',ax=['subttl','lttl'],ay=['$\ell_{\rm max}$','$\ell_{\rm min} = '+strn(llmin)+'$'],ellsize=ellsize
        ;; endif

        count3 = count3 + 1
     endfor

  endif
  

end
