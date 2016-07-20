pro dipole_dir_lminmax,res, mcres, mres,dres, dirroot=dirroot,method_mask=method_mask,binvec=binvec,$
                     nside=nside,lmax=lmax,overwrite=overwrite,aps=aps,llmin=llmin,llmax=llmax,$
                     mctheta=mctheta,mcphi=mcphi,mcdipamp=mcdipamp,nproc=nproc,dskip=dskip,maxlblock=maxlblock,minlblock=minlblock,noplot=noplot,unnorm=unnorm
                     


  if not keyword_set(aps) then aps = 0 & daps=0
  if not keyword_set(dirroot) then dirroot = '../output/ddx9/'
  if not keyword_set(method_mask) then method_mask = ''
  if not keyword_set(nproc) then nproc=125
  if not keyword_set(dskip) then dskip=1
  if not keyword_set(llmin) then llmin=2
  if not keyword_set(llmax) then llmax=1500
  if not keyword_set(maxlblock) then maxlblock=93
  if not keyword_set(minlblock) then minlblock=0

  nband=15
  dipole_nside=512l
  lastb=nband-1

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

  mctheta = dblarr(nband,nsim+1,2) 
  mcphi = mctheta
  mcdipamp = mctheta


  allband = 'oneband/'+method_mask+'/lminmax_func/'
  dirout = dirroot+allband
  spawn, 'mkdir -p '+dirout


;;compute dipdir as a function of lmin upto lastb=14

  alltheta=dblarr(nband,nsim+1)
  allphi=alltheta
  dipamp=alltheta
  

  
  ;;-----------function of lmin---------------------------------------------  

  falltheta = dirout+'alltheta'+'_lmax'+strn(llmax)+'.unf'
  fallphi = dirout+'allphi'+'_lmax'+strn(llmax)+'.unf'
  fdipamp = dirout+'dipamp'+'_lmax'+strn(llmax)+'.unf'

  if not file_test(falltheta,/regular) or keyword_set(overwrite) then begin ; 
     nbj = maxlblock  ;;the last dl=16 band to consider


     get_cross_lmin_dipdir,bins,res,dres,mcres,mres,lmax=lmax,nbj=nbj,$
                           nside=nside,nband=nband,aps=aps,nproc=nproc,nnbin=nbin2bin,$
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
  

  mctheta(0:nband-1,0:nsim,0) = alltheta(0:nband-1,0:nsim)
  mcphi(0:nband-1,0:nsim,0) = allphi(0:nband-1,0:nsim)
  mcdipamp(0:nband-1,0:nsim,0) = dipamp(0:nband-1,0:nsim)

  ;;------------function of lmax--------------------------------------------

  falltheta = dirout+'alltheta'+'_lmin'+strn(llmin)+'.unf'
  fallphi = dirout+'allphi'+'_lmin'+strn(llmin)+'.unf'
  fdipamp = dirout+'dipamp'+'_lmin'+strn(llmin)+'.unf'


  if not file_test(falltheta,/regular) or keyword_set(overwrite) then begin ; 
     nbj = minlblock  ;;the first dl=16 band to consider

     get_cross_lmax_dipdir,bins,res,dres,mcres,mres,lmax=lmax,nbj=nbj,$
                           nside=nside,nband=nband,aps=aps,nproc=nproc,nnbin=nbin2bin,$
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
  

  mctheta(0:nband-1,0:nsim,1) = alltheta(0:nband-1,0:nsim)
  mcphi(0:nband-1,0:nsim,1) = allphi(0:nband-1,0:nsim)
  mcdipamp(0:nband-1,0:nsim,1) = dipamp(0:nband-1,0:nsim)
  ;;--------------------------------------------------------

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
  ;;------------- plot dipdir for fixed lmin and function of lmax
  theta = mctheta(0:nband-1,0,0)
  phi = mcphi(0:nband-1,0,0)
  dipamp = mcdipamp(0:nband-1,0,0)
     

  val = [2,val_max[0:nband-2]]
  map = make_disk_to_mollview(theta, phi,val,2d0,dipole_nside)
  ;;plot dipole directions 
  plot_mollview,map,2,max([2,lastb*100]),file_ps=figdir+'asymm_'+method_mask+'_lmin_func_lmax'+strn(llmax)+'.ps',$
                /show_dipdir,title=' ',mtheta=theta,mphi=phi,midval=0.5,ctitle='subttl',tkey='lttl',ax=['subttl','lttl'],ay=['$\ell_{\rm min}$','$\ell_{\rm max} = '+strn(llmax)+'$'],val=val


  ;;------------- plot dipdir for fixed lmax and function of lmin
  theta = mctheta(0:nband-1,0,1)
  phi = mcphi(0:nband-1,0,1)
  dipamp = mcdipamp(0:nband-1,0,1)
     

  val = val_max[0:nband-1]
  map = make_disk_to_mollview(theta, phi,val,2d0,dipole_nside)
  ;;plot dipole directions 
  plot_mollview,map,min(val),max([2,lastb*100]),file_ps=figdir+'asymm_'+method_mask+'_lmax_func_lmin'+strn(llmin)+'.ps',$
                /show_dipdir,title=' ',mtheta=theta,mphi=phi,midval=0.5,ctitle='subttl',tkey='lttl',ax=['subttl','lttl'],ay=['$\ell_{\rm min}$','$\ell_{\rm max} = '+strn(llmax)+'$'],val=val

  endif

end
