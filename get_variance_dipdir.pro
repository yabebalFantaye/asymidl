pro get_variance_dipdir,res,dres,mcres,mres,dipole_nside=dipole_nside,lmax=lmax,mask=mask,$
               nside=nside,firstb=firstb,nband=nband,aps=aps,nproc=nproc,nnbin=nnbin,$
               theta=theta,phi=phi,bres=bres,nbj=nbj,load=load,unnorm=unnorm,dipamp=dipamp

if n_params() lt 4  then begin
    print, 'USAGE: '
    print, 'get_dipdir,bins,res,dres,dnres,mcres,mres,mnres,nside=nside,lmax=lmax,$'
    print, '               nside_pixwin=nside_pixwin,nband=nband,aps=aps,nproc=nproc,$'
    print, '               theta=theta,phi=phi,bres=bres'
    return
endif


if not keyword_set(nproc) then nproc = 4
if not keyword_set(lmax) then lmax = 2500l
if not keyword_set(nside) then nside=2048l
if not keyword_set(dipole_nside) then dipole_nside=512l 
if not keyword_set(firstb) then firstb = 0l
if not keyword_set(nband) then nband = 15l
if not keyword_set(aps) then aps = 0d0
if not keyword_set(nnbin) then nnbin = 6l
if not keyword_set(nbj) then nbj = 0l
if not keyword_set(unnorm) then unnorm = 1l
nell = lmax+1l


if n_params() eq 2 then begin
   mcres = res
   dres=res
   mres=res
   nproc=2
   unnorm=2
endif

nbins=0

nspots = n_elements(mcres[*,0,0])
nrad = n_elements(mcres[0,*,0])
nsim = n_elements(mcres[0,0,*])
print, 'get_variance: nsim, nrad',nsim, nrad

if not keyword_set(mask) then mask = fltarr(nspots,nrad)
print, '************** dipole fit mask, nspots, nrad, nsim = ',total(mask),nspots, nrad, nsim

nbins=nrad
bins=indgen(nrad)+1
dl = long(floor(bins(nbins-1)- bins(nbins-2)))

l=lindgen(lmax)
pp=pixwin(nside)
pp=pp(0:lmax,0)

spawn,'mkdir -p temp'

wunf, double(res),'temp/res.unf'
wunf, double(dres),'temp/dres.unf'
wunf, double(mres),'temp/mres.unf'
wunf, double(mcres), 'temp/mcres.unf'
wunf, double(pp), 'temp/pixwin.unf'
wunf, long(l), 'temp/ell.unf'
wunf, double(mask), 'temp/mask.unf'
wunf,double(aps),'temp/Aps.unf'
wunf,long([[unnorm],[nbj],[dl/2l],[nrad], [nbins], [nell],[dipole_nside],[nsim],[nspots],[nnbin]]),'temp/params.unf'

spawn, 'setup-intel-compilers'
if not keyword_set(load) then begin
    spawn,'pwd',dir
    str='cd '+dir+'; mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME $HEALPIX3/bin/get_variance_dipdir'
    print, '************************* spawning **********************'
    print, '** --', str, ' -- **'
    print, '*********************************************************'
    SPAWN,str
endif 
if keyword_set(load) then print, 'loading file ..'
bres=dblarr(nrad,nspots)
theta=dblarr(nrad,nsim+1)
phi=theta
dipamp=theta
runf, bres, 'temp/bres.unf'


ntasks=nproc

nsim_pp=lindgen(ntasks)
nsim_pp[*]=(nsim+1)/ntasks
if (((nsim+1) mod ntasks) ne 0) then nsim_pp(0:((nsim+1) mod ntasks)-1) = nsim_pp(0:((nsim+1) mod ntasks)-1)+1


for ii=0,nproc-1 do begin
   theta_temp = dblarr(nrad,nsim_pp(ii))
   phi_temp=theta_temp
    if unnorm ne 0 then dipamp_temp=theta_temp

   runf, theta_temp, 'temp/theta_'+strn(ii)+'.unf'
   runf, phi_temp, 'temp/phi_'+strn(ii)+'.unf'
   if unnorm ne 0 then   runf, dipamp_temp, 'temp/dipamp_'+strn(ii)+'.unf'


   if ii eq 0 then begin
      s0 = 0
      sn = nsim_pp(0)-1
   endif else begin
      s0 = round(total(nsim_pp(0:ii-1)))
      sn = round(total(nsim_pp(0:ii))-1)
   endelse
   ;print, ii, s0,sn

   theta[*,s0:sn] = theta_temp
   phi[*,s0:sn] = phi_temp
   if unnorm ne 0 then dipamp[*,s0:sn] = dipamp_temp
endfor



end
