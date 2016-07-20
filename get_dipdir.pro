pro get_dipdir,bins,res,dres,dnres,mcres,mres,mnres,dipole_nside=dipole_nside,lmax=lmax,$
               nside=nside,nband=nband,aps=aps,nproc=nproc,nnbin=nnbin,$
               theta=theta,phi=phi,bres=bres,nbj=nbj,load=load

if n_params() lt 6 then begin
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
if not keyword_set(nband) then nband = 15l
if not keyword_set(aps) then aps = 0d0
if not keyword_set(nnbin) then nnbin = 6l
if not keyword_set(nbj) then nbj = 0l
nell = lmax+1l

nbins = n_elements(bins)
dl = long(floor(bins(nbins-1)- bins(nbins-2)))

nspots = n_elements(mcres[0,*,0])
nsim = n_elements(mcres[0,0,*])

l=lindgen(lmax)
pp=pixwin(nside)
pp=pp(0:lmax,0)

spawn,'mkdir -p temp'

wunf, double(res),'temp/res.unf'
wunf, double(dres),'temp/dres.unf'
wunf, double(mres),'temp/mres.unf'
wunf, double(mnres),'temp/mnres.unf'
wunf, double(dnres),'temp/dnres.unf'
wunf, double(mcres), 'temp/mcres.unf'
wunf, double(pp), 'temp/pixwin.unf'
wunf, long(l), 'temp/ell.unf'
wunf, long(bins), 'temp/bins.unf'
wunf,double(aps),'temp/Aps.unf'
wunf,long([[nbj],[dl/2l],[nband], [nbins], [nell],[dipole_nside],[nsim],[nspots],[nnbin]]),'temp/params.unf'

spawn, 'setup-intel-compilers'
if not keyword_set(load) then begin
    spawn,'pwd',dir
    str='cd '+dir+'; mpirun -np '+strtrim(string(nproc),2)+' $HOME/software/Healpix_2.15/bin/get_dipdir'
    print, '************************* spawning **********************'
    print, '** --', str, ' -- **'
    print, '*********************************************************'
    SPAWN,str
endif 


bres=dblarr(nband,nspots)
theta=dblarr(nband,nsim+1)
phi=theta

runf, bres, 'temp/bres.unf'
runf, theta, 'temp/theta.unf'
runf, phi, 'temp/phi.unf'


end
