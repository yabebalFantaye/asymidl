pro get_dipdir_idl,bins,res,dres,dnres,mcres,mres,mnres,dipole_nside=dipole_nside,lmax=lmax,$
               nside=nside,nband=nband,aps=aps,nproc=nproc,nsim_use=nsim_use,$
               theta=theta,phi=phi,bres=bres,load=load

if n_params() lt 6 then begin
    print, 'USAGE: '
    print, 'get_dipdir_idl,bins,res,dres,dnres,mcres,mres,mnres,nside=nside,lmax=lmax,$'
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
nell = lmax+1l

nbins = n_elements(bins)
dl = long(floor(bins(nbins-1)- bins(nbins-2)))

nspots = n_elements(mcres[0,*,0])
nsim = n_elements(mcres[0,0,*])

l=lindgen(lmax)
pp=pixwin(nside)
pp=pp(0:lmax,0)

if not keyword_set(nsim_use) then nsim_use = nsim+1

bres=dblarr(nband,nspots)
theta=dblarr(nband,nsim_use+1)
phi=theta


; print, '------------------'
; print, 'res: ', double(res(0:10,0:1))
; print, '------------------'
; print, 'mres: ', mres(0:10,0:1)
; print, '------------------'
; print, 'dres: ', dres(0:10,0:1)
; print, '------------------'
; print, 'mnres: ', mnres(0:10,0:1)
; print, '------------------'
; print, 'dnres: ',dnres(0:10,0:1)
; print, '------------------'
; return


;;--------------------------
for sim=0,nsim_use do begin
    for i=0,nband-1 do begin
        for j=0,nspots-1 do begin
            b1=i*6l
            b2=(i+1)*6l-1l
            if (i ge 2) then b2=b2+1
            if (i gt 2) then b1=b1+1
            if (i ge 6) then b2=b2+1
            if (i gt 6) then b1=b1+1
            if (i ge 10) then b2=b2+1
            if (i gt 10) then b1=b1+1
            if (i ge 14) then b2=b2+1
            if (i gt 14) then b1=b1+1

            if (sim eq 0) then begin
                sres=mres(b1:b2,j)+aps/pp(bins(b1:b2))^2
                bres(i,j)=mean(res(b1:b2,j)-sres)
                dres0=dres(b1:b2,j)^2+(2d0*sres*mnres(b1:b2,j))/(2d0*(l(bins(b1:b2))+dl/2d0)+1d0)+dnres(b1:b2,j)^2/2d0
                ddres=sqrt(total(dres0)/double(b2-b1+1d0))
                bres(i,j)=bres(i,j)/ddres

            endif else begin
                sres=mres(b1:b2,j)
                bres(i,j)=mean(mcres(b1:b2,j,sim-1)-sres)
                dres0=dres(b1:b2,j)^2
                ddres=sqrt(total(dres0)/double(b2-b1+1d0))
                bres(i,j)=bres(i,j)/ddres
            endelse        
        endfor

        ;;make map of cl-distrib.
        map=reform(bres(i,*))
        ;;get dipole 
         ;;if(sim eq 0) then print, map(0:100)
        map2alm,map,alm,lmax=1
        alm2map,alm,dp,nside=dipole_nside
        ii=where(dp eq max(dp))
        ii=ii[0]
        pix2ang_ring,512l,ii,th,ph

        theta(i,sim)=th
        phi(i,sim)=ph
        ;;if(sim eq 0) then print, 'ii,th, ph: ',ii,th,ph

    endfor
endfor



end
