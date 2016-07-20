pro bin_res,bins,res,bres,lmax=lmax,nband=nband,nsim_use=nsim_use,aps=aps,nnbin=nnbin,bbins=bbins,error=error,maxnbins=maxnbins


if n_params() lt 3 then begin
    print, 'USAGE: '
    print, 'bin_res,bins,res,bres,lmax=lmax,nside=nside,nband=nband,nsim_use=nsim_use'
    return
endif



if not keyword_set(lmax) then lmax = 2500l
if not keyword_set(nband) then nband = 15l
if not keyword_set(aps) then aps = 0d0
if not keyword_set(maxnbins) then maxnbins = 1500l/16l

nbin2bin=6l
if keyword_set(nnbin) then nbin2bin = nnbin

;;print, 'nbin2bin: ', nbin2bin

nell = lmax+1l

nbins = n_elements(bins)
dl = long(floor(bins(nbins-1)- bins(nbins-2)))

nspots = n_elements(res[0,*])


l=lindgen(lmax)
nsim = n_elements(res[0,0,*])

bres=dblarr(nband,nspots,nsim)
bbins = dblarr(nband)

;;--------------------------
for k=0,nsim-1 do begin
for i=0,nband-1 do begin
    for j=0,nspots-1 do begin
        b1=i*nbin2bin
        b2=(i+1)*nbin2bin-1l
        if nbin2bin eq 6l then begin
           if (i ge 2) then b2=b2+1
           if (i gt 2) then b1=b1+1
           if (i ge 6) then b2=b2+1
           if (i gt 6) then b1=b1+1
           if (i ge 10) then b2=b2+1
           if (i gt 10) then b1=b1+1
           if (i ge 14) then b2=b2+1
           if (i gt 14) then b1=b1+1
        endif

        nb = 1

        b22 = min([b2,maxnbins-1])

        if keyword_set(error) then nb = 1./sqrt(b22-b1+1.)
        ;;print, 'i, b1, b2, b22: ',i, b1,b2, b22


        bbins[i] = mean(bins[b1:b22])
        bres[i,j,k]=mean(res[b1:b22,j,k])*nb
    endfor
endfor
endfor



end
