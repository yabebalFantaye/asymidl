pro sims_mf_peaks,froot, nsim,nj, max_vec=max_vec,rms_vec=rms_vec,peakmaps=peakmaps,in_nside=in_nside,nside=nside,bsim=bsim,nomap=nomap

isim0=0 ;;initial 
if keyword_set(bsim) then isim0=bsim

max_vec = dblarr(nsim,nj+1,2)
rms_vec = dblarr(nsim,nj+1,2)

;if not keyword_set(nside) then nside=128


for isim=isim0,nsim-1 do begin
   if isim mod 10 eq 0 then print, 'reading sim ',isim
   for ifwhm=0,nj do begin
      peakmaps=read_mf_peaks(froot, isim, ifwhm, avg=avg,rms=rms,minmax=minmax,in_nside=in_nside,nside=nside,nomap=nomap)
      max_vec[isim,ifwhm,*]=minmax
      rms_vec[isim,ifwhm,*]=rms
   endfor
endfor

print, 'simulation mean max values for j=9...'+strn(9+nj), mean(max_vec[*,*,1],dimension=1)

end
