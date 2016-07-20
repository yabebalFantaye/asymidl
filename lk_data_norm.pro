pro lk_data_norm, sig, dersig, Ak=Ak


sig0_new = sqrt(2d0)*sig 
omega=[1d0,2d0,!pi]

nj=n_elements(sig)
nf=3 ;;0-area, 1-length, 2-genus
Ak = dblarr(nj,nf)

inv2pi = 1d0/(2d0*!pi)

for k=0,nf-1 do begin
   omega_frac = omega[2]/(omega[2-k]*omega[k])

   Ak(*,k) = inv2pi^((k+1d0)/2d0)*omega_frac*(dersig/sig0_new)^k
endfor




end
