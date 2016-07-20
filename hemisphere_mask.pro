pro hemisphere_mask,theta, phi,nside,maskn,masks,deg=deg,shown=shown,shows=shows,show=show

if n_params() lt 3 then begin
   print, 'Usage: hemisphere_mask,ang,nside,mask_north,mask_south'
endif


maskn = fltarr(nside2npix(nside))
maskn[*] = 0.0
masks = maskn

if keyword_set(deg) then begin
   theta = theta*!dtor
   phi = phi*!dtor
endif

ANG2VEC,  theta, phi, vec
rad = !pi/2d0
query_disc,nside,vec,rad,listpix
maskn(listpix) = 1.0
masks = 1.0 - maskn


if keyword_set(show) then begin
   shown=1
   shows=1
endif

if keyword_set(shown) then    mollview, maskn
if keyword_set(shows) then    mollview, masks


end
