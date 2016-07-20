function galactic_mask,nside,file=file

npix=nside2npix(nside)
mask = fltarr(npix)+1.0

listpix=lindgen(npix)
pix2ang_ring, nside,listpix,th, phi

iii = where(th gt (!pi/2-10*!dtor) and th lt (!pi/2+10*!dtor))
mask[iii]=0.0

if keyword_set(file) then write_fits_map,file,mask,/ring

return,mask

end
