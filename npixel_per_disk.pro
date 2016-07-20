function npixel_per_disk,nside,nspots,rval

if n_params() lt 3 then begin
   print, 'usage: e.g. n_pixels=npixel_per_disk(nside,nspots,rval)'
   ;;return
   on_error,2
endif

nrad = n_elements(rval)

print, nspots, nrad

n_pixels = lindgen(nspots,nrad)

for j=0,nspots-1 do begin           
   for i=0,nrad-1 do begin              
      
      pix2vec_ring,nside,j,vec
      query_disc,nside,vec,rval[i]*!dtor,listpix
      n_pixels[j,i] = n_elements(listpix)
   endfor
endfor

fname='/mn/owl1/d3/yabebalf/planck/common_files/n_pixels_nside'+strn(nside)+'_nspot'+strn(nspots)+'_nrad'+strn(nrad)+'.unf'
print, 'saving N_pixels to: '
print, fname
wunf, n_pixels,fname

return, n_pixels


end
