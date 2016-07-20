function multipeak_clean,nside,ipix,pixrad=pixrad,ipeaks=ipeaks,vpeaks=vpeaks
;; This function first sets a disk around each ipix pixel, and then
;; groups any ipix falling with in a radius of PIXRAD as one.
;;
;;Input:
;; NSIDE - HELAPIX nside
;; IPIX -  Healpix pixel numbers of peak locations
;; PIXRAD - the radius in degree within which peaks are merged
;;OUTPUT:
;; IPEAKS - IPIX with close peaks merged
;;INPUT OUTPUT:
;; VPEAKS - as input the values at IPIX locations. As output the
;;          values at IPEAKS
;;

  ;;create the logical indices of peaks. Disks will be labeled 
  ;;by these indices
  idpeaks=lindgen(n_elements(ipix))
  
  ;;define a map to contain disks
  map0 = lonarr(nside2npix(nside))
  map0[*]=-100

  ;;make disks around ipix with values being the logical index of ipix  
  ;;map0 values are now either index of peak or -100
  map0 = make_pix_disk(ipix,idpeaks,pixrad, nside,map=map0,/unique,bad_value=-100)

  ;;now we want to look at the values of the map at ipix locations
  ;;if they are unique that means all ipix peaks are well separated.
  ;;those of ipix values very close to each other will have similar 
  ;;value. Idpeaks contains repeated values for the later case
  idpeaks = map0[ipix]

  ;;Not get the uniq idpeaks
  idpeaks_uniq=unique(idpeaks,count,/sort)  

  ;;
  ipeaks=ipix[idpeaks_uniq]
  if keyword_set(vpeaks) then vpeaks = vpeaks[idpeaks_uniq]

  return, idpeaks_uniq

end
