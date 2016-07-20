function read_mf_peaks, froot, isim,ifwhm, avg=avg,rms=rms,minmax=minmax,in_nside=in_nside,$
                        nside=nside,nomap=nomap,pixrad=pixrad,iminima=iminima,vminima=vminima,$
                        imaxima=imaxima,vmaxima=vmaxima,isaddle=isaddle,vsaddle=vsaddle,$
                        nested=nested,nopmap=nopmap,hbins=hbins, inring=inring,$
                        vhot=vhot, vcold=vcold

  if n_params() lt 3 then begin
     print, 'usage: read_mf_peaks, froot, isim,ifwhm, avg=avg,rms=rms,minmax=minmax,peakmaps=peakmaps,in_nside=in_nside,nside=nside,nomap=nomap'
     return,0
  endif

  unit=99
  status = FSTAT(unit)
  if (status.open ne 0) then begin
     WHILE(status.open ne 0) do begin
        unit=unit-1
        status = FSTAT(unit)
     ENDWHILE
  endif

  order_in='ring'
  if keyword_set(nested) then order_in='nest'


  minval=0d0
  maxval=0d0
  avgm0=0d0
  avgm1=0d0
  rmsm0=0d0
  rmsm1=0d0

  npeaks=0L

  print, 'isim, ifwhm',isim, ifwhm  
  openr,unit,froot+'.unf_ipix_val_map'+strn(isim)+'_j'+strn(ifwhm),/f77

  ;;read min, max, avg and rms for map and its 1st derivative
  readu,unit,minval, maxval
  readu,unit,avgm0, rmsm0
  readu,unit,avgm1, rmsm1

  avg=[avgm0,avgm1]
  rms=[rmsm0,rmsm1]
  minmax=[minval, maxval]


  peakmaps=0

  if not keyword_set(nomap) then begin

     peaks_nside=0L
     if not keyword_set(in_nside) then  readu,unit,peaks_nside
     if keyword_set(in_nside) then  peaks_nside=in_nside

     ;;read MINIMA length and read pixel list and values
     readu,unit,npeaks
     iminima=lonarr(npeaks)
     vminima=dblarr(npeaks)
     readu,unit,iminima
     readu,unit,vminima


     ;;read MAXIMA length and read pixel list and values
     readu,unit,npeaks
     imaxima=lonarr(npeaks)
     vmaxima=dblarr(npeaks)
     readu,unit,imaxima
     readu,unit,vmaxima



     ;;read SADDLE length and read pixel list and values
     readu,unit,npeaks
     isaddle=lonarr(npeaks)
     vsaddle=dblarr(npeaks)
     readu,unit,isaddle
     readu,unit,vsaddle


     ntotmaxima=n_elements(imaxima)
     
     if not keyword_set(nopmap) then begin

        if not keyword_set(pixrad) then begin
           peakmaps=fltarr(nside2npix(peaks_nside),3)
           peakmaps[*,*]=!healpix.bad_value
           
           peakmaps[iminima,0]=vminima
           peakmaps[imaxima,1]=vmaxima
           peakmaps[isaddle,2]=vsaddle
        endif
        
        if keyword_set(pixrad) then begin
           bad_value=!healpix.bad_value
           peakmaps=fltarr(nside2npix(peaks_nside),4)
           peakmaps[*,*]=bad_value
           peakmaps[*,0] = make_pix_disk(iminima,vminima,pixrad,peaks_nside,bad_value=bad_value)
           ;;peakmaps[imaxima,1] = make_pix_disk(imaxima,vmaxima,pixrad,peaks_nside,bad_value=bad_value,/square)
           ;;peakmaps[isaddle,2] = make_pix_disk(isaddle,vsaddle,pixrad,peaks_nside,bad_value=bad_value,/triangle)
           ;;peakmaps[*,3] = total(peakmaps, 2)
        endif

        if keyword_set(nside) then begin
           print, 'ud_grading peak map to nside ',nside
           peakmaps_old=peakmaps
           peakmaps=0
           ud_grade,peakmaps_old,peakmaps, nside_out=nside,order_in=order_in,order_out='ring'
        endif
        
     endif
     

     ;;transform the peak pixels into ring order
     if keyword_set(inring) then begin
        nest2ring, peaks_nside, isaddle, ipring
        isaddle=ipring

        nest2ring, peaks_nside, imaxima, ipring
        imaxima=ipring

        nest2ring, peaks_nside, iminima, ipring
        iminima=ipring
     endif


     ;;return maximas above or below a certain threshold
     if keyword_set(vhot) or keyword_set(vcold) then begin        
        if not keyword_set(vcold) then vcold=max(vmaxima)  ;;maximum threshod
        if not keyword_set(vhot) then vhot=min(vmaxima)    ;;minimum threshold
        
        id_hot=where(vmaxima ge vhot and vmaxima le vcold ,count)
        if count ne 0 then begin
           vmaxima=vmaxima[id_hot]
           imaxima=imaxima[id_hot]

           print, 'read_mf_peaks: ntotmaxima=, nmaxima= above sigma=',ntotmaxima,n_elements(imaxima),vhot
        endif else begin
           print, 'WARNING: There are no maximas above sigma=',vhot
           vmaxima=[] & imaxima=[]
        endelse
     endif


     ;;clean peaks that fall within pixrad degree radius
     if keyword_set(pixrad) then begin
        print, 'read_mf_peaks: peaks within pixrad= will be merged as one.',pixrad
        idpeaks = multipeak_clean(peaks_nside,imaxima,pixrad=pixrad,ipeaks=imaxima,vpeaks=vmaxima)
     endif
           
     sort_indx=sort(vmaxima)
     vmaxima=vmaxima[sort_indx]
     imaxima=imaxima[sort_indx]


     if keyword_set(hbins) then begin
        vminima=(histogram(vminima,nbins=n_elements(hbins),min=min(hbins),max=max(hbins)))/n_elements(vminima)
        vmaxima=(histogram(vmaxima,nbins=n_elements(hbins),min=min(hbins),max=max(hbins)))/double(n_elements(vmaxima))
        vsaddle=(histogram(vsaddle,nbins=n_elements(hbins),min=min(hbins),max=max(hbins)))/n_elements(vsaddle)
     endif

     if keyword_set(nopmap) then peakmaps=imaxima

  endif

  close,unit


  return,peakmaps

end
