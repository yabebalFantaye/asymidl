pro get_bin_dipdir,bincenter,binmin,binmax,show=show,bvec=bvec,bminmax=bminmax,nnbin=nnbin


  if not keyword_set(nnbin) then  begin
     nnbin=6l
     iimax = 25l
     dl_band = 100l
  endif else begin
     if nnbin eq 6l then begin
        iimax = 25l
        dl_band = 100l
     endif else begin
        iimax = 150l/nnbin
        dl_band = nnbin*16l
     endelse
endelse

bincenter = 0*indgen(iimax)
binmin=bincenter
binmax=binmin
bvec=intarr(iimax,2)
ell = lindgen(2500)+2

for ii=1,iimax do begin
   i = ii-1

    b1=i*nnbin
    b2=(i+1)*nnbin-1
    if nnbin eq 6l then begin
       if (i ge 2) then b2=b2+1
       if (i gt 2) then b1=b1+1
       if (i ge 6) then b2=b2+1
       if (i gt 6) then b1=b1+1
       if (i ge 10) then b2=b2+1
       if (i gt 10) then b1=b1+1
       if (i ge 14) then b2=b2+1
       if (i gt 14) then b1=b1+1
    endif
    bvec[i,0] = b1
    bvec[i,1] = b2
    binmin[i] = ell[b1*16l]
    binmax[i] = ell[b2*16l]-1
    bincenter[i] = ceil(mean(ell[(b1*16l):(b2*16l)]))

;    binmin[i] = floor(min(2+b1*16l+indgen(b2*16l-b1*16l+1)))
;    bincenter[i] = ceil(mean(2+b1*16l+indgen(b2*16l-b1*16l+1)))
;    binmax[i]=ceil(max(2+b1*16l+indgen(b2*16l-b1*16l+1)))
 endfor


bminmax = intarr(iimax,2)
for ii=1,iimax do begin
   
   lmi = (ii-1)*dl_band
   lma = ii*dl_band
   if ii eq 1 then lmi=2
   
   i = ii-1
   bminmax[i,0] = lmi
   bminmax[i,1] = lma
endfor




if keyword_set(show) then begin
    print, '-------------------------'
    print, 'bincenter: ',bincenter[0:i]
    print, '-------------------------'
    print, 'binmin: ',binmin[0:i]
    print, '-------------------------'
    print, 'binmax: ',binmax[0:i]
    print, '-------------------------'
endif

end
