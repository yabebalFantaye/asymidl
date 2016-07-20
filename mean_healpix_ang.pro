function mean_healpix_ang, th, phi,mth=mth,mph=mph,astro=astro,mvec=mvec,nband=nband,median=median
;;
;;This function makes disks on a healpix map of nside which subtends a
;;radious of rad and sets the value of it to be val
;;
;;

  if n_params() lt 2 then begin
     print, 'USAGE: thph=mean_healpix_ang(th, ph,mth=mth,mph=mph,astro=astro)'
     on_error, 0
  endif

  nn = n_elements(th[*,0])
  nsim = n_elements(th[0,*])



  thph=dblarr(2,nsim)
  mvec=dblarr(nsim, 3)
  vec_arr=dblarr(nn, 3)

  nband_use=1
  get_bin_dipdir,bincenter,binmin,binmax,bvec=bvec

  if keyword_set(nband) then begin
     nband_use=nband
     thph=dblarr(nband,nsim,2)
  endif

  if nband_use eq 1 then begin
     bvec[0,0]=0
     bvec[0,1]=nn-1
  endif

  for j=0,nsim-1 do begin

     for k=0,nband_use-1 do begin
        vec_mean=dblarr(3)
        if j eq 0 then print, 'averaging binmin to binmax: ',bvec[k,0],bvec[k,1]

        count=0
        vec_arr[*,*]=0.
        for i=bvec[k,0],min([bvec[k,1],nn-1]) do begin
           
                                ;print, th(i), phi(i)
           if keyword_set(astro) then begin
              ANG2VEC,  (!pi/2d0-th(i,j))*!rtod, phi(i,j)*!rtod, vec,/astro
           endif else begin
              ANG2VEC,  th(i,j), phi(i,j), vec
           endelse
           
           vec_arr[count,*]=vec
           count=count+1
        endfor
        count=count-1

        vec_mean = mean(vec_arr[0:count,*],dimension=1)
        if keyword_set(median) then begin
           if j eq 0 then print, 'using median to average ..'
           vec_mean = median(vec_arr[0:count,*],dimension=1,/even,/double)
        endif
        vec=vec_mean
        vec_mean= vec/sqrt(vec[0]^2+vec[1]^2+vec[2]^2)
        
        ;;if j eq 0 then print,'vec_mean', vec_mean
        ;;if j eq 1 then help, mvec
        mvec[j,*]=reform(vec_mean)
        
        VEC2ANG,  reform(vec_mean),mth,mph ;,astro=astro
        
        if not keyword_set(nband) then begin 
           thph[0,j] = mth
           thph[1,j] = mph
        endif else begin
           thph[k,j,0] = mth
           thph[k,j,1] = mph
        endelse

     endfor
  endfor


  return, thph

end
