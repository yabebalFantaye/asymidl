pro join_spotpar_files,fname,ntasks,ndim,mat,fnamex=fnamex,matx=matx,$
                    fout=fout,fmatx=fmatx,prec=prec

  if n_params() lt 4 then begin
     print, 'USAGE: join_nsim_files,fname,ntasks,ndim,mat,fnamex=fnamex,matx=matx,prec=prec'
     return
  endif

if not keyword_set(prec) then prec=4

;;get dimension of the big matrix
  case ndim of

     3: begin
        nsim = n_elements(mat[0,0,*])
        n0 = n_elements(mat[*,0,0])
        nspots = n_elements(mat[0,*,0])
     end

     else: begin
        print, 'ndim ~= 4 dimension not implemented, change asymidl/join_spotpar_files.pro'
        on_error, 1
     end
  endcase

print, 'nsim,nspots,nbins',nsim, nspots, n0

;;get how the big matrix is divided at its last index
  nsim_pp=lindgen(ntasks)
  nsim_pp[*]=nspots/ntasks
  if (((nspots) mod ntasks) ne 0) then nsim_pp(1:((nspots) mod ntasks)-1) = nsim_pp(1:((nspots) mod ntasks)-1)+1


;;load the files at each task and put it in a big matrix
  for ii=0,ntasks-1 do begin

     if prec eq 4 then begin
        case  ndim of
           3: mat_temp = fltarr(n0,nsim_pp(ii),nsim)
        endcase
     endif else begin
        case  ndim of
           3: mat_temp = dblarr(n0,nsim_pp(ii),nsim)
        endcase
     endelse

     runf, mat_temp, fname+'_'+strn(ii)


     if ii eq 0 then begin
        s0 = 0
        sn = nsim_pp(0)-1
     endif else begin
        s0 = round(total(nsim_pp(0:ii-1)))
        sn = round(total(nsim_pp(0:ii))-1)
     endelse
                                ;print, ii, s0,sn

     case ndim of
        3: begin

           count=0
           for jj=s0,sn do begin
              mat[*,jj,*] = mat_temp[*,count,*]
              count = count+1
           endfor

        end

     endcase
  endfor


  if not keyword_set(fout) then fout=fname
  print, 'saving whole matrix file to: ',fout
  wunf, mat,fout




end
