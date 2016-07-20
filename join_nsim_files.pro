pro join_nsim_files,fname,ntasks,ndim,mat,fnamex=fnamex,matx=matx,$
                    fout=fout,fmatx=fmatx,prec=prec

  if n_params() lt 4 then begin
     print, 'USAGE: join_nsim_files,fname,ntasks,ndim,mat,fnamex=fnamex,matx=matx,prec=prec'
     return
  endif

if not keyword_set(prec) then prec=4

;;get dimension of the big matrix
  case ndim of
     2: begin 
        nsim = n_elements(mat[0,*])
        n0 = n_elements(mat[*,0])
     end
     3: begin
        nsim = n_elements(mat[0,0,*])
        n0 = n_elements(mat[*,0,0])
        n1 = n_elements(mat[0,*,0])
     end
     4: begin 
        nsim = n_elements(mat[0,0,0,*])
        n0 = n_elements(mat[*,0,0,0])
        n1 = n_elements(mat[0,*,0,0])
        n1 = n_elements(mat[0,0,*,0])
     end
     else: begin
        print, '>4 dimension not implemented, change asymidl/join_nsim_files.pro'
        on_error, 1
     end
  endcase


;;get how the big matrix is divided at its last index
  nsim_pp=lindgen(ntasks)
  nsim_pp[*]=nsim/ntasks
  if (((nsim) mod ntasks) ne 0) then nsim_pp(1:((nsim) mod ntasks)-1) = nsim_pp(1:((nsim) mod ntasks)-1)+1


;;load the files at each task and put it in a big matrix
  for ii=0,ntasks-1 do begin

     if prec eq 4 then begin
        case  ndim of
           2: mat_temp = fltarr(n0,nsim_pp(ii))
           3: mat_temp = fltarr(n0,n1,nsim_pp(ii))
           4: mat_temp = fltarr(n0,n1,n2,nsim_pp(ii))
        endcase
     endif else begin
        case  ndim of
           2: mat_temp = dblarr(n0,nsim_pp(ii))
           3: mat_temp = dblarr(n0,n1,nsim_pp(ii))
           4: mat_temp = dblarr(n0,n1,n2,nsim_pp(ii))
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
        2: begin
           mat[*,s0:sn] = mat_temp

           if keyword_set(fnamex) then begin
              matx_temp = mat_temp
              runf, matx_temp, fnamex+'_'+strn(ii)
              matx_temp[*,s0:sn] = matx_temp
           endif
        end
        3: begin
           mat[*,*,s0:sn] = mat_temp

           if keyword_set(fnamex) then begin
              matx_temp = mat_temp
              runf, matx_temp, fnamex+'_'+strn(ii)
              matx_temp[*,*,s0:sn] = matx_temp
           endif
        end
        4: begin
           mat[*,*,*,s0:sn] = mat_temp

           if keyword_set(fnamex) then begin
              matx_temp = mat_temp
              runf, matx_temp, fnamex+'_'+strn(ii)
              matx_temp[*,*,*,s0:sn] = matx_temp
           endif
        end

     endcase
  endfor

  if not keyword_set(fout) then fout=fname
  print, 'saving whole matrix file to: ',fout
  wunf, mat,fout


  if keyword_set(fnamex) then begin
     if keyword_set(fmatx) then fmatx=fnamex
     print, 'saving whole matrix file to: ',fmatx
     wunf, matx,fmatx
  endif


end
