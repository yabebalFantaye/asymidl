pro read_asymfast_res,res,nsim=nsim,split=split,nbins=nbins,froot=froot,nspots=nspots

   ;; split=10l
   ;; nsim=500l
   ;; nbins=107

  if not (keyword_set(split) and keyword_set(nbins) and keyword_set(froot) and keyword_set(nsim) and keyword_set(nspots) ) then begin
     print, 'you must set: split, nbins, froot, nsim, nspots - keywords'
     return
  endif

;nsim=50
   nsim_use=nsim

   res0=dblarr(nbins,nspots,split)
   res=dblarr(nbins,nspots,nsim)
   nblock=nsim/split


   for i=0l,nblock-1l do begin
      simno=trim(string(i),2)
      if (strlen(simno) eq 1) then simno='000'+simno
      if (strlen(simno) eq 2) then simno='00'+simno
      if (strlen(simno) eq 3) then simno='0'+simno
      file=froot+simno+'.unf'

;print,file
      runf,res0,file
      res(*,*,i*split:(i+1l)*split-1l)=res0
   endfor

end
