pro write_asymfast_spot_res,res,nsim=nsim,split=split,nbins=nbins,froot=froot,nspots=nspots,ispot=ispot,fkey=fkey

 

if not keyword_set(nsim) then nsim=7300
if not keyword_set(split) then split=100
if not keyword_set(nbins) then nbins=62
if not keyword_set(froot) then froot="/mn/owl1/d3/yabebalf/wmap/9yr/allyrs/w9_output/asymfast/ap2_nside16/r90/res_asymfast_hemis_avgvw_w9kq85_apod2_lmax1200_62bins_nsim1000_"
if not keyword_set(foutroot) then foutroot="/mn/owl1/d3/yabebalf/wmap/parameter_estimation/importance/data/res_spectra/res_vw_asymfast_hemis_lmax1200_62bins_nsim100_"
if not keyword_set(ispot) then ispot=[2054,998]
if not keyword_set(fkey) then fkey=['asym_','antiasym_']
if not keyword_set(nspots) then nspots=3072

   ;; split=10l
   ;; nsim=500l
   ;; nbins=107

  if not (keyword_set(split) and keyword_set(nbins) and keyword_set(froot) and keyword_set(nsim) and keyword_set(nspots) and keyword_set(ispot) ) then begin
     print, 'you must set: split, nbins, froot, nsim, nspots, ispot - keywords'
     return
  endif
  if not keyword_set(fkey) then begin
     print, 'warning fkey not set, assuming fkey=asymdir_'
     fkey='asymdir_'
  endif

;nsim=50
   nsim_use=nsim


   res0=dblarr(nbins,nspots,split)
   nblock=nsim/split

for ijk=0,n_elements(ispot)-1 do begin

   key = fkey[ijk]
   ind_spot=ispot[ijk]

   for i=0l,nblock-1l do begin
      simno=trim(string(i),2)
      if (strlen(simno) eq 1) then simno='000'+simno
      if (strlen(simno) eq 2) then simno='00'+simno
      if (strlen(simno) eq 3) then simno='0'+simno

      file=froot+simno+'.unf'

      fileout=foutroot+key+strn(i)+'.unf'

      runf,res0,file
      res=res0[*,ind_spot,*]

      print, 'writing: '+fileout
      wunf,reform(res),fileout
   endfor

endfor
   help, reform(res)

end
