pro get_listpix_disk,nside_map,nside_spots,file=file,nproc=nproc,unit=unit,overwrite=overwrite,radvec=radvec,isnest=isnest
;; This code produces nlist and listpix in *.unf file

if n_params() lt 2 then begin
    print, 'USAGE: '
    print, 'get_listpix_disk,nside_map, nside_spots,file=file,nproc=nproc,unit=unit,radvec=radvec'
    return
endif


if not keyword_set(nproc) then nproc=150
if not keyword_set(radvec) then radvec=[16,18,20,22,30,50,70,90]  ;1,2,4,6,8,10,12,14

ordering='ring'
excutable = 'get_listpix_ring_disk'
if keyword_set(isnest) then begin
   ordering='nest'
   excutable = 'get_listpix_nest_disk'
endif

nproc=min([nproc,nside2npix(nside_spots)])

rdir='/mn/owl1/d3/yabebalf/planck/common_files/listpix_'+ordering+'/'
spawn,'mkdir -p '+rdir

nrad=n_elements(radvec)
nspots=nside2npix(nside_spots)


wunf, long(nside_map), rdir+'nside_map.unf'
wunf, long(nside_spots), rdir+'nside_spots.unf'


for irad=0,nrad-1 do begin

   rad=radvec[irad]

   file='/mn/owl1/d3/yabebalf/planck/common_files/listpix_'+ordering+'_nin'+strn(nside_map)+'_nout'+strn(nside_spots)+'_rad'+strn(floor(rad))+'.unf'
   
   
   print, 'idl: calling '+excutable
   print, 'idl: nside_map, nside_spots, rad, nproc: ',nside_map, nside_spots, rad, nproc

   wunf,double(rad),rdir+'rad_listpix.unf' ;;in degrees

   if keyword_set(overwrite) or file_test(file,/regular) eq 0 then begin
      spawn,'pwd',dir
      str='cd '+dir+'; time mpirun -np '+strtrim(string(nproc),2)+' -hostfile $HOME/mpihosts --prefix $MPI_HOME $HOME/software/Healpix_2.15/bin/'+excutable
      print, '************************* spawning **********************'
      print, '** --', str, ' -- **'
      print, '*********************************************************'
      SPAWN,str
   endif

endfor

if keyword_set(unit) then begin
   print, 'listpix file: '+file+' is open. Close unit='+strn(unit)+' at the end.'
   openr, unit, file, /f77
endif

end
