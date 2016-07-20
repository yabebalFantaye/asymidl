function  read_listpix,fname,nlistvec=nlistvec,ntot=ntot,nrad=nrad,map=map,nside=nside,nrout=nrout,$
                       tmap=tmap,ipix=ipix,listtmap=listtmap,ibrad=ibrad,ifrad=ifrad,nsspot=nsspot

if not keyword_set(ipix) then ipix=8
if not keyword_set(nsspot) then nsspot=16

if n_params() lt 1 then begin
   fname='../../common_files/listpix_ring_nin2048_nout16_nrad16/listpix'
   print, 'using fname,ipix='+fname,ipix
endif
fname=fname+'_'+strn(ipix)+'.unf'

print,'reading: '+fname

openr,12,fname,/f77
ntot=0l
nrad=0l
readu,12, nrad, ntot

nlistvec=lindgen(nrad)
readu,12,nlistvec

print,'ntot, total(nlistve):',ntot,total(nlistvec)
ntot=total(nlistvec)
listpix=lindgen(ntot)
readu,12,listpix

close,12


print, 'nrad, ntot:',nrad,ntot
print, 'nlistvec=',nlistvec
help, listpix

;rval=[50,70,90]
;rval=[10,12,14,16,18,20,22,30,50,70,90]
rval=[1,2,4,6,8,10,12,14,16,18,20,22,30,50,70,90]
;rval=[70,90]

nlistvec1 = floor(total(nlistvec,/cumulative))
print, 'cumulative nlistvec:',nlistvec1

nlistvec2=nlistvec*0


if keyword_set(nrout) then nrad=nrout
if not keyword_set(nside) then nside=16



map=fltarr(nside2npix(nside))
tmap=map
nlistvec2=0
nlistvec3=0

if not keyword_set(ibrad) then begin
   irad=0   

   outlist=listpix[0:nlistvec[0]-1]
   map[ outlist ]=10.
   nlistvec3 = nlistvec1[0]
   ;if (nlistvec3 lt 100) then print, 'mod irad0:',outlist(sort(outlist))

   tmap=make_pix_disk(ipix,10,rval[0],nside,nsspot=nsspot,bad_value=0)
   nlistvec2=n_elements(where(tmap eq 10))
   listtmap=where(tmap eq 10)
   ;if (nlistvec2 lt 100) then print,'orig irad0:',listtmap(sort(listtmap))

   print, 'ipix='+strn(ipix)+' irad='+strn(irad+1)+', total(diff)',total((listtmap(sort(listtmap))-outlist(sort(outlist))))

   ibrad = 1
endif

if not keyword_set(ifrad) then ifrad=nrad-1

for irad=ibrad,ifrad do begin

   np=nlistvec1[irad]
   nlistvec3 = [nlistvec3,np]
   outlist=listpix[0:np-1]
   map[outlist]=map[outlist]+irad+20

   zmap=make_pix_disk(ipix,irad+20,rval[irad],nside,nsspot=nsspot,bad_value=0)
   listtmap=where(zmap eq irad+20)
   nlistvec2=[nlistvec2,n_elements(where(zmap eq irad+20))]
   tmap=tmap + zmap
   print, 'ipix='+strn(ipix)+' irad='+strn(irad+1)+', total(diff)',total((listtmap(sort(listtmap))-outlist(sort(outlist))))

endfor

print,'modified nlistvec',nlistvec3
print,'correct nlistvec',nlistvec2

return,outlist
end
