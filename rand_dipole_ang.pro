pro rand_dipole_ang,nbins,nsim,ang_seps,outdir=outdir

nn=nbins
nsm=nsim

cth=acos(randomu(3745654,nsm,nn)*2d0-1d0)
cph=randomu(5111654,nsm,nn)*2d0*!pi 

if keyword_set(outdir) then begin
    outdir = outdir+'/'+'bins'+strn(nn)+'_nsim'+strn(nsim)+'/'
    figdir = outdir+'/'+'figures/'
    spawn,'mkdir -p '+outdir
    spawn,'mkdir -p '+figdir
endif

cor=dblarr(nn,nn)
dcor=cor
allcor=dblarr(nn,nn,nsm)

print, 'compute angle separation and dispersion of dipole directions ...'
for i=0,nn-1 do begin
    for j=i,nn-1 do begin
        for k=0l,nsm-1l do begin
            ANG2VEC, cth(k,i), cph(k,i), vec1
            ANG2VEC, cth(k,j), cph(k,j), vec2
            if (i ne j) then cor(i,j)=cor(i,j)+acos(FLOAT(total(vec1*vec2)))*180d0/!pi
            if (i ne j) then dcor(i,j)=dcor(i,j)+(acos(FLOAT(total(vec1*vec2)))*180d0/!pi)^2
            if (i ne j) then allcor(i,j,k)=acos(FLOAT(total(vec1*vec2)))*180d0/!pi
        endfor
    endfor
endfor
print, 'for loops for sep angle and disp done'

cor=cor/double(nsm)
dcor=dcor/double(nsm) ;; error^2 = sum(err^2/n^2)
dcor=sqrt(dcor-cor^2)

if keyword_set(outdir) then begin
    print, 'write sep angles and dispersions to a file'
    wunf, cor ,outdir+'dipole_dir_diff_cor.unf'
    wunf, dcor ,outdir+'dipole_dir_diff_dcor.unf'
    wunf, allcor ,outdir+'dipole_dir_diff_allcor.unf'
endif

nnx=nn
nn0=0

v=dblarr(nn,nn)

;;set one to lower diagonal matrix and 0 to diagonal and upper
;;diagonal parts
count = 0l
for i=nn0,nnx-2 do begin
    for j=i+1,nnx-1 do begin
        v(i,j)=1.
        count = count+1l
    endfor
endfor


print, 'compute the mean of separation angles per simulation '
mm=dblarr(nsm)
for i=0l,nsm-1l do begin
    mm(i)=total((allcor(*,*,i))*v)/total(v)
endfor

;;put separation angles per simulation in a column and add the mean at
;;the end
ang_seps = dblarr(count+1,nsm)
ij_inds=intarr(count+1,2)

for k=0,nsm-1l do begin
    count=0l
    for i=nn0,nnx-2 do begin
        for j=i+1,nnx-1 do begin
            ang_seps(count, k)=allcor[i,j,k]
            ij_inds[count, 0:1] = [i,j]
            count = count+1l
        endfor
    endfor
    ang_seps(count, k)=mm(k)  ;;save mean at the end
    ij_inds[count, 0:1] = [100,100]
endfor

if keyword_set(outdir) then begin

    wunf, ang_seps ,outdir+'dipole_dir_ang_seps_ncol'+strn(count)+'.unf'


    for i=0,count do begin
        
        hist = Histogram(ang_seps[i,*], min=MIN(ang_seps[i,*]), max=Max(ang_seps[i,*]))
        xbins = FINDGEN(N_ELEMENTS(hist)) + MIN(ang_seps[i,*]) 
        
        ps_start, file=figdir+'hist_dipole_dir_angsep_binij_'+strn(i)+'.ps' 
        
        str_bin =  'bin-'+strn(ij_inds[i,0])+' to bin-'+strn(ij_inds[i,1])+' dipole dtheta'
        
        if i eq count then str_bin = ' mean of all bins dipole dtheta'
        
        Plot, xbins, hist, YRANGE = [MIN(hist), MAX(hist)+1], PSYM = 10, $ ;
          title='Uniform generated dipole directions: '+str_bin, $
          YTITLE = 'Density per Bin', charsize=1, xtitle=textoidl('\Delta \theta (deg)')
        
        oplot,ang_seps[i,nsm-1]+0.*indgen(max(hist)+1), indgen(max(hist)+1),thick=3
        
        ps_end
    endfor

endif

end
