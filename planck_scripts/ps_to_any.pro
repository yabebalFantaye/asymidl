pro ps_to_any,feps=feps,png=png,cl=cl,large=large,medium=medium,latex=latex,charsize=charsize,$
                  xsize=xsize,ysize=ysize,axisx=axisx,axisy=axisy,mollview=mollview,fadd=fadd,nops=nops,dops=dops


    if keyword_set(feps) then begin
       if not file_test(feps,/regular) then message,feps+' not produced!'
    endif

   SZ = 88d
   if keyword_set(medium) then SZ=120d
   if keyword_set(large) then SZ=180d
   if not keyword_set(xsize) then xsize=SZ/25.4d
   if not keyword_set(ysize) then  ysize=SZ/10d/1.5d 

   if not keyword_set(charsize) then  charsize=1d
   if not keyword_set(dops) then  nops=1
    ;
    ;stop

   if keyword_set(fadd) then begin
      fadd = '_latexify' 
   endif else fadd = ''

   if keyword_set(feps) then  decompose,feps,d,path,name,ext,v
    ;
   if keyword_set(feps) and keyword_set(latex) then begin
       ;; height and width in cm

       fout = path+name+fadd+ext
       fout_pdf = path+name+fadd+'.pdf'
       fout_png = path+name+fadd+'.png'

       axis0=''
       axis1=''

       if keyword_set(cl) then begin
          axis0 = ['yttl','xttl','yttl2']
          axis1 = ['$\ell (\ell+1)C_\ell/2\pi$\quad $\left[\mu{\rm K}^2\right]$','$\ell$','$\Delta C_\ell/C_\ell [\%]$']
       endif

       if keyword_set(mollview) then begin
          axis0 = ['dipole','Central Multipole','-120','-60','0','60','120','-45','45','50','1450']
          axis1 = ['dipole','Central Multipole','240$^\circ$','300$^\circ$','0$^\circ$','60$^\circ$','120$^\circ$','-45$^\circ$','45$^\circ$','50','1450']          
       endif

       if keyword_set(axisx) then axis0=axisx
       if keyword_set(axisy) then axis1=axisy


       scale = DBLARR(N_ELEMENTS(axis0)) + charsize
       ls_latexify, feps, axis0,axis1 , scale*8d/11d, outname=fout, height=ysize*2.54d, width=xsize*2.54d, FDIR='' ;, /full

       ;IF Dname EQ 'X' THEN begin
       print, 'ps2pdf '+fout 
       SPAWN, 'cd '+path+' && ps2pdf '+fout 
       spawn,' pdfcrop '+fout_pdf+' '+fout_pdf ;	If in linux, convert the postscript to a pdf file straight away.
       ;endif
    endif

    ;;if  keyword_set(feps) and keyword_set(png) then spawn,'convert -alpha off -density 300 '+fout_pdf+' -flatten '+fout_png
    if  keyword_set(feps) and keyword_set(png) then begin
       spawn,'convert -alpha off -density 300 '+fout_pdf+' -flatten '+fout_png
       print, 'figures saved in: '+fout_png
       if keyword_set(nops) then spawn,'rm '+feps
    endif else  print, 'figures saved in: '+fout

end
