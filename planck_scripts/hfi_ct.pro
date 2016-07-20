;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  Utility files needed to get Planck colourtable loaded, and Healpix mollview plots to work.
;
;
;
;
;
;    Credit to Stephane Columbi for his work on the colour table!!  
;    I just took his code and deleted some of the uneeded bits, and put a wrapper function around it.
;
;
;
;
;
;
;    I have lifted three healpix routines and modified them:  mollview -> LS_mollview, proj2out -> LS_proj2out, and oplot_graticule -> LS_oplot_graticule
;    
;    The version of epstopdf on magique3 does not allow changing the pdf resolution, so I have pointed to a current version (2.18) available from CTAN, 
;    this file must be downloaded and pointed to for this portion of the code to work properly.
;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;============================================================
pro hfi_ct,ascii=ascii,plotlum=plotlum,verbose=verbose, CTDIR=CTDIR, CTfile=CTfile, LOAD=LOAD,NColors=ncolors, Bottom=bottom
;============================================================
; Create the official Planck colour table.
; ascii=filename : to output the color table in a file.
;          This will create a 3 column ascii file, with
;          R for the 1st one, G for the 2nd one and B for the
;          3rd one.
; /plotlum : to plot the luminance of the color table. 
; /verbose : sets verbose mode
; CTDIR:  is the directory to place the colortable file in (it copies the current colortable file, overwrites table 41, then saves it as CTfile
; CTfile:  is the resultant/modified CT file, default is colors1.tbl
; LOAD: Also loadds the modified colour table.
; 
; For more questions about this routine, ask 
; Locke Spencer (LockeSpencerAtastro.cf.ac.uk), or Stephane Colombi (colombiATiap.fr)
;============================================================
;
OS = !VERSION.OS
IF STRCMP(OS,'Win',3) THEN WIN = 1 ELSE WIN = 0
;IF WIN EQ 1 THEN spawn, 'chdir', CurDIR ELSE spawn, 'pwd', CurDIR ; assume all non windows machines accept the `pwd' command
cd,current=CurDir
ds = path_sep()
CurDir += path_sep()

IF N_ELEMENTS(CTDIR) EQ 0 THEN CTDIR = CurDIR
IF N_ELEMENTS(CTfile) EQ 0 THEN CTfile = 'colors1.tbl'
;
; Check for the file color1.tbl in the current directory.
; 
CT_found = file_test(CTDIR+CTfile)
CTsuf = ds+'resource'+ds+'colors'+ds
CTfile_orig = 'colors1.tbl'
IF WIN EQ 1 THEN cpstr = 'copy ' ELSE cpstr = 'cp ' 
IF CT_found EQ 0 THEN spawn, cpstr+!DIR+CTsuf+CTfile_orig+' '+CTDIR+CTfile , CTmsg
;
   npoints=7
;
Rim=lonarr(npoints)
Gim=lonarr(npoints)
Bim=lonarr(npoints)

 Rim(0)=0
 Gim(0)=0
 Bim(0)=255
 
 Rim(1)=0
 Gim(1)=112
 Bim(1)=255

 Rim(2)=0
 Gim(2)=221
 Bim(2)=255

 Rim(3)=255
 Gim(3)=237
 Bim(3)=217

 Rim(4)=255
 Gim(4)=180
 Bim(4)=0

 Rim(5)=255
 Gim(5)=75
 Bim(5)=0

 Rim(6)=100
 Gim(6)=0
 Bim(6)=0

if (keyword_set(verbose)) then begin
   print, 'luminance=',luminance(Rim,Gim,Bim)
   print, 'R=',Rim
   print, 'G=',Gim
   print, 'B=',Bim
endif

nc=256
R=bytarr(nc)
G=bytarr(nc)
B=bytarr(nc)

ipos=floor(((findgen(npoints))/(npoints-1))*255)
ipos(0)=0
ipos(npoints-1)=255

for i=1,npoints-1 do begin
   ip1=ipos(i-1)
   ip2=ipos(i)
   R1=Rim(i-1)
   R2=Rim(i)
   G1=Gim(i-1)
   G2=Gim(i)
   B1=Bim(i-1)
   B2=Bim(i)
   for j=ip1,ip2 do begin
      R(j)=floor(R1+(R2-R1)*double(j-ip1)/double(ip2-ip1))
      G(j)=floor(G1+(G2-G1)*double(j-ip1)/double(ip2-ip1))
      B(j)=floor(B1+(B2-B1)*double(j-ip1)/double(ip2-ip1))
   endfor
endfor

if (keyword_set(ascii)) then begin
   openw, unit, ascii,/get_lun
   for i=0l,255l do begin
      printf,unit, R(i),G(i),B(i)
   endfor
   free_lun,unit
endif
if (keyword_set(verbose)) then begin
   print, 'R final=',R
   print, 'G final=',G
   print, 'B final=',B
endif
;stop

modifyct,41,'parchment1',R,G,B, FILE=CTDIR+CTfile

if (keyword_set(plotlum)) then begin
   loadct,0
   window,0
   plot,luminance(R,G,B),xtitle='index',ytitle='luminance',xr=[0,255],/xs
endif

IF KEYWORD_SET(LOAD) THEN loadct, 41, FILE=CTDIR+CTFILE,NColors=ncolors, Bottom=bottom

end
