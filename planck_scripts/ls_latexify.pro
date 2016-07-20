;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
;
;
;+
;PURPOSE
; to use PSfrag to do latex syntax correctly
;SYNTAX
; LS_latexify, filename, tags, tex, scale
;INPUTS
; filename: name of the eps file
; tags: name of placeholder text string in the .eps file
; tex: the tex you want to replace it
; scale: the scale you want the latex to be [default to 1] (i.e. latex 11pt font) 
;NOTES:
; requires the following latex packages properly installed:
;   geometry, graphicx, psfrag
; requires you to use font=0 in generating the IDL .eps figures.
;
; Follows from discussion by Saurav
;   http://www.astrobetter.com/idl-psfrag/
; Writes and deletes files called LS_HFIfig_temp.*
; The LS version is based on code produced by R. da Silva, UCSC, 2011.
; 
;EXAMPLE
;
;   SET_PLOT, 'PS'
;   !P.FONT=0
;   DEVICE, FILENAME='figure.eps'
;   plot,findgen(10),xtitle='xtitle'
;   DEVICE, /CLOSE
;
;   LS_latexify, 'figure.eps', 'xttl', 'position   [$\mu$m]', FDIR='/some/directory/where/your/file/is/'
; 
; latexify.pro written by R. da Silva, UCSC, 1-4-11 
; LS_latexify written by L.D. Spencer (based on latexify.pro mentioned in previous line).
; 
;   
;   This program is free software: you can redistribute it and/or modify
;   it under the terms of the GNU General Public License as published by
;   the Free Software Foundation, either version 3 of the License, or
;   (at your option) any later version.
;
;   This program is distributed in the hope that it will be useful,
;   but WITHOUT ANY WARRANTY; without even the implied warranty of
;   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;   GNU General Public License for more details.
;
;   A copy of the GNU General Public License is available at 
;   <http://www.gnu.org/licenses/>.
;   
;   Copyright Locke D. Spencer, 2013
;   
; 
; 
;-
pro LS_latexify, filename, tag, tex, scale, outname=outname, keep=keep,$
      height=height, width=width, full=full, FDIR=FDIR, POSN=POSN, PSPOSN=PSPOSN
  ;
  ;IF N_ELEMENTS(POSN) EQ 0 THEN
  ;
  if ~keyword_set(scale) then scale=replicate(1, n_elements(tag))
  scale = strcompress(string(scale), /remove_all)
  ;
  if keyword_set(outname) EQ 0 then begin
    outname='LS_latex_temp.eps'
    noname=1
 endif
  decompose,filename,d,path,name_orig,ext,v
  name = name_orig+'_temp'

  openw, lun, path+name+'.tex', /get_lun
  printf, lun, '\documentclass{article}'
  printf, lun, '\usepackage{geometry, graphicx, psfrag}'
  printf, lun, '\renewcommand{\familydefault}{\sfdefault}'
  printf, lun,'\usepackage{helvet}'
  printf, lun,'\pagestyle{empty}'
  printf, lun,'\geometry{paperwidth='+strcompress(string(width+0.7), /remove_all)+'cm,'+$
                      'paperheight='+strcompress(string(height+0.1), /remove_all)+'cm,margin=1pt}'
  printf, lun,'\begin{document}'

  for i=0, n_elements(tag)-1 do $
        printf, lun,'\psfrag{'+tag[i]+'}[cc][cc]['+scale[i]+']{{'+tex[i]+'}}'
  printf, lun,'\includegraphics[width='+$
       strcompress(string(width), /remove_all)+'cm,height='+strcompress(string(height), /remove_all)+'cm]{'+name_orig+ext+'}'
  printf, lun,'\end{document}'
  close, lun
  free_lun, lun
  ;
  spawn, 'cd '+path+' && latex '+name+'.tex', msg0
  spawn, 'cd '+path+' && dvips -o '+name+'.ps '+name+'.dvi', msg1
  spawn, 'cd '+path+' && ps2epsi '+name+'.ps '+name+'.epsi', msg2
  spawn, " perl -ne 'print unless /^%%BeginPreview/../^%%EndPreview/'<"+path+name+'.epsi > '+outname, msg4  
  

  spawn, " perl -ne 'print unless /^%%BeginPreview/../^%%EndPreview/'<"+path+name+'.ps > '+outname, msg4

  print, 'perl latex done. Entering .. '+path
  print, path+name+'.ps'+' exists = ',file_test(path+name+'.ps',/regular)
  print, 'copy '+path+name+'.ps    to    '+outname
  spawn, 'cd '+path+' && cp '+name+'.ps '+outname,msgx

  ;
  ; Check to see if the file *.eps file already exists, if so rename the old one so the perl script will work.
  ; 
  if keyword_set(keep) then begin
     flCheck = file_test(outname)
     decompose,outname,d,pathout,nameout,extout,v
     IF flCheck THEN spawn, 'cd '+pathout+' && mv -f '+nameout+extout+' old_'+nameout+extout, msg3
  endif
  ;spawn, " && perl -ne 'print unless /^%%BeginPreview/../^%%EndPreview/' <"+FDIR+'LS_HFIfig_temp.epsi > '+FDIR+outname, msg4


  ;
  ; Clean up the temporary files now.
  ;
  if keyword_set(noname) then begin
    spawn, 'mv LS_latex_temp.eps '+filename
    outname=filename
  endif
  ;stop
  spawn, 'cd '+path+' && rm -f '+name+'*'
  ;

end
