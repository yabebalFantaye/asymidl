pro planck_mollview_window,out,CTDIR,CTFILE,GR,W,PX,large=large,medium=medium,NColors=ncolors, Bottom=bottom

  FDIR='temp/'
  spawn, 'mkdir -p '+FDIR

         ;
                                ;  Set the new colortable
  ;if keyword_set(hfi_ct) then begin
     CTDIR = FDIR               ;  The CTDIR needs to be a directory that IDL can read/write to.  
     CTFILE = 'Planck_CT.tbl'

     hfi_ct, CTDIR=CTDIR, CTFILE=CTFILE, /LOAD ,NColors=ncolors, Bottom=bottom
  ;endif
;; this tells me the location of the revised colour table, file CTFILE in directory CTDIR, 
  ;; these are also inputs if you have a colourtable file already.  The LOAD keyword then loads the colour table after it has been created.
  ;; Make an outline for the lattitude and longitude lines. This is needed to draw an ellipse around the map [if desired]. 
                                
  Ngrat = 181d                  ; number of points for the additional graticule curves.
  out_ = {COORD:'G',RA:DBLARR(Ngrat),DEC:DBLARR(Ngrat), LINESTYLE:0, PSYM:0, SYMSIZE:0} ; the outline structure accepted by mollview.
                                ; 
  Nout = 2                      ; the outline is done as two half-curves.
  out = REPLICATE(out_,Nout)
  struct_out = out
                                ;
  RA = DBLARR(Ngrat)            ; held constant while DEC changes  ;  The RA and DEC are in healpix/mollview notation.
  DEC = DINDGEN(Ngrat) - 90d    ; bottom to top, -90 to +90 deg.
                                ;
  out[0].RA = RA - 180d         ;  The half curve at -180 deg.
  out[0].DEC = DEC
  out[0].LINESTYLE=0
                                ;
  out[1].RA = RA + 180d         ;  The half curve at +180 deg.
  out[1].DEC = DEC
  out[1].LINESTYLE=0
                                ;
                                ;
  GR_88 = [60,45]               ;  The graticule spacing for the 88 mm figure.
  GR_120 = [60,30]              ;  The graticule spacing for the 120 mm figure.
  GR_180 = [60,30]              ;  The graticule spacing for the 180 mm figure.
                                ;
  W_88 = 8.8d                   ; cm
  W_120 = 12d                   ; cm
  W_180 = 18d                   ; cm
                                ;
  NmPref = '' ;;'PlanckFig_map_columbi1_IDL_'
  NmSuf = 'mm'
  Nm_88  = NmPref + '88' + NmSuf
  Nm_120 = NmPref + '120'+ NmSuf
  Nm_180 = NmPref + '180'+ NmSuf
                                ;
  FigRes = 300d                 ; in dpi.  The paper figures should be at least 300 dpi, according to the A&A author guide.  For Maps, 600 dpi would be better.  
                                ;
  PX_88  =  8.8d/2.54d*FigRes
  PX_120 = 12.0d/2.54d*FigRes
  PX_180 = 18.0d/2.54d*FigRes
                                ;
  !P.FONT = 0
                                ;


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 88 mm figure ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  SZ    = 88d
  SZstr = '88'
  PX  = PX_88
  GR  = GR_88
  W   =  W_88
  Nm  = Nm_88


  if keyword_set(medium) then begin
     SZ   =   120d
     SZstr = '120'
     PX  = PX_120
     GR  = GR_120
     W   =  W_120
     Nm  = Nm_120
  endif

  if keyword_set(large) then begin
     SZ   =   180d
     SZstr = '180'
     PX  = PX_180
     GR  = GR_180
     W   =  W_180
     Nm  = Nm_180
  endif

end
