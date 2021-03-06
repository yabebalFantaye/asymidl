;
;+
; NAME:
; LS_square
;
; PURPOSE:
; This procedure defines a square to be used as the user-defined plot symbol (PSYM=8) 
;
; CATEGORY:
; Utilities
;
; CALLING SEQUENCE:
; LS_square
;
; INPUTS:
; SZ - the size of the symbol (default=1)
; TK - the symbol thickness (default=1)
; FILL - whether or not to fill the square (default = empty)
;
; KEYWORD PARAMETERS:
; see INPUTS above
;
; OUTPUTS:
; This procedure sets the IDL PSYM[=8] user-defined plot symbol.
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;   Written by: L.D.Spencer 2013/Jan/04
;-
PRO LS_square, SZ=SZ, TK=TK, FILL=FILL,color=color
  ;
  ; This routine allows a circle to be plotted using the IDL user-selected PSYM (PSYM=8)
  ; It can be filled or open.
  ;
  IF ~KEYWORD_SET(SZ) THEN SZ=1.  ; default to one
  IF ~KEYWORD_SET(TK) THEN TK=1.  ; default to one
  IF N_ELEMENTS(FILL) EQ 0 THEN FILL = 0 
  ;
  X = [-.5,-.5,.5,.5,-.5]
  Y = [-.5,.5,.5,-.5,-.5]
  SZ = SZ*2.
  ;
  X = X*SZ
  Y = Y*SZ
  ;
  USERSYM, X, Y, THICK=TK, FILL=FILL,color=color
  ;
END
;
;
