;2c;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;+
;
; NAME:
;    REMOVE_TAGS
;       
; PURPOSE:
;    remove the specified tags from input structure
;
; CALLING SEQUENCE:
;    remove_tags, oldstruct, tagnames, newstruct
;
; INPUTS: 
;    oldstruct: the original structure
;    tagnames: the names of tags to be removed (can be an array)
;
; OPTIONAL INPUTS:
;    NONE.
;
; KEYWORD PARAMETERS:
;    NONE.
;       
; OUTPUTS: 
;    newstruct: the new structure without tags.
;
; OPTIONAL OUTPUTS:
;    NONE
;
; CALLED ROUTINES:
;    
; 
; PROCEDURE: 
;    
; 
;
; REVISION HISTORY:
;    ????? Judith Racusin
;    25-OCT-2000 Modified to handle arbitrary tag types. Also error 
;          handling. Erin Scott Sheldon
;       
;                                      
;-                                       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO remove_tags, struct, tagnames, newstruct

  IF n_params() EQ 0 THEN BEGIN 
      print,'Syntax - remove_tags, oldstruct, tagnames, newstruct'
      print
      print,'Use doc_library,"remove_tags"  for more help.'  
      return
  END

  ;; Figure out which tags get removed

  tags=tag_names(struct)
  n=n_elements(tags)
  tagnames=strupcase(tagnames)
  nt=n_elements(tagnames)
  IF nt EQ 1 THEN BEGIN
      t=where(tags NE tagnames[0],nw) 
      CASE nw OF
        0: BEGIN ; there are no more tags after the one is removed
          ;
          print,'-------------------------------------------------------------'
          message,'This would remove all tags! Dummy tag "foo" added to structure',/inf
          print,'-------------------------------------------------------------'
          newstruct=create_struct('foo',' ')
          return
          ;
        END
        n: BEGIN ;  THe numbwer of tags is the same with and without, therefore no match
          print,'-----------------------------------------------------'
          message,'Tag did not match, structure unchanged',/inf
          print,'-----------------------------------------------------'
          newstruct = struct
          return
        END
        ELSE: ; Maybe something needed here, but think is OK blank for now.
      ENDCASE
;      IF nw EQ n THEN BEGIN
;          print,'-----------------------------------------------------'
;          message,'Tag did not match, structure unchanged',/inf
;          print,'-----------------------------------------------------'
;          newstruct = struct
;          return
;      ENDIF ELSE BEGIN  ;  There is only one tag, and it is desired to be removed.  add a dummy tag
;      ENDELSE
  ENDIF ELSE BEGIN 
      match,tags,tagnames,m
      IF m[0] EQ -1 THEN BEGIN
          print,'-------------------------------------------------'
          message,'No tags matched, structure unchanged',/inf
          print,'-------------------------------------------------'
          newstruct=struct
          return
      ENDIF 
      nm=n_elements(m)
      IF nm EQ n THEN BEGIN 
          print,'-------------------------------------------------------------'
          message,'This would remove all tags! structure unchanged',/inf
          print,'-------------------------------------------------------------'
          newstruct=struct
          return
      ENDIF 
      t=lindgen(n)
      remove, m, t
  ENDELSE 
      
  ;; create new structure
  tags=tags[t]
  n=n_elements(tags)

  newstruct=create_struct(tags[0],struct[0].(t[0]))
  
  FOR i=1L, n-1 DO newstruct = create_struct(temporary(newstruct), $
                                             tags[i], struct[0].(t[i]) )

  newstruct=replicate( temporary(newstruct), n_elements(struct) )
  struct_assign,struct,newstruct

  return
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;
;
;
;+
; NAME:
; LS_DecRound
;
; PURPOSE:
; This procedure rounds a number to a user-defined number of decimal places and outputs this number as a string (and other optional formats).
; It also rounds to a defined amount of significant digits in scientific notation with the /SCI keyword. 
;
; CATEGORY:
; Utilities
;
; CALLING SEQUENCE:
; LS_DecRound, Num, DEC=3, STR=strNum, etc. 
;
; INPUTS:
;   Num - the number (or array of numbers) to be rounded and converted to a string.
;   Dec - the number of decimal places / digits to include, DEC=2 gets  123.45, DEC=-2 gets 100, etc.
;   STR - the string output, appropriately rounded / truncated.
;
; KEYWORD PARAMETERS:
;   SCI - display the number in scientific notation
;   NOSCI - force the exclusion of scientific notation presentation
;   RNDSCI - round in scientific notation (i.e. DEC = sigFigs rather than decimal number), but present as other wise requested (NOSCI keyword)
;   WASSCI - output to indicate if the number would nomrally be presented in scientific notation in a `print, num' statement.
;   ALLOWSCI - allow scientific notation if it is the default, but do not force it.
;
; OUTPUTS:
; This procedure sets the STR keyword with a string (or array of strings) containing the desired number of decimal places / significant figures.
;
; RESTRICTIONS:
;
;
; MODIFICATION HISTORY:
;   Written by: L.D.Spencer 2011/2012
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
;-
PRO LS_DecRound, Num, DEC=DEC, DBL=DBL, FLT=FLT, STR=STR, SCI=SCI, NOSCI=NOSCI, RNDSCI=RNDSCI, WASSCI=WASSCI, ALLOWSCI=ALLOWSCI
;
; This procedure takes a number (float or double) and rounds it to a certain number of decimal places, returning a float, double, or string.
;
; The string will be truncated at the specified number of decimal places to exclude any trailing zeroes.
;
; KEYWORDS:
;           DEC   -   Number of decimal places to round to
;           DBL   -   Output variable for rounded double
;           FLT   -   Output variable for rounded float
;           STR   -   Output variable for rounded string 
;           
;           SCI   -   Set to 1 to force scientific notation output
;           NOSCI -   Force output to be in standard notation (i.e. not scientific)   ;;;;    SCI takes precedence over NOSCI
;           ALLOWSCI- allow scientific notation if the string conversion does this, but do not externally force it.  
;           RNDSCI  - Round the number in scientific notation rather than as a decimal.  
;           WASSCI-   Indicates if it would have been scientific notation anyways
;
IF N_ELEMENTS(NUM) EQ 0 THEN NUM = !C_ ; 1.23456789
IF N_ELEMENTS(DEC) EQ 0 THEN DEC = 1  ; default is 1 decimal place
IF N_ELEMENTS(SCI) EQ 0 THEN SCI = 0
IF N_ELEMENTS(NOSCI) EQ 0 THEN NOSCI = 0
IF N_ELEMENTS(ALLOWSCI) EQ 0 THEN ALLOWSCI = 1  ; default is to allow scientific as it comes up.
;IF N_ELEMENTS(RNDSCI) EQ 0 THEN RNDSCI = 0 ; If the number comes up as sci normally, then choose to round based on this if this keyword is not set.
IF N_ELEMENTS(RNDSCI) EQ 0 THEN RNDIND = 1 ELSE RNDIND = 0
IF N_ELEMENTS(RNDSCI) EQ 0 THEN RNDSCI = 0
;
IF ((SCI EQ 1) AND (NOSCI EQ 1)) THEN NOSCI = 0   ; I should generate an error here also, but I'll leave it for later. 
;
DecNum = 10d^(DOUBLE(DEC))
;
Nels = N_ELEMENTS(Num)
Nds = N_ELEMENTS(DEC)
;
DBL = DBLARR(Nels)
FLT = FLTARR(Nels)
STR = STRARR(Nels)
;
WasSCI = INTARR(Nels)
;
FOR i = 0, Nels - 1 DO BEGIN
  ;
  IF Nds EQ 1 THEN BEGIN
    DEC_ = DEC
    DecNum_ = DecNum 
  ENDIF ELSE BEGIN
    DEC_ = DEC[i]
    DecNum_ = DecNum[i]
  ENDELSE
  ;
  Numi = Num[i]
  stri = STRING(Numi)
  ;
  ;   First round the number as required.   ; RNDIND means that RNDSCI was not set globally, so set it individually...
  IF RNDIND EQ 1 THEN BEGIN ; check if the number comes up to be scientific by itself, then set the RNDSCI keyword for each number in the loop.
    HasExp = STRPOS(stri,'e')
    IF HasExp LT 0 THEN RNDSCI = 0 ELSE RNDSCI = 1
  ENDIF
  IF RNDSCI EQ 1 THEN BEGIN
    ; get the base number and exponents separated
    strsci = STRING(Numi, FORMAT='(e)')
    Lsci = STRLEN(strsci)
    SciExp = STRPOS(strsci,'e')
    eSTR = STRMID(strsci,(SciExp + 1),(Lsci - (SciExp + 1)))
    eVal = DOUBLE(eSTR)
    baseVal = DOUBLE(STRMID(strsci,0,SciExp - 1))
    DBLsci = DOUBLE(ROUND(baseVal*DecNum_))/DecNum_*(10d^eVal)
    RNDi = DBLsci
    ;FLTsci = FLOAT(DBLsci)
    ;stringsci = STRMID(STRTRIM(STRING(DOUBLE(ROUND(baseVal*DecNum))/DecNum),2),0,2 + DEC)+'e'+eSTR
  ENDIF ELSE BEGIN
    ; check to see if the number can be rounded normally (i.e. is it within the limits of a 64 bit integer?)
    NotTooBig = DOUBLE(Numi*DecNum_) LT DOUBLE(9223372036854775807ll)
    NotTooSmall=DOUBLE(Numi*DecNum_) GT DOUBLE(-9223372036854775808ll)
    IF ((NotTooBig) AND (NotTooSmall)) THEN BEGIN
      DBLi = DOUBLE(ROUND(Numi*DecNum_, /L64))/DecNum_
      RNDi = DBLi
    ENDIF ELSE BEGIN  ; will do sci round anyways
      ; get the base number and exponents separated
      strsci = STRING(Numi, FORMAT='(e)')
      Lsci = STRLEN(strsci)
      SciExp = STRPOS(strsci,'e')
      eSTR = STRMID(strsci,(SciExp + 1),(Lsci - (SciExp + 1)))
      eVal = DOUBLE(eSTR)
      baseVal = DOUBLE(STRMID(strsci,0,SciExp - 1))
      DBLsci = DOUBLE(ROUND(baseVal*DecNum_))/DecNum_*(10d^eVal)
      RNDi = DBLsci
    ENDELSE
  ENDELSE
  ; Have now rounded the number as requested by the user (natural vs. scientific rounding)
  ; Now must format the output as requested (mostly a problem for the string)
  ;
  IF RNDi LT 0d THEN NEGstr = '-' ELSE NEGstr = ''
  IF RNDi LT 0d THEN NEGspc = 1 ELSE NEGspc = 0
  ;
  IF DEC_ LT 8 THEN strsci = STRTRIM(STRING(FLOAT(RNDi), FORMAT='(e)'),2) ELSE strsci = STRTRIM(STRING(RNDi, FORMAT='(e)'),2)
  Lsci = STRLEN(strsci)
  SciExp = STRPOS(strsci,'e')
  eSTR = STRMID(strsci,(SciExp + 1),(Lsci - (SciExp + 1)))
  eVal = DOUBLE(eSTR)
  baseVal = STRMID(strsci,0,SciExp - 1)
  IF RNDSCI EQ 1 THEN BEGIN ; round to DEC number of decimals in the scientific notation, easy...    
    stringi = STRMID(baseVal,0,2 + NEGspc + DEC_)+'e'+eSTR
  ENDIF ELSE BEGIN  ; display the sci notation of the number rounded in normal format, i.e. figure out how many zeroes there should be.
    DECsci = FIX(eVal + DEC_)
    IF DECsci LT 8 THEN strsci = STRTRIM(STRING(FLOAT(RNDi), FORMAT='(e)'),2) ELSE strsci = STRTRIM(STRING(RNDi, FORMAT='(e)'),2)
    Lsci = STRLEN(strsci)
    SciExp = STRPOS(strsci,'e')
    eSTR = STRMID(strsci,(SciExp + 1),(Lsci - (SciExp + 1)))
    eVal = DOUBLE(eSTR)
    baseVal = STRMID(strsci,0,SciExp - 1)
    stringi = STRMID(baseVal,0,2 + NEGspc + DECsci)+'e'+eSTR
  ENDELSE
  STRsci = stringi
  ;
  RNDlen = NEGspc + 1 + FIX(eVAL) + 1 + DEC_
  IF eVAL LT 0d THEN RNDlen = NEGspc + 1 + FIX(ABS(eVAL)) - 1 + 1 + DEC_
  IF eVAL GT 0 THEN BEGIN
    IF RNDlen LT 8 THEN BEGIN 
      strnormi=STRTRIM(STRING(FLOAT(RNDi)),2)
    ENDIF ELSE BEGIN
      IF RNDlen LT 10 THEN BEGIN
        strnormi=STRTRIM(STRING((RNDi)),2)
      ENDIF ELSE BEGIN
        strnormi=STRTRIM(STRING(RNDi, FORMAT='(d64)'),2)  ; this may need fixed still, but it should get me out of my current fix...
      ENDELSE
    ENDELSE
  ENDIF ELSE BEGIN
    IF RNDlen LT 8 THEN BEGIN 
      strnormi=STRTRIM(STRING(FLOAT(RNDi), FORMAT='(d)'),2)
    ENDIF ELSE BEGIN
      IF RNDlen LT 10 THEN BEGIN
        strnormi=STRTRIM(STRING((RNDi), FORMAT='(d)'),2)
      ENDIF ELSE BEGIN
        strnormi=STRTRIM(STRING(RNDi, FORMAT='(d64)'),2)  ; this may need fixed still, but it should get me out of my current fix...
      ENDELSE
    ENDELSE
  ENDELSE
  Strdot = STRPOS(strnormi,'.')
  IF RNDSCI EQ 1 THEN BEGIN ; the number was rounded in sci notation, but is to be displayed in normal format, determine how to do this.
      ; get the sci exponent and DEC number to determine the correct display
      NumDEC = FIX(eVal - DEC_) 
      IF NumDEC GE 0 THEN strnormi = STRMID(strnormi,0,Strdot - 1) ELSE strnormi = STRMID(strnormi,0, strdot + 1 + ABS(NumDEC))
  ENDIF ELSE BEGIN  ; the number was rounded in decimal format, and must now be displayed in decimal format, no problem.
    NumDEC = DEC_ ; FIX(-DEC_)
    IF NumDEC LE 0 THEN strnormi = STRMID(strnormi,0,Strdot) ELSE strnormi = STRMID(strnormi,0,Strdot + DEC_ + 1)
    ;strnormi = STRMID(strnormi,0,Strdot + DEC_ + 1)
    ;stop
  ENDELSE
  ; 
  stri = STRING(RNDi)
  HasExp = STRPOS(stri,'e')
  DBLi = RNDi
  FLTi = FLOAT(RNDi)
  IF HasExp GT 0 THEN BEGIN ; The number is in sci format
    WASSCI[i] = 1
    STRi = STRsci
;    Lsci = STRLEN(stri)
;    SciExp = STRPOS(stri,'e')
;    eSTR = STRMID(stri,(SciExp + 1),(Lsci - (SciExp + 1)))
;    eVal = DOUBLE(eSTR)
;    baseVal = STRMID(stri,0,SciExp - 1)
;    IF RNDSCI EQ 1 THEN stringi = STRMID(baseVal,0,2 + DEC_)+'e'+eSTR $
;    ELSE stringi = STRMID(baseVal,0,2 + DEC_)+'e'+eSTR  ; FIXME this should be displayed differently I think as it rounds the string to sci decimals rather than natural ones, i.e. the two roundings may be different.
  ENDIF ELSE BEGIN  ; The number is not in sci format
    WasSCI[i] = 0
    STRi = STRnormi
;    stringi = STRTRIM(stri,2)
;    ; get rid of end zeroes
;    Strdot = STRPOS(stringi,'.')
;    IF RNDSCI EQ 0 THEN stringi = STRMID(stringi,0,Strdot + DEC_ + 1) $
;    ELSE stringi = STRMID(stringi,0,Strdot + DEC_ + 1) ; FIXME If I did scientific rounding then this needs to be fixed.
  ENDELSE
  ;STRi = stringi
  ;
  ;   Now select the string to include
  IF NOSCI EQ 1 THEN STRi = strnormi
  IF SCI EQ 1 THEN STRi = strsci
  ;
  ;stop
  DBL[i] = DBLi
  FLT[i] = FLTi
  STR[i] = STRi
ENDFOR
 
END
;
;


;+
;PURPOSE
; to rotate the orientation of y-axis plot tick-mark-labels for Planck/HFI figures by 90 deg. ccw.
; This should work for output to a file (i.e. .eps) or to the screen (e.g. for tvread to .png)
; returns  the same result as the plot command, with the y-axis labels rotated.
;SYNTAX
; HFI_plot, x, y, ..., [DECMODX=], [DECMODY=], [Y_DX=], [Y_DY=], [X_DX=], [X_DY=]
;INPUTS
; x, y, similar to the plot command
; _extra: extra keywords for plot()
;KEYWORDS
; DECMODX - set this to an integer to increase (or decrease) the number of digits in the x-axis tick mark labels (default=0).
; DECMODY - set this to an integer to increase (or decrease) the number of digits in the y-axis tick mark labels (default=0).
; Y_DX    - set this to manually shift the y-axis tick labels horizontally, in data coordinates.
; Y_DY    - set this to manually shift the y-axis tick labels vertically, in data coordinates.
; X_DX    - set this to manually shift the x-axis tick labels horizontally, in data coordinates.
; X_DY    - set this to manually shift the x-axis tick labels vertically, in data coordinates.
;
;Written by L.D.Spencer, Dec. 2012 / Jan. 2013
;
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
;-

;+
;PURPOSE
; to rotate the orientation of y-axis plot tick-mark-labels for
;Planck/HFI figures by 90 deg. ccw.
; This should work for output to a file (i.e. .eps) or to the screen
;(e.g. for tvread to .png)
; returns  the same result as the plot command, with the y-axis labels
;rotated.
;SYNTAX
; HFI_plot, x, y, ..., [DECMODX=], [DECMODY=], [Y_DX=], [Y_DY=],
;[X_DX=], [X_DY=]
;INPUTS
; x, y, similar to the plot command
; _extra: extra keywords for plot()
;KEYWORDS
; DECMODX - set this to an integer to increase (or decrease) the
;           number of digits in the x-axis tick mark labels
;(default=0).
; DECMODY - set this to an integer to increase (or decrease) the
;           number of digits in the y-axis tick mark labels
;(default=0).
; Y_DX    - set this to manually shift the y-axis tick labels
;           horizontally, in data coordinates.
; Y_DY    - set this to manually shift the y-axis tick labels
;           vertically, in data coordinates.
; X_DX    - set this to manually shift the x-axis tick labels
;           horizontally, in data coordinates.
; X_DY    - set this to manually shift the x-axis tick labels
;           vertically, in data coordinates.
;
;Written by L.D.Spencer, Dec. 2012 / Jan. 2013
;
;
;   
;   This program is free software: you can redistribute it and/or
;modify
;   it under the terms of the GNU General Public License as published
;by
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
;-
;+
;PURPOSE
; to rotate the orientation of y-axis plot tick-mark-labels for
;Planck/HFI figures by 90 deg. ccw.
; This should work for output to a file (i.e. .eps) or to the screen
;(e.g. for tvread to .png)
; returns  the same result as the plot command, with the y-axis labels
;rotated.
;SYNTAX
; HFI_plot, x, y, ..., [DECMODX=], [DECMODY=], [Y_DX=], [Y_DY=],
;[X_DX=], [X_DY=]
;INPUTS
; x, y, similar to the plot command
; _extra: extra keywords for plot()
;KEYWORDS
; DECMODX - set this to an integer to increase (or decrease) the
;           number of digits in the x-axis tick mark labels
;(default=0).
; DECMODY - set this to an integer to increase (or decrease) the
;           number of digits in the y-axis tick mark labels
;(default=0).
; Y_DX    - set this to manually shift the y-axis tick labels
;           horizontally, in data coordinates.
; Y_DY    - set this to manually shift the y-axis tick labels
;           vertically, in data coordinates.
; X_DX    - set this to manually shift the x-axis tick labels
;           horizontally, in data coordinates.
; X_DY    - set this to manually shift the x-axis tick labels
;           vertically, in data coordinates.
;
;Written by L.D.Spencer, Dec. 2012 / Jan. 2013
;
;
;   
;   This program is free software: you can redistribute it and/or
;modify
;   it under the terms of the GNU General Public License as published
;by
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
;-
;+
;PURPOSE
; to rotate the orientation of y-axis plot tick-mark-labels for
;Planck/HFI figures by 90 deg. ccw.
; This should work for output to a file (i.e. .eps) or to the screen
;(e.g. for tvread to .png)
; returns  the same result as the plot command, with the y-axis labels
;rotated.
;SYNTAX
; HFI_plot, x, y, ..., [DECMODX=], [DECMODY=], [Y_DX=], [Y_DY=],
;[X_DX=], [X_DY=]
;INPUTS
; x, y, similar to the plot command
; _extra: extra keywords for plot()
;KEYWORDS
; DECMODX - set this to an integer to increase (or decrease) the
;           number of digits in the x-axis tick mark labels
;(default=0).
; DECMODY - set this to an integer to increase (or decrease) the
;           number of digits in the y-axis tick mark labels
;(default=0).
; Y_DX    - set this to manually shift the y-axis tick labels
;           horizontally, in data coordinates.
; Y_DY    - set this to manually shift the y-axis tick labels
;           vertically, in data coordinates.
; X_DX    - set this to manually shift the x-axis tick labels
;           horizontally, in data coordinates.
; X_DY    - set this to manually shift the x-axis tick labels
;           vertically, in data coordinates.
;
;Written by L.D.Spencer, Dec. 2012 / Jan. 2013
;
;
;   
;   This program is free software: you can redistribute it and/or
;modify
;   it under the terms of the GNU General Public License as published
;by
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
;-







PRO hfi_plot, x, y, _EXTRA=_EXTRA, DECMODX=DECMODX, DECMODY=DECMODY, $
  Y_DX=Y_DX, Y_DY=Y_DY, X_DX=X_DX, X_DY=X_DY, $
  YTTL_DX=YTTL_DX, YTTL_DY=YTTL_DY, XTTL_DX=XTTL_DX, XTTL_DY=XTTL_DY, $
  DEBUG=DEBUG, XTICK_GET=XT, YTICK_GET=YT, XLGNAT=XLGNAT, YLGNAT=YLGNAT, FIXMINUS=FIXMINUS
;
; This script is to plot to a virtual device to get the y-axis labels, and allow them to be re-plotted with a 90^o anti-clockwise rotation.
; ;
; 
;  If no keywords are passed then it doesn't like it, so make a structure to use as a placeholder
IF N_ELEMENTS(_EXTRA) EQ 0 THEN _EXTRA = {HFI:0}
;
IF N_ELEMENTS(y) EQ 0 THEN BEGIN  ;  I think that x was forgotten, generate a new one.
  x_ = DINDGEN(N_ELEMENTS(x))
  y = x
  x = x_
ENDIF
IF N_ELEMENTS(DECMODX) EQ 0 THEN DECMODX = 0
IF N_ELEMENTS(DECMODY) EQ 0 THEN DECMODY = 0
IF N_ELEMENTS(Y_DX) EQ 0 THEN Y_DX = 0d
IF N_ELEMENTS(Y_DY) EQ 0 THEN Y_DY = 0d
IF N_ELEMENTS(X_DX) EQ 0 THEN X_DX = 0d
IF N_ELEMENTS(X_DY) EQ 0 THEN X_DY = 0d
IF N_ELEMENTS(YTTL_DX) EQ 0 THEN YTTL_DX = 0d
IF N_ELEMENTS(YTTL_DY) EQ 0 THEN YTTL_DY = 0d
IF N_ELEMENTS(XTTL_DX) EQ 0 THEN XTTL_DX = 0d
IF N_ELEMENTS(XTTL_DY) EQ 0 THEN XTTL_DY = 0d
IF N_ELEMENTS(XLGNAT) EQ 0 THEN XLGNAT = 0
IF N_ELEMENTS(YLGNAT) EQ 0 THEN YLGNAT = 0
; stop ;   
IF N_ELEMENTS(y) EQ 0 THEN y = DINDGEN(5) + 1d
Ny = N_ELEMENTS(y)
IF N_ELEMENTS(x) EQ 0 THEN x = dindgen(Ny)
;y = (x + 3);/1d-1
;
; check for some alternate names to my keywords:  I shortened XRANGE to XR, XSTYLE to XS, etc.
_EXTRA_orig = _EXTRA
;
IF TAG_EXIST(_EXTRA, 'XRANGE') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'XR', _EXTRA.XRANGE)
  _EXTRA = temporary(_E)
ENDIF
IF TAG_EXIST(_EXTRA, 'XSTYLE') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'XS', _EXTRA.XSTYLE)
  _EXTRA = temporary(_E)
ENDIF
IF TAG_EXIST(_EXTRA, 'XSTY') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'XS', _EXTRA.XSTY)
  _EXTRA = temporary(_E)
ENDIF
IF TAG_EXIST(_EXTRA, 'XST') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'XS', _EXTRA.XST)
  _EXTRA = temporary(_E)
ENDIF
;
;
IF TAG_EXIST(_EXTRA, 'YRANGE') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'YR', _EXTRA.YRANGE)
  _EXTRA = temporary(_E)
ENDIF
IF TAG_EXIST(_EXTRA, 'YSTYLE') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'YS', _EXTRA.YSTYLE)
  _EXTRA = temporary(_E)
ENDIF
IF TAG_EXIST(_EXTRA, 'YSTY') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'YS', _EXTRA.YSTY)
  _EXTRA = temporary(_E)
ENDIF
IF TAG_EXIST(_EXTRA, 'YST') EQ 1 THEN BEGIN
  _E = create_struct(_EXTRA, 'YS', _EXTRA.YST)
  _EXTRA = temporary(_E)
ENDIF
;
;IF TAG_EXIST(_EXTRA, 'YRANGE') EQ 1 THEN _EXTRA = create_struct(temporary(_EXTRA), 'YR', _EXTRA.YRANGE)
;IF TAG_EXIST(_EXTRA, 'YSTYLE') EQ 1 THEN _EXTRA = create_struct(temporary(_EXTRA), 'YS', _EXTRA.YSTYLE)
;IF TAG_EXIST(_EXTRA, 'YSTY')   EQ 1 THEN _EXTRA = create_struct(temporary(_EXTRA), 'YS', _EXTRA.YSTY)
;IF TAG_EXIST(_EXTRA, 'YST')    EQ 1 THEN _EXTRA = create_struct(temporary(_EXTRA), 'YS', _EXTRA.YST)
;
IF TAG_EXIST(_EXTRA, 'YLOG') EQ 1 THEN YLG = _EXTRA.YLOG ELSE YLG = 0
IF TAG_EXIST(_EXTRA, 'XLOG') EQ 1 THEN XLG = _EXTRA.XLOG ELSE XLG = 0
CHSZ_Tag = TAG_EXIST(_EXTRA, 'CHARSIZE')
IF CHSZ_TAG THEN BEGIN
  CHSZ = _EXTRA.CHARSIZE
  CHSZ_ = _EXTRA.CHARSIZE 
ENDIF ELSE BEGIN
  CHSZ = 0
  CHSZ_ = 1d
ENDELSE
;
YS_ = (TAG_EXIST(_EXTRA, 'YS')); OR (TAG_EXIST(_EXTRA, 'YSTYLE')))
XS_ = (TAG_EXIST(_EXTRA, 'XS'))
;
YR_ = (TAG_EXIST(_EXTRA, 'YR'))
XR_ = (TAG_EXIST(_EXTRA, 'XR'))
;
XTICKINT = (TAG_EXIST(_EXTRA,'XTICKINTERVAL'))
YTICKINT = (TAG_EXIST(_EXTRA,'YTICKINTERVAL'))
;
XTCKS = (TAG_EXIST(_EXTRA,'XTICKS'))
YTCKS = (TAG_EXIST(_EXTRA,'YTICKS'))
;
XTKV  = (TAG_EXIST(_EXTRA, 'XTICKV'))
YTKV  = (TAG_EXIST(_EXTRA, 'YTICKV'))

;
IF YS_ THEN BEGIN
  YS_bin = STRING(_EXTRA.YS, FORMAT='(B0)')
  YS_4 = STRMID('000000'+YS_bin, 2, 1, /REVERSE_OFFSET)
  IF YS_4 EQ '1' THEN YST = _EXTRA.YS ELSE YST = _EXTRA.YS + 4
ENDIF ELSE BEGIN
  YST = 4
ENDELSE
;
;
IF XS_ THEN BEGIN
  XS_bin = STRING(_EXTRA.XS, FORMAT='(B0)')
  XS_4 = STRMID('000000'+XS_bin, 2, 1, /REVERSE_OFFSET)
  IF XS_4 EQ '1' THEN XST = _EXTRA.XS ELSE XST = _EXTRA.XS + 4
ENDIF ELSE BEGIN
  XST = 4
ENDELSE
;
IF YR_ THEN YRN = _EXTRA.YR ELSE YRN=[min(y),max(y)]
IF XR_ THEN XRN = _EXTRA.XR ELSE XRN=[min(x),max(x)]
;
IF XTICKINT THEN XTINT = _EXTRA.XTICKINTERVAL ELSE XTINT = 0d
IF YTICKINT THEN YTINT = _EXTRA.YTICKINTERVAL ELSE YTINT = 0d
;
IF YTCKS THEN YTS = _EXTRA.YTICKS ELSE YTS=0
IF XTCKS THEN XTS = _EXTRA.XTICKS ELSE XTS=0
;
IF YTKV THEN YVS = _EXTRA.YTICKV ELSE YVS=0
IF XTKV THEN XVS = _EXTRA.XTICKV ELSE XVS=0
;
;stop
;
print, 'nx, xtickv',XTS,XVS
print, 'ny, ytickv',YTS,YVS



plot, x, y, YLOG=YLG, XLOG=XLG, /NODATA, /NOERASE, YS=YST, XS=XST, YTICK_GET = YT, CHARSIZE=CHSZ, XTICK_GET=XT, XR=XRN, YR=YRN, $
  XTICKINTERVAL=XTINT, YTICKINTERVAL=YTINT, XTICKS=XTS, YTICKS=YTS, XTICKV=XVS, YTICKV=YVS
; The above will plot no axes, and no lines.  It just gets the ytick values.
; , YS=4, XS=4;, YTICKFORMAT='(A1)' ; just to get Y tick values...



; FIXME: I need to check the width of the character string to be sure that it doesn't overlap with neighboring data points.
;        If so then take every other data point.
Dev_to_DataY = (!Y.CRANGE[1] - !Y.CRANGE[0])/DOUBLE(!P.CLIP[3] - !P.CLIP[1])
Dev_to_DataX = (!X.CRANGE[1] - !X.CRANGE[0])/DOUBLE(!P.CLIP[2] - !P.CLIP[0])


CH_YSZ = !D.Y_CH_SIZE*CHSZ_ ; MAX([!D.X_CH_SIZE,!D.Y_CH_SIZE])*CHSZ_
CH_XSZ = !D.X_CH_SIZE*CHSZ_

IF TAG_EXIST(_EXTRA,'YCHARSIZE') THEN CH_YSZ = CH_YSZ*_EXTRA.YCHARSIZE
IF TAG_EXIST(_EXTRA,'XCHARSIZE') THEN CH_XSZ = CH_XSZ*_EXTRA.XCHARSIZE
;
LS_DecRound, YT, DEC=NumDEC, STR=YT_STR, SCI=SCI, NOSCI=NOSCI, RNDSCI=RNDSCI, WASSCI=WASSCI, ALLOWSCI=ALLOWSCI ; FIXME, need to check all of this. 
LS_DecRound, XT, DEC=NumDEC, STR=XT_STR, SCI=SCI, NOSCI=NOSCI, RNDSCI=RNDSCI, WASSCI=WASSCI, ALLOWSCI=ALLOWSCI ; FIXME, need to check all of this. 

CH_lenY = STRLEN(YT_STR)
CH_lenX = STRLEN(XT_STR) 



print, 'hfi_plot: len(YSTR), len(XSTR)',CH_lenY, CH_lenX



;
NumY = N_ELEMENTS(YT)

;center the label at the tick mark 
;;YT = YT - CH_lenY*CH_YSZ*Dev_to_DataY

YTICKNAME = STRARR(NumY) + ' '

Ytick_Xval = !X.crange[0]   
Ytick_Xval = DBLARR(NumY) + Ytick_Xval
Ytick_Yval = YT
;
NumX = N_ELEMENTS(XT)

;center the label at the tick mark 
;;XT = XT - CH_lenX*CH_XSZ*Dev_to_DataX/2d

XTICKNAME = STRARR(NumX) + ' '

Xtick_Yval = !Y.crange[0] ;  - (!Y.crange[1] - !Y.crange[0])/100d
Xtick_Yval = DBLARR(NumX) + Xtick_Yval
Xtick_Xval = XT

print, '!XY.crange: ',!X.crange,!Y.crange
print, 'hfi_plot: Xtick_Xval, Xtick_Yval',Xtick_Xval, Xtick_Yval
print, 'hfi_plot: Ytick_Xval, Ytick_Yval',Ytick_Xval, Ytick_Yval


;
IF YLG THEN BEGIN
  ;  DO some things differently as it is a log plot
  YTICKNAME = STRARR(NumY) + '  '
  YT_exp = ALOG10(YT)
  YT_exp_DEC = 0 + DECMODY
  LS_DecRound, YT_exp, DEC=YT_exp_DEC, STR=YT_exp_str_
;  YT_str = '10!U'+STRTRIM(STRING(YT_exp),2)+'!N'
  YT_str = '10!U'+YT_exp_str_+'!N'
  IF KEYWORD_SET(FIXMINUS) THEN BEGIN
    FOR Minus_i = 0, N_ELEMENTS(YT_STR)-1 DO BEGIN
      YT_STR[Minus_i] = StrJoin(StrSplit(YT_STR[Minus_i], '-', /Regex, /Extract, /Preserve_Null), '!M-!X!U')
    ENDFOR
  ENDIF
  IF KEYWORD_SET(YLGNAT) THEN BEGIN
    ;  Do not want the axis labels in 10^n notation.
    YT_EXP_DEC = CEIL((-1d)*YT_EXP) + DECMODY
    LS_DecRound, YT, DEC=YT_EXP_DEC, STR=YT_STR, SCI=SCI, NOSCI=NOSCI, RNDSCI=RNDSCI, WASSCI=WASSCI, ALLOWSCI=ALLOWSCI
    ;
    IF KEYWORD_SET(FIXMINUS) THEN BEGIN
      FOR Minus_i = 0, N_ELEMENTS(YT_STR)-1 DO BEGIN
        YT_STR[Minus_i] = StrJoin(StrSplit(YT_STR[Minus_i], '-', /Regex, /Extract, /Preserve_Null), '!M-!X')
      ENDFOR
    ENDIF
    NegYTicks = WHERE(YT LT 0d, Nyneg)
    IF Nyneg GT 0 THEN YT_STR[NegYTicks] = YT_STR[NegYTicks]+' '
  ENDIF
ENDIF ELSE BEGIN
  ; 
  dY = (!Y.crange[1] - !Y.crange[0])/DOUBLE(NumY - 1d)
  ;
  ; FIXME, include a check for /YLOG being set!!
  NumDec = CEIL(ALOG10(dY)*(-1d)) + DECMODY ;+ 1
  ;If NumDec LT 0 THEN NumDec = 0
  LS_DecRound, YT, DEC=NumDEC, STR=YT_STR, SCI=SCI, NOSCI=NOSCI, RNDSCI=RNDSCI, WASSCI=WASSCI, ALLOWSCI=ALLOWSCI ; FIXME, need to check all of this. 
  ;
  IF KEYWORD_SET(FIXMINUS) THEN BEGIN
    FOR Minus_i = 0, N_ELEMENTS(YT_STR)-1 DO BEGIN
      YT_STR[Minus_i] = StrJoin(StrSplit(YT_STR[Minus_i], '-', /Regex, /Extract, /Preserve_Null), '!M-!X')
    ENDFOR
  ENDIF
  NegYTicks = WHERE(YT LT 0d, Nyneg)
  IF Nyneg GT 0 THEN YT_STR[NegYTicks] = YT_STR[NegYTicks]+' '
  ;
ENDELSE
;


IF XLG THEN BEGIN
  ;  DO some things differently as it is a log plot
  ;XT_exp = FIX(ALOG10(XT))
  XT_exp = ALOG10(XT)
  ;XT_str = '10!U'+STRTRIM(STRING(XT_exp),2)+'!N'
  XT_exp_DEC = 0 + DECMODX
  LS_DecRound, XT_exp, DEC=XT_exp_DEC, STR=XT_exp_str_
  XT_str = '10!U'+XT_exp_str_+'!N'
  IF KEYWORD_SET(FIXMINUS) THEN BEGIN
    FOR Minus_i = 0, N_ELEMENTS(XT_STR)-1 DO BEGIN
      XT_STR[Minus_i] = StrJoin(StrSplit(XT_STR[Minus_i], '-', /Regex, /Extract, /Preserve_Null), '!M-!X!U')
    ENDFOR
  ENDIF
  IF KEYWORD_SET(XLGNAT) THEN BEGIN
    ;  Do not want the axis labels in 10^n notation.
    XT_EXP_DEC = CEIL((-1d)*XT_EXP) + DECMODX
    LS_DecRound, XT, DEC=XT_EXP_DEC, STR=XT_STR, SCI=SCI, NOSCI=NOSCI, RNDSCI=RNDSCI, WASSCI=WASSCI, ALLOWSCI=ALLOWSCI
    IF KEYWORD_SET(FIXMINUS) THEN BEGIN
      FOR Minus_i = 0, N_ELEMENTS(XT_STR)-1 DO BEGIN
        XT_STR[Minus_i] = StrJoin(StrSplit(XT_STR[Minus_i], '-', /Regex, /Extract, /Preserve_Null), '!M-!X')
      ENDFOR
    ENDIF
    NegXTicks = WHERE(XT LT 0d, Nxneg)
    IF Nxneg GT 0 THEN XT_STR[NegXTicks] = XT_STR[NegXTicks]+' '
  ENDIF
ENDIF ELSE BEGIN

  ;
  dX = (!X.crange[1] - !X.crange[0])/DOUBLE(NumX - 1d)
  ;
  ; FIXME, include a check for /YLOG being set!!
  NumDec = CEIL(ALOG10(dX)*(-1d)) + DECMODX ;+ 1
  ;If NumDec LT 0 THEN NumDec = 0
  LS_DecRound, XT, DEC=NumDEC, STR=XT_STR, SCI=SCI, NOSCI=NOSCI, RNDSCI=RNDSCI, WASSCI=WASSCI, ALLOWSCI=ALLOWSCI ; FIXME, need to check all of this. 
  ;
  IF KEYWORD_SET(FIXMINUS) THEN BEGIN
    FOR Minus_i = 0, N_ELEMENTS(XT_STR)-1 DO BEGIN
      XT_STR[Minus_i] = StrJoin(StrSplit(XT_STR[Minus_i], '-', /Regex, /Extract, /Preserve_Null), '!M-!X')
    ENDFOR
  ENDIF
  NegXTicks = WHERE(XT LT 0d, Nxneg)
  IF Nxneg GT 0 THEN XT_STR[NegXTicks] = XT_STR[NegXTicks]+' '
  ;
ENDELSE


;stop
;
_EXTRA_orig = _EXTRA
IF TAG_EXIST(_EXTRA,'XTITLE') THEN BEGIN
  remove_tags, _EXTRA, 'XTITLE', _EXTRA_
  _EXTRA = _EXTRA_
ENDIF 
;
IF TAG_EXIST(_EXTRA,'YTITLE') THEN BEGIN
  remove_tags, _EXTRA, 'YTITLE', _EXTRA_
  _EXTRA = _EXTRA_
ENDIF 
; 
ytickname_given=0
IF TAG_EXIST(_EXTRA,'YTICKNAME') THEN BEGIN
   ytickname_given=1
  ;YTICKNAME=_EXTRA.YTICKNAME
  ;print,'ystr before',YTICKNAME
  remove_tags, _EXTRA, 'YTICKNAME', _EXTRA_
  _EXTRA = _EXTRA_

ENDIF 

   xtickname_given=0
IF TAG_EXIST(_EXTRA,'XTICKNAME') THEN BEGIN
   xtickname_given=1
  ;XTICKNAME=_EXTRA.XTICKNAME
  ;print,'xstr before',XTICKNAME
  remove_tags, _EXTRA, 'XTICKNAME', _EXTRA_
  _EXTRA = _EXTRA_
ENDIF 



print,'hfi_plot: ystr to show',YT_STR
print,'hfi_plot: xstr to show',XT_STR

;
plot, x, y, _EXTRA=_EXTRA, YTICKNAME=YTICKNAME, XTICKNAME=XTICKNAME, YTITLE=' ', XTITLE=' ' ;  
;
_EXTRA = _EXTRA_orig




;
ITERY = 1
ITERX = 1


CH_strtsY = Ytick_Yval 
IF YLG THEN BEGIN 
  ;
  ;CH_lenY = (CH_lenY - 4d - 2d)*0.62d + 2d
  IF KEYWORD_SET(YLGNAT) THEN CH_lenY = CH_lenY ELSE CH_lenY = (CH_lenY - 4d - 2d)*0.62d + 2d
  IF KEYWORD_SET(FIXMINUS) THEN BEGIN ; the non-printing characters are changing the string length calculations so they need removed from the count.
    ;
    NegExpY = WHERE(ALOG10(Ytick_Yval) LT 0d, NnegExpY)
    IF ((NnegExpY GT 0) AND (~KEYWORD_SET(YLGNAT))) THEN CH_lenY[NegExpY] = CH_lenY[NegExpY] - 6d*0.62d
    ;
  ENDIF
  ;
  ;Dev_to_DataY = (ALOG10(!Y.CRANGE[1]) - ALOG10(!Y.CRANGE[0]))/DOUBLE(!P.CLIP[3] - !P.CLIP[1])
  Dev_to_DataY = (!Y.CRANGE[1] - !Y.CRANGE[0])/DOUBLE(!P.CLIP[3] - !P.CLIP[1])
  ;
  CH_strtsY = ALOG10(Ytick_Yval) - (CH_lenY)*CH_YSZ*Dev_to_DataY/2d
  ;CH_strtsY = ALOG10(Ytick_Yval - (CH_lenY)*CH_YSZ*Dev_to_DataY/2d)
  CH_strtsY = 10d^CH_strtsY + Y_DY
  ;
  CH_endsY  = ALOG10(Ytick_Yval) + (CH_lenY)*CH_YSZ*Dev_to_DataY/2d  
  ;CH_endsY  = ALOG10(Ytick_Yval + (CH_lenY)*CH_YSZ*Dev_to_DataY/2d)
  CH_endsY = 10d^CH_endsY + Y_DY
  ;
ENDIF ELSE BEGIN
  ;
  Dev_to_DataY = (!Y.CRANGE[1] - !Y.CRANGE[0])/DOUBLE(!P.CLIP[3] - !P.CLIP[1])
  ;
  print, 'y dev_to_data',Dev_to_DataY

  IF KEYWORD_SET(FIXMINUS) THEN BEGIN ; the non-printing characters are changing the string length calculations so they need removed from the count.
    ;
    NegY = WHERE(Ytick_Yval LT 0d, NnegY)
    IF (NnegY GT 0) THEN CH_lenY[NegY] = CH_lenY[NegY] - 4d
    ;
  ENDIF
  ;

  CH_strtsY= Ytick_Yval - CH_lenY*CH_YSZ*Dev_to_DataY + Y_DY ;;make the fist point remain unchanged
  ;
  CH_endsY  = Ytick_Yval + CH_lenY*CH_YSZ*Dev_to_DataY/2d + Y_DY 
  ;
ENDELSE

print, 'line 1218 hfi_plot'


CH_strtsX = Xtick_Xval
IF XLG THEN BEGIN 
  ;
  ;CH_lenX = (CH_lenX - 4d - 2d)*0.62d + 2d
  IF KEYWORD_SET(XLGNAT) THEN CH_lenX = CH_lenX ELSE CH_lenX = (CH_lenX - 4d - 2d)*0.62d + 2d
  IF KEYWORD_SET(FIXMINUS) THEN BEGIN ; the non-printing characters are changing the string length calculations so they need removed from the count.
    ;
    NegExpX = WHERE(ALOG10(Xtick_Xval) LT 0d, NnegExpX)
    IF ((NnegExpX GT 0) AND (~KEYWORD_SET(XLGNAT))) THEN CH_lenX[NegExpX] = CH_lenX[NegExpX] - 6d*0.62d
    ;
  ENDIF
  ;
  ;Dev_to_DataX = (ALOG10(!X.CRANGE[1]) - ALOG10(!X.CRANGE[0]))/DOUBLE(!P.CLIP[2] - !P.CLIP[0])
  Dev_to_DataX = (!X.CRANGE[1] - !X.CRANGE[0])/DOUBLE(!P.CLIP[2] - !P.CLIP[0])
  ;
  CH_strtsX = ALOG10(Xtick_Xval) - CH_lenX*CH_XSZ*Dev_to_DataX/2d
  ;CH_strtsX = ALOG10(Xtick_Xval - CH_lenX*CH_XSZ*Dev_to_DataX/2d)
  CH_strtsX = 10d^CH_strtsX + X_DX
  ;
  CH_endsX  = ALOG10(Xtick_Xval) + CH_lenX*CH_XSZ*Dev_to_DataX/2d  
  ;CH_endsX  = ALOG10(Xtick_Xval + CH_lenX*CH_XSZ*Dev_to_DataX/2d)
  CH_endsX  = 10d^CH_endsX + X_DX
  ; 
ENDIF ELSE BEGIN
  ;
  Dev_to_DataX = (!X.CRANGE[1] - !X.CRANGE[0])/DOUBLE(!P.CLIP[2] - !P.CLIP[0])
  ;
  IF KEYWORD_SET(FIXMINUS) THEN BEGIN ; the non-printing characters are changing the string length calculations so they need removed from the count.
    ;
    NegX = WHERE(Xtick_Xval LT 0d, NnegX)
    IF (NnegX GT 0) THEN CH_lenX[NegX] = CH_lenX[NegX] - 4d
    ;
  ENDIF
  ;
  
  CH_strtsX = Xtick_Xval - CH_lenX*CH_XSZ*Dev_to_DataX/2d + X_DX ;;make the fist point remain unchanged
  ;
  CH_endsX  = Xtick_Xval + CH_lenX*CH_XSZ*Dev_to_DataX/2d + X_DX 
  ;
ENDELSE
;


print, 'hfi_plot: CH_strtsY, CH_endsY',CH_strtsY, CH_endsY
print, 'hfi_plot: CH_strtsX, CH_endsX',CH_strtsX, CH_endsX


;
IF NumY GT 1 THEN OlapY = CH_strtsY[1:NumY - 1] - CH_endsY[0:NumY - 2] ELSE OlapY = 1d
If NumX GT 1 THEN OlapX = CH_strtsX[1:NumX - 1] - CH_endsX[0:NumX - 2] ELSE OlapX = 1d
;
NegOlapY = WHERE(OlapY LE 0d, NnegOlapY)
NegOlapX = WHERE(OlapX LE 0d, NnegOlapX)

;
;;if not ytickname_given then begin
   IF NnegOlapY GT 0 THEN ITERY = 2 ELSE ITERY = 1

;  Check that the last number for either axis does not go outside of the plot range.  
   IF TAG_EXIST(_EXTRA,'YMARGIN') THEN YMAR = _EXTRA.YMARGIN ELSE YMAR = !Y.MARGIN
;
   YBOX = !Y.CRANGE + [-1d,1d]*DOUBLE(!D.Y_CH_SIZE)*DOUBLE(YMAR)/DOUBLE(!D.Y_SIZE)/!Y.S[1]
   IF YLG THEN YBOX = 10d^(!Y.CRANGE + [-1d,1d]*DOUBLE(!D.Y_CH_SIZE)*DOUBLE(YMAR)/DOUBLE(!D.Y_SIZE)/!Y.S[1])
   IF CH_strtsY[0] LT YBOX[0] THEN YT_str[0] = ' '
   IF CH_strtsY[0] LT YBOX[0] THEN CH_strtsY_ = CH_strtsY[1:*] ELSE CH_strtsY_ = CH_strtsY
   IF CH_strtsY[0] LT YBOX[0] THEN CH_endsY_  = CH_endsY[1:*] ELSE CH_endsY_ = CH_endsY
   IF CH_strtsY[0] LT YBOX[0] THEN Yfirst = 1 ELSE Yfirst = 0
;;
   IF CH_endsY[NumY - 1] GT YBOX[1] THEN YT_str[NumY - 1] = ' '
   IF CH_endsY[NumY - 1] GT YBOX[1] THEN CH_endsY_  = CH_endsY_[0:N_ELEMENTS(CH_endsY_) - 2]
   IF CH_endsY[NumY - 1] GT YBOX[1] THEN CH_strtsY_ = CH_strtsY_[0:N_ELEMENTS(CH_strtsY_) - 2]
   
;;endif else CH_strtsY_ = CH_strtsY


;;if not xtickname_given then begin
   IF NnegOlapX GT 0 THEN ITERX = 2 ELSE ITERX = 1
;  Check that the last number for either axis does not go outside of the plot range.  
   IF TAG_EXIST(_EXTRA,'XMARGIN') THEN XMAR = _EXTRA.XMARGIN ELSE XMAR = !X.MARGIN
;;
   XBOX = !X.CRANGE + [-1d,1d]*DOUBLE(!D.X_CH_SIZE)*DOUBLE(XMAR)/DOUBLE(!D.X_SIZE)/!X.S[1]
   IF XLG THEN XBOX = 10d^(!X.CRANGE + [-1d,1d]*DOUBLE(!D.X_CH_SIZE)*DOUBLE(XMAR)/DOUBLE(!D.X_SIZE)/!X.S[1])
   IF CH_strtsX[0] LT XBOX[0] THEN XT_str[0] = ' '
   IF CH_strtsX[0] LT XBOX[0] THEN CH_strtsX_ = CH_strtsX[1:*] ELSE CH_strtsX_ = CH_strtsX
   IF CH_strtsX[0] LT XBOX[0] THEN CH_endsX_  = CH_endsX[1:*] ELSE CH_endsX_ = CH_endsX
   IF CH_strtsX[0] LT XBOX[0] THEN Xfirst = 1 ELSE Xfirst = 0
;;
   IF CH_endsX[NumX - 1] GT XBOX[1] THEN XT_str[NumX - 1] = ' '
   IF CH_endsX[NumX - 1] GT XBOX[1] THEN CH_endsX_  = CH_endsX_[0:N_ELEMENTS(CH_endsX_) - 2]
   IF CH_endsX[NumX - 1] GT XBOX[1] THEN CH_strtsX_ = CH_strtsX_[0:N_ELEMENTS(CH_strtsX_) - 2]
   

;;endif
;
;stop
;


;;Xtick_Yval = Xtick_Yval - CH_YSZ*Dev_to_DataY*1.5d   ; Need to subtract one character as xyouts takes the bottom of the letter location.
IF YLG THEN Xtick_Yval = 10d^(Xtick_Yval - CH_YSZ*Dev_to_DataY*1.5d) ELSE Xtick_Yval = Xtick_Yval - CH_YSZ*Dev_to_DataY*1.5d ; Need to subtract one character as xyouts takes the bottom of the letter location.

;;IF YLG THEN Xtick_Yval = 10d^(Xtick_Yval - CH_YSZ*Dev_to_DataY*1.5d) ELSE
;;Ytick_Xval = Ytick_Xval - CH_YSZ*Dev_to_DataX*0.5d  ; It is the coords for the bottom of the characters
IF XLG THEN Ytick_Xval = 10d^(Ytick_Xval - CH_YSZ*Dev_to_DataX*0.5d) ELSE Ytick_Xval = Ytick_Xval - CH_YSZ*Dev_to_DataX*0.5d ; It is the coords for the bottom of the characters
;;IF XLG THEN Ytick_Xval = 10d^Ytick_Xval

;; if xtickname_given eq 0 then begin
;; endif

;
;;Check to see if there is overlap between the first x and first y labels
IF (Ytick_Xval + Y_DX)[0] GT (CH_strtsX_)[0] THEN BEGIN  ; the first y axis tick label may be printed over the first x axis tick label
  ;
  IF YLG THEN XticksTop = (10d^( ALOG10(Xtick_Yval) + CH_YSZ*Dev_to_DataY))[0] ELSE XticksTop = (Xtick_Yval + CH_YSZ*Dev_to_DataY)[0]
  IF CH_strtsY_[0] LT XticksTop THEN YT_str[Yfirst] = ' '
  ;
ENDIF



;; print,'Yx', (Ytick_Xval + Y_DX)[0:*:ITERY]
;; print,'Yy', (Ytick_Yval + Y_DY)[0:*:ITERY]
;; print, 'Yname', YT_STR[0:*:ITERY]

;; print,'Xx', (Xtick_Xval + X_DX)[0:*:ITERX]
;; print,'Xy', (Xtick_Yval + X_DY)[0:*:ITERX]
;; print, 'Xname', XT_STR[0:*:ITERX]

;YABEBAL COMMENTED OUT THIS
IF ~TAG_EXIST(_EXTRA,'YTICKNAME') THEN xyouts, (Ytick_Xval + Y_DX)[0:*:ITERY], (Ytick_Yval + Y_DY)[0:*:ITERY], YT_STR[0:*:ITERY], ALIGNMENT=0.5, ORIENTATION=90d, /DATA, _EXTRA=_EXTRA ELSE xyouts, (Ytick_Xval + Y_DX), (Ytick_Yval + Y_DY), _EXTRA.YTICKNAME, ALIGNMENT=0.5, ORIENTATION=90d, /DATA, _EXTRA=_EXTRA
IF ~TAG_EXIST(_EXTRA,'XTICKNAME') THEN xyouts, (Xtick_Xval + X_DX)[0:*:ITERX], (Xtick_Yval + X_DY)[0:*:ITERX], XT_STR[0:*:ITERX], ALIGNMENT=0.5, ORIENTATION=0d, /DATA, _EXTRA=_EXTRA ELSE xyouts, (Xtick_Xval + X_DX), (Xtick_Yval + X_DY), _EXTRA.XTICKNAME, ALIGNMENT=0.5, ORIENTATION=0d, /DATA, _EXTRA=_EXTRA
;
;;YABEBAL ADDED THIS
 ;; xyouts, (Ytick_Xval + Y_DX)[0:*:ITERY], (Ytick_Yval + Y_DY)[0:*:ITERY], YT_STR[0:*:ITERY], ALIGNMENT=0.5, ORIENTATION=90d, /DATA, _EXTRA=_EXTRA 
 ;; xyouts, (Xtick_Xval + X_DX)[0:*:ITERX], (Xtick_Yval + X_DY)[0:*:ITERX], XT_STR[0:*:ITERX], ALIGNMENT=0.5, ORIENTATION=0d, /DATA, _EXTRA=_EXTRA 



; Now position the axis labels, including the Y[/X]TTL_DX[/Y] KEYWORDS.
; First determine where I think they should be anyways.
; 
;w_pix = !D.X_SIZE - ((!X.range)[0] + (!X.range)[1])*!D.X_CHSIZE    ;   height of plot range in device coordinates
;w_dat = !X.crange[1] - !X.crange[0] ; height in data units.   ; height of plot-region in data units
;txtWidth = (!D.X_CH_SIZE)*w_dat/w_pix        ; The approximate width of one character...with an extra 2 characters on either side
;;
;h_pix = !D.Y_SIZE - ((!Y.range)[0] + (!Y.range)[1])*!D.Y_CH_SIZE    ;   height of plot range in device coordinates
;h_dat = !Y.crange[1] - !Y.crange[0] ; height in data units.   ; height of plot-region in data units
;txtHeight = (!D.Y_CH_SIZE)*h_dat/h_pix        ; Calculate the average character height.
;
;  Now figure out where xttl and yttl should be, I want yttl centred vertically, and offset 2.5 character spaces from the axis.
;  I want xttl centred horizontally, and vertically offset 2.5 character spaces.
;  or perhaps three spaces...
;  
IF YLG THEN yttl_y = 10d^((!Y.crange[1] + !Y.crange[0])/2d) ELSE yttl_y = ((!Y.crange[1] + !Y.crange[0])/2d) 
IF XLG THEN yttl_x = 10d^(!X.CRANGE[0])*10d^((-1d)*CH_XSZ*Dev_to_DataX*3.5d) ELSE yttl_x = (!X.CRANGE[0]) - CH_XSZ*Dev_to_DataX*3d
;
IF XLG THEN xttl_x = 10d^((!X.crange[1] + !X.crange[0])/2d) ELSE xttl_x = ((!X.crange[1] + !X.crange[0])/2d)
IF YLG THEN xttl_y = 10d^(!Y.CRANGE[0])*10d^((-1d)*CH_YSZ*Dev_to_DataY*2.75d) ELSE xttl_y = (!Y.CRANGE[0]) - CH_YSZ*Dev_to_DataY*2.75d
;
;stop
;

IF TAG_EXIST(_EXTRA,'XTITLE') THEN xyouts, xttl_x + xttl_dx, xttl_y + xttl_dy, ALIGNMENT=0.5, /DATA, _EXTRA.XTITLE, _EXTRA=_EXTRA
IF TAG_EXIST(_EXTRA,'YTITLE') THEN xyouts, yttl_x + yttl_dx, yttl_y + yttl_dy, ALIGNMENT=0.5, ORIENTATION=90d, /DATA, _EXTRA.YTITLE, _EXTRA=_EXTRA
;
IF KEYWORD_SET(DEBUG) THEN BEGIN
  print, '    '
  print, 'YTICKS were supposed to be printed at x=', YTICK_XVAL[0],' and y=',YTICK_YVAL
  print, 'The Y ticks were supposed to be :', YT_STR
  IF TAG_EXIST(_EXTRA,'YTITLE') THEN print, 'YTITLE of "',_EXTRA.YTITLE,'" was supposed to be printed at x=',yttl_x + yttl_dx,' and y=',yttl_y + yttl_dy
  print, '    '
  print, 'XTICKS were supposed to be printed at x=', XTICK_XVAL,' and y=',XTICK_YVAL[0]
  print, 'The X ticks were supposed to be :', XT_STR
  IF TAG_EXIST(_EXTRA,'XTITLE') THEN print, 'XTITLE of "',_EXTRA.XTITLE,'" was supposed to be printed at x=',xttl_x + xttl_dx,' and y=',xttl_y + xttl_dy
  print, '     '
ENDIF


;stop
;
;CLP = !P.CLIP
;Xsz = !D.X_SIZE
;Ysz = !D.Y_SIZE
;;
;Y1 = CLP[1]
;Y2 = CLP[3]
;Y_ = DOUBLE(Y2 - Y1)/DOUBLE(Ysz)
;;
;X1 = CLP[0]
;X2 = CLP[2]
;X_ = DOUBLE(X2 - X1)/DOUBLE(Xsz)
;;
;OrthScl = DOUBLE(Y2 - Y1)/DOUBLE(X2 - X1)
;;
;T3D, /RESET, ROTATE=[0,0,90]
;T3D, TRANSLATE=[DOUBLE(Y2)/DOUBLE(Ysz),(-1d)*DOUBLE(X1)/DOUBLE(Xsz),0d]
;T3D, SCALE = [1d,Y_/X_,1d]
;T3D, SCALE = [OrthScl*Y_/X_,1d,1d]
;T3D, TRANSLATE=[DOUBLE(X1)/DOUBLE(Xsz),0d,0d]
;T3D, TRANSLATE=[0d,DOUBLE(Y1)/DOUBLE(Ysz),0d]
;;
;axis, XAXIS=1, XRANGE=[!Y.crange], color=2, XS=1, /T3D, XTICKV=YT, CHARSIZE=1d/OrthScl*Y_/X_/1.15d
;;xyouts, Ytick_Xval, Ytick_Yval, STRTRIM(STRING(FIX(YT)),2), ALIGNMENT=0.5, ORIENTATION=90d, /DATA
;;
;stop
;;
END
