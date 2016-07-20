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
