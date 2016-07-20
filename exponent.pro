   FUNCTION Exponent, axis, index, number

     ; A special case.
     IF number EQ 0 THEN RETURN, '0'

     ; Assuming multiples of 10 with format.
     ex = String(number, Format='(e8.1)')
     pt = StrPos(ex, '.')

     preval=StrMid(ex, 0, pt+2)
     first = string(round(float(preval)))

;     print, 'exponent:, ex, pt,', ex, pt, preval, first

     sign = StrMid(ex, pt+3, 1)
     thisExponent = StrMid(ex, pt+4)

     ; Shave off leading zero in exponent
     WHILE StrMid(thisExponent, 0, 1) EQ '0' DO thisExponent = StrMid(thisExponent, 1)

     ; Fix for sign and missing zero problem.
     IF (Long(thisExponent) EQ 0) THEN BEGIN
        sign = ''
        thisExponent = '0'
     ENDIF

     ;print, first, sign, thisExponent

     ; Make the exponent a superscript.
     IF sign EQ '-' THEN BEGIN
        RETURN, first + 'E' + sign + thisExponent ; + '!N'
     ENDIF ELSE BEGIN
        RETURN, first + 'E' + thisExponent ; + '!N'
     ENDELSE



   END
