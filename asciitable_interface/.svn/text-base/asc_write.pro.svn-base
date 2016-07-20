;**************************************************************;
;+
;*NAME:
;
;    ASC_WRITE 
;
;*PURPOSE:
;
;    Creates an ASCII table file in which each column represents one
;    of the specified input vectors. Up to 15 parameters are allowed.
;
;*CALLING SEQUENCE:
;
;    ASC_WRITE,FILENAME,P1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15
;
;*PARAMETERS:
;
;    FILENAME  (REQ) (I) (0) (S)
;        Output file name. If no extension is specified, .txt is used.
;
;    P1  (REQ) (I) (1) (F)
;        first input vector to be converted to an ASCII file. The
;        number of elements in P1 determines the size of the ASCII table.
;
;    P2-P15  (OPT) (I) (1) (F)
;        input vectors to be converted to an ASCII file. Each vector will
;        be stored as a column in an ASCII table. Vectors must have the
;        same number of elements as P1 (or more). 
;
;*PROCEDURE:
;
;    ASC_WRITE uses the number of parameters in the procedure call
;    to determine the number of entries to write to each line in
;    the output ASCII file. Parameters are written out using the PRINTF
;    command, in the same order as they are specified in the procedure 
;    call.
;
;*EXAMPLE:
;
;    ASC_WRITE,'MOD1',W,F
;         writes vectors to file MOD1.txt with two parameters per line
;         (i.e., a wavelength and a flux value).
;
;*SUBROUTINES:
;
;    DECOMPOSE
;    PARCHECK
;
;*NOTES:
;
;    Up to 15 parameters can be written to each row or record.
;    There is currently no limit to the size of the input vectors,
;     and there is no checking of the input vector sizes. The number
;     of records written is determined solely from the number of elements
;     in P1. If any other input vector is smaller than P1, an error will 
;     occur. Vectors larger than P1 will be truncated in the output file.
;    The default formats defined by IDL are used to determine the number
;     of significant figures. Therefore, some numbers may be rounded-off.
;
;    tested with IDL version 2.1.0 (sunos sparc)	13 Sep 91
;    tested with IDL version 2.1.0 (vms vax)    	13 Sep 91
;    tested with IDL version 2.1.0 (ultrix mipsel)	N/A
;    
;*MODIFICATION HISTORY:
;
;    written by RWT 12/28/90 
;    1/15/91 RWT use DECOMPOSE
;    1/22/91 PJL transferred to SUN/UNIX, added PARCHECK
;    6/19/91 PJL cleaned up; tested on SUN and VAX; updated prolog
;    7/01/91 RWT increase to 10 parameters
;    7/18/91 RWT increase to 15 parameters and add format keyword for
;            compatibility with IDL version 2. tested on VAX and SUN
;    8/22/91 RWT add a blank between parameters to avoid numbers running
;            together
;    9/13/91 RWT make internal parameter i a longword integer to be
;            compatible with new keyword NEL allowed in ASC_READ.
;-
;****************************************************************
 pro asc_write,filename,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,header=header
;
; check parameters
;
 npar = n_params(0)
 if npar eq 0 then begin
    print,' ASC_WRITE,FILENAME,P1,p2,p3,p4,...p15'
    retall
 endif  ; npar
 parcheck,npar,indgen(15)+2,'ASC_WRITE'
 decompose,filename,d,u,n,e,v
 if (e eq '') then e = '.txt'    ; .txt is default extension
 filename = d + u + n + e + v
;
; initialize variables
;
 a = ''
 i = 0L
 npar = npar - 1
 nrow = n_elements(p1) - 1

;
; open output file 
;
 get_lun,un
 openw,un,filename

;
; Optionally write header
;
 if(n_elements(header) ne 0) then begin
   printf,un,header
 endif
 
;
; write out vectors using printf command
;
 case npar of
    1: while (i le nrow) do begin
          x1 = p1(i) 
          printf,format='(1A)',un,x1
          i = i+1
       endwhile  ; 1
    2: while (i le nrow) do begin
          x1 = p1(i)
          x2 = p2(i) 
          printf,format='(A," ",A)',un,x1,x2
          i = i+1
       endwhile  ; 2
    3: while (i le nrow) do begin
          x1 =p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          printf,format='(3(A," "))',un,x1,x2,x3
          i = i+1
       endwhile  ; 3
    4: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          printf,format='(4(A," "))',un,x1,x2,x3,x4
          i = i+1
       endwhile  ; 4
    5: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          printf,format='(5(A," "))',un,x1,x2,x3,x4,x5
          i = i+1
       endwhile  ; 5
    6: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          printf,format='(6(A," "))',un,x1,x2,x3,x4,x5,x6
          i = i+1
       endwhile  ; 6
    7: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          printf,format='(7(A," "))',un,x1,x2,x3,x4,x5,x6,x7
          i = i+1
       endwhile  ; 7
    8: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          printf,format='(8(A," "))',un,x1,x2,x3,x4,x5,x6,x7,x8
          i = i+1
       endwhile  ; 8
    9: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          x9 = p9(i)
          printf,format='(9(A," "))',un,x1,x2,x3,x4,x5,x6,x7,x8,x9
          i = i+1
       endwhile  ; 9
   10: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          x9 = p9(i)
          x10 = p10(i)
          printf,format='(10(A," "))',un,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
          i = i+1
       endwhile  ; 10
   11: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          x9 = p9(i)
          x10 = p10(i)
          x11 = p11(i)
          printf,format='(11(A," "))',un,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11
          i = i+1
       endwhile  ; 11
   12: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          x9 = p9(i)
          x10 = p10(i)
          x11 = p11(i)
          x12 = p12(i)
          printf,format='(12(A," "))',un,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
          i = i+1
       endwhile  ; 12
   13: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          x9 = p9(i)
          x10 = p10(i)
          x11 = p11(i)
          x12 = p12(i)
          x13 = p13(i)
          printf,format='(13(A," "))', $
              un,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13
          i = i+1
       endwhile  ; 13
   14: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          x9 = p9(i)
          x10 = p10(i)
          x11 = p11(i)
          x12 = p12(i)
          x13 = p13(i)
          x14 = p14(i)
          printf,format='(14(A," "))', $
                 un,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14
          i = i+1
       endwhile  ; 14
   15: while (i le nrow) do begin
          x1 = p1(i) 
          x2 = p2(i) 
          x3 = p3(i) 
          x4 = p4(i) 
          x5 = p5(i) 
          x6 = p6(i) 
          x7 = p7(i) 
          x8 = p8(i)
          x9 = p9(i)
          x10 = p10(i)
          x11 = p11(i)
          x12 = p12(i)
          x13 = p13(i)
          x14 = p14(i)
          x15 = p15(i)
          printf,format='(15(A," "))', $
                 un,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15
          i = i+1
       endwhile  ; 15
    else: print,'Invalid number of parameters'
 endcase  ; npar
;
 free_lun,un
 print,byte(npar),' Columns of ',fix(n_elements(p1)),' points written'
 return
 end  ; asc_write
