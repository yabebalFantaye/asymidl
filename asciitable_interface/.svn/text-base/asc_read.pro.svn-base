;****************************************************************
;+
;*NAME:
;
;    ASC_READ
;
;*PURPOSE:
;
;    Reads an ASCII file containing a table of numbers (up to 15 columns)
;    and outputs the values as floating point vectors.
;
;*CALLING SEQUENCE:
;
;    ASC_READ,FILENAME,P1,p2,p3,p4...,p15,sl=f,ll=l,nel=n,comment_symbol=s
;
;*PARAMETERS:
;
;    FILENAME  (REQ) (I) (0) (S)
;        Input file name. If no extension is specified, .txt is assumed.
;
;    P1  (REQ) (O) (1) (F)
;        Output vector containing data from the first column of the 
;        input table file.
;
;    P2-P15  (OPT) (O) (1) (F)
;        Output vectors containing data from additional columns in table
;        (up to 15 are currently allowed).
;
;    SL (OPT) (KEY) (0) (BILF)
;        Keyword desribing first line of data to be read.
;        (Allows column headings to be skipped.)
;
;    LL (OPT) (KEY) (0) (BILF)
;        Keyword describing last line of data to be read.
;        
;    NEL (OPT) (KEY) (0) (BILF)
;        Keyword describing number of elements in output vector(s).
;        Default = 2500.
;
;   COMMENT_SYMBOL (OPT)
;        Set to eg '#' to skip lines that begin with this symbol.
;
;*PROCEDURE:
;
;    ASC_READ uses the number of parameters in the procedure call
;    to determine the number of entries contained in each line
;    of the input ASCII file. It then skips the appropriate 
;    number of lines and starts extracting data on the following line. 
;    Data is extracted until either an EOF is reached, NEL entries are
;    read, or the specified rows are read. All 0 values at end of
;    initialized arrays are removed.
;
;*EXAMPLE:
;
;    ASC_READ,'MOD1',SL=5,W,F
;         assumes file MOD1.txt contains two parameters per line
;         (i.e., a wavelength and a flux value) starting on the 5th
;         line from the beginning of the file, and outputs all values
;         found as W and F.
;    ASC_READ,'SWP12345',H,W,F,E
;         assumes file SWP12345.txt contains 4 parameters per line
;         (H,W,F,E) with no column headings, and extracts values 
;         as 4 vectors.
;    ASC_READ,'SWP12345',H,W,F,E,SL=5,LL=30
;         skips 4 lines and then reads rows 5 through 30 from the input 
;         ASCII file.
;    ASC_READ,'RADEC.TXT',CAM,IMAGE,RA,DEC,NEL=10000
;         reads first 10,000 entries from input file. In this particular
;         example, the camera number, image number, ra, and dec fields
;         from an INGRES database table were converted to ASCII and copied
;         to a disk file called RADEC.TXT using an SQL script. ASC_READ 
;         could then be used to read the resulting file.
;
;*SUBROUTINES:
;
;    DECOMPOSE
;    PARCHECK
;
;*NOTES:
;
;    Assumes each line of file has 1 to 15 parameters,
;     and each line has the same number of parameters.
;    Default array size is currently 2500 points.
;    The number of entries read depends upon the smaller value of
;     NEL, (LL - SL), or (EOF - SL).
;    Uses the default formats defined by IDL to determine the number
;     of significant figures. Therefore, some numbers may be rounded-off.
;
;    tested with IDL version 2.1.0 (sunos sparc)	13 Sep 91
;    tested with IDL version 2.3.2 (vms vax)		04 Dec 92
;    tested with IDL version 2.1.0 (ultrix mispel)	N/A
;    
;*MODIFICATION HISTORY:
;
;    written by RWT 11/12/90
;    12-11-90 RWT add SKIPL parameter
;    1-15-91 RWT add DECOMPOSE
;    1-22-91 PJL transferred to sun/unix, added PARCHECK
;    6-19-91 PJL cleaned up; tested on SUN and VAX; updated prolog
;    6-28-91 RWT increase to 10 parameters
;    7-18-91 RWT increase to 15 parameters & replace parameter SKIP with
;            keywords SL and LL.
;    9-13-91 RWT add keyword NEL.
;   12-17-91 PJL corrected typo in prolog
;    3-02-92 RWT correct parcheck parameter to allow specifying just
;            1 output parameter (as reported by PJL).
;   12-04-92 RWT allow sl to be > 32768 (i.e., make J a longword integer)
;    7-07-93 RWT Correct error that occurs when 11 parameters are specified
;   30-11-06 Sam Leach, SISSA, modifies to accept comment_symbol.
;   35-04-09 Sam Leach, SISSA, adding silent keyword.
;-
;****************************************************************
pro asc_read,fname,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15, $
             sl=f,ll=l,nel=num,comment_symbol=comment_symbol,$
             outcomment=outcomment,silent=silent
;
; check parameters
;
 npar = n_params(0)
 if npar eq 0 then begin
    print,' ASC_READ,FILENAME,P1,p2,p3,p4,...p15,sl=n,ll=m,nel=num'
    retall
 endif  ; npar
 if n_elements(comment_symbol) eq 0 then comment_symbol='#'
 parcheck,npar,indgen(18)+2,'ASC_READ'
 decompose,fname,d,u,n,e,v
 if (e eq '') then e = '.txt'       ; use .txt if no extension found
 filename = d + u + n + e + v
;
; initialize output arrays
;
 x1 = dblarr(1)
 x2 = x1 & x3 = x1 & x4 = x1 & x5 = x1
 x6 = x1 & x7 = x1 & x8 = x1 & x9 = x1 & x10 = x1
 if keyword_set(num) then vnum = num else vnum = 2500
 p1 = dblarr(vnum)
 p2 = p1 & p3 = p1 & p4 = p1 & p5 = p1 
 p6 = p1 & p7 = p1 & p8 = p1 & p9 = p1 & p10 = p1
 outcomment=make_array(vnum,value= ' ')
 if npar gt 11 then begin
    x11 = x1 & x12 = x1 & x13 = x1 & x14 = x1 & x15 = x1
    p11 = p1 & p12 = p1 & p13 = p1 & p14 = p1 & p15 = p1
    endif
 a = ''
 i = 0L
 npar = npar - 1
;
; determine # lines to skip & # of rows to read
;
 if keyword_set(f) then skipl = (f - 1) > 0 else skipl = 0
 if keyword_set(l) then nrow = (l - skipl) else nrow = vnum 
;
; skip lines containing column headings
;
 get_lun,un
 openr,un,filename
 if (skipl ge 1) then for j=0L,skipl-1 do readf,un,a
;
; read rest of file, store scalar values in array
;
 case npar of
    1: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1
	    p1(i) = x1
            i = i+1
          endif
	endwhile  ; 1
    2: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2
	    p1(i) = x1
	    p2(i) = x2
	    i = i+1
	  endif
       endwhile  ; 2
    3: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
;            help,a
	    reads,a,x1,x2,x3
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    i = i+1
	  endif
       endwhile  ; 3
    4: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    i = i+1
	  endif
	endwhile  ; 4
    5: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    i = i+1
	  endif
       endwhile  ; 5
    6: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin	    
	    reads,a,x1,x2,x3,x4,x5,x6
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    i = i+1
	  endif
       endwhile  ; 6
    7: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    i = i+1
	  endif
       endwhile  ; 7
    8: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7,x8
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    i = i+1
	  endif
       endwhile  ; 8
    9: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7,x8,x9
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    p9(i) = x9
	    i = i+1
	  endif
       endwhile  ; 9
   10: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    p9(i) = x9
	    p10(i) = x10
	    i = i+1
	  endif
       endwhile  ; 10
   11: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    p9(i) = x9
	    p10(i) = x10
	    p11(i) = x11
	    i = i+1
	  endif
	endwhile  ; 11
   12: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    p9(i) = x9
	    p10(i) = x10
	    p11(i) = x11
	    p12(i) = x12
	    i = i+1
	  endif
	endwhile  ; 12
   13: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,un,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    p9(i) = x9
	    p10(i) = x10
	    p11(i) = x11
	    p12(i) = x12
	    p13(i) = x13
	    i = i+1
	  endif
       endwhile  ; 13
   14: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    p9(i) = x9
	    p10(i) = x10
	    p11(i) = x11
	    p12(i) = x12
	    p13(i) = x13
	    p14(i) = x14
	    i = i+1
	  endif
       endwhile  ; 14
   15: while ( (not eof(un)) and (i lt nrow) ) do begin
          readf,un,a
	  pos=strpos(a,comment_symbol,0)
          if (pos[0] gt 0) then outcomment(i) = strmid(a,pos[0]+1)
	  if (pos[0] ne 0) then begin
	    reads,a,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15
	    p1(i) = x1
	    p2(i) = x2
	    p3(i) = x3
	    p4(i) = x4
	    p5(i) = x5
	    p6(i) = x6
	    p7(i) = x7
	    p8(i) = x8
	    p9(i) = x9
	    p10(i) = x10
	    p11(i) = x11
	    p12(i) = x12
	    p13(i) = x13
	    p14(i) = x14
	    p15(i) = x15
	    i = i+1
	  endif
       endwhile  ; 15
    else: print,'Invalid number of parameters'
 endcase  ; npar
;
; remove any padded zeroes & finish up
;
 p1 = p1(0:i-1)
 p2 = p2(0:i-1)
 p3 = p3(0:i-1)
 p4 = p4(0:i-1)
 p5 = p5(0:i-1)
 p6 = p6(0:i-1)
 p7 = p7(0:i-1)
 p8 = p8(0:i-1)
 p9 = p9(0:i-1)
 p10 = p10(0:i-1)
 outcomment = outcomment(0:i-1)
 if (npar gt 10) then begin
    p11 = p11(0:i-1)
    p12 = p12(0:i-1)
    p13 = p13(0:i-1)
    p14 = p14(0:i-1)
    p15 = p15(0:i-1)
    endif
 free_lun,un
 if (i lt 32767) then i = fix(i)

 if(not keyword_set(silent)) then print,byte(npar),' Columns of',i,' points extracted'
 return
 end  ; asc_read
