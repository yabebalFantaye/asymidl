PRO WINDOW_TO_PS,WINDOW,FNAME
;************************************************************************
;+
;*NAME: 
;   WINDOW_TO_PS
;  
;*CLASS:
;   IMAGE PROCESSING
;  
;*CATEGORY:
;  
;*PURPOSE:
;   Outputs the contents of a window to a postscript file that can be printed
; on a postscipt printer.
;  
;*CALLING SEQUENCE:
;   WINDOW_TO_PS, WINDOW, FNAME
;  
;*PARAMETERS:
;   WINDOW - (OPT) (I)(0)(B) - the window number that is to be printed
;   FNAME  - (OPT) (I)(0)(S) - the file name the window is to be saved under
;   
;*INTERACTIVE INPUT:
;   none
;  
;*FILES USED:
;   Creates a .ps file
;  
;*SYSTEM VARIABLES USED:
;   !D
;  
;*SUBROUTINES CALLED:
;   none
;  
;*SIDE EFFECTS:
;   creates .ps file
;  
;*RESTRICTIONS:
;   Works only with 'X' window displays on IDL V2. Tested on VAXstation 3100
;  
;*NOTES:
;   Some Postscipt options can be modified to suit the users needs. See the 
; IDL V2 User's Guide for avaliable options.
;   
;  
;*PROCEDURE:
;   The procedure reads the image off the window with the TVRD option, the
; plotting device 'PS' is selected and the image is "printed" to a postscript
; file. The user then can print the '.ps' file to a postscript printer.
;
;*EXAMPLES:
;   [Enter IDL V2]
;   [Load the image and display it on a window, say window 0]
;   IDL>window_to_ps, 0, test.ps
;   [THE PROGRAM WILL THEN SAVE THE WINDOW TO A POSTSCIPT FILE]
;   [RUNNING TIME 1/4 TO 3 MINUTES]
;   IDL>
;   [PRINT THE POSTSCRIPT FILE TEST.PS TO A POSTSCRIPT PRINTER]  
;     
;*MODIFICATION HISTORY:
;   Written and documented by C. Scott Merkle, Aug 1990  C.A.S.A. - R.D.A.F
;-
;************************************************************************
IF N_PARAMS(0) GT 1 THEN WSET, WINDOW
IF N_PARAMS(0) LT 2 THEN FNAME = 'IDL.PS'
X_SIZE=!D.X_SIZE
Y_SIZE=!D.Y_SIZE
PRINT, 'SAVING SCREEN TO ARRAY'
SCREEN=TVRD(0,0,X_SIZE,Y_SIZE)
SET_PLOT,'PS'
DEVICE,BITS_PER_PIXEL=8,/LANDSCAPE,FILENAME=FNAME
PRINT,'SAVING SCREEN TO FILE'
TV,SCREEN
DEVICE,/CLOSE
SET_PLOT,'X'
RETURN
END

