      SUBROUTINE PPMECH
 
C     *************************************************************************
C
C     PPMECH
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     27-JUL-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     PREPROCESSES CHEMICAL MECHANISM DATA FOR PPCHEM
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_chemin.h'
      INCLUDE 'com_ppccom.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      CHARACTER*80 ALINE,ENDLIN
      CHARACTER*50 ERRSTR


C     BEGIN
C     =====

C     =========================================================================

C     READ THE MECHANISM FILE
C     -----------------------
      OPEN(UNIT=NCMECH,FILE=FNMECH,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1000)

C     FIVE LINES OF COMMENT ARE IGNORED
C     ---------------------------------
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE

C     =========================================================================

C     PROCESS THE SPECIES LIST
C     ------------------------
      CALL PPMSPC

C     =========================================================================

C     PROCESS THE THIRD BODY LIST
C     ---------------------------
      CALL PPMBDY

C     =========================================================================

C     PROCESS THE MECHANISM LIST
C     --------------------------
      CALL PPMMEC

C     =========================================================================

C     APPLY CONVERSION FACTORS TO THE REACTION RATE DATA
C     ------------------------
      CALL PPMRCF

C     =========================================================================

C     READ THE LAST LINE AND CLOSE THE MECHANISM FILE
C     -----------------------------------------------
      READ(NCMECH,'(A)',END=1010)ENDLIN
      WRITE(NCREPT,'(A)')ENDLIN
      CLOSE(NCMECH)
     
C     =========================================================================

C     MECHANISM FILE ERROR HANDLER
      GOTO 1020
1000  ERRSTR = 'error opening mechanism file'
      CALL ERHAND(ERRSTR,IEFATL)
1010  ERRSTR = 'unexpected end of mechanism file'
      CALL ERHAND(ERRSTR,IEFATL)
1020  CONTINUE
C     END OF MECHANISM FILE ERROR HANDLER

C     =========================================================================


      RETURN
      END
