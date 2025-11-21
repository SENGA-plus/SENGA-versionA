      PROGRAM PPCHEM
 
C     *************************************************************************
C
C     PPCHEM
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
C     PREPROCESSES CHEMICAL DATA FOR SENGA
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_chemin.h'
      INCLUDE 'com_ppccom.h'
C     -------------------------------------------------------------------------


C     BEGIN
C     =====

C     =========================================================================

C     ANNOUNCEMENT
C     ------------
      WRITE(NCREPT,*)'========================================'
      WRITE(NCREPT,*)
      WRITE(NCREPT,*)'PPCHEM: chemistry preprocessor for SENGA'
      WRITE(NCREPT,*)
      WRITE(NCREPT,*)'========================================'
      WRITE(NCREPT,*)

C     =========================================================================

C     GET THE FILENAMES AND CHECK THE FILES
C     -------------------------------------
      CALL FCHECK

C     =========================================================================

C     INITIALISE THE LISTS AND TABLES
C     -------------------------------
      CALL INLIST

C     =========================================================================

C     PROCESS THE MECHANISM FILE
C     --------------------------
      CALL PPMECH

C     =========================================================================

C     PROCESS THE SPECIES DATA
C     ------------------------
      CALL PPSPEC

C     =========================================================================

C     CREATE THE LISTS
C     ----------------
      CALL LISTEM

C     =========================================================================

C     WRITE OUT THE LISTS
C     -------------------
      CALL OUTPUT

C     =========================================================================

      STOP
      END
