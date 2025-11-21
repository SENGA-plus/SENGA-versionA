      SUBROUTINE FCHECK
 
C     *************************************************************************
C
C     FCHECK
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     CHECKS MECHANISM AND SPECIES FILES FOR PPCHEM
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
      CHARACTER*50 ERRSTR


C     BEGIN
C     =====

C     =========================================================================

C     GET THE FILENAMES AND CHECK THE FILES
C     -------------------------------------
      WRITE(NCREPT,*)'Please enter the name of the mechanism file:'
      READ(NCINPT,'(A)')FNMECH
      WRITE(NCREPT,'(A)')FNMECH
      WRITE(NCREPT,*)
      OPEN(UNIT=NCMECH,FILE=FNMECH,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1010)
      CLOSE(NCMECH)

      WRITE(NCREPT,*)
     +  'Please enter the name of the thermochemical data file:'
      READ(NCINPT,'(A)')FNSPEC
      WRITE(NCREPT,'(A)')FNSPEC
      WRITE(NCREPT,*)
      OPEN(UNIT=NCSPEC,FILE=FNSPEC,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1020)
      CLOSE(NCSPEC)

C     FILENAME ERROR HANDLER
      GOTO 1030
1010  ERRSTR = 'Error opening mechanism file:'
      CALL ERHAND(ERRSTR,IEFATL)
1020  ERRSTR = 'Error opening thermochemical data file:'
      CALL ERHAND(ERRSTR,IEFATL)
1030  CONTINUE
C     END OF FILENAME ERROR HANDLER

C     =========================================================================
      

      RETURN
      END
