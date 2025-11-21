      SUBROUTINE FNDCHK
 
C     *************************************************************************
C
C     FNDCHK
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-SEP-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     CHECKS FILES FOR PPDIFF
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_diffin.h'
      INCLUDE 'com_ppdcom.h'
C     -------------------------------------------------------------------------


C     BEGIN
C     =====

C     =========================================================================

C     GET THE FILENAMES AND CHECK THE FILES
C     -------------------------------------
      WRITE(NCREPT,*)'Please enter the name of the species data file:'
      READ(NCINPT,'(A)')FNSPEC
      WRITE(NCREPT,'(A)')FNSPEC
      WRITE(NCREPT,*)
      OPEN(UNIT=NCSPEC,FILE=FNSPEC,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1000)
      CLOSE(NCSPEC)

      WRITE(NCREPT,*)
     +      'Please enter the name of the thermochemical data file:'
      READ(NCINPT,'(A)')FNTHMO
      WRITE(NCREPT,'(A)')FNTHMO
      WRITE(NCREPT,*)
      OPEN(UNIT=NCTHMO,FILE=FNTHMO,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1010)
      CLOSE(NCTHMO)

      WRITE(NCREPT,*)
     +      'Please enter the name of the diffusion data file:'
      READ(NCINPT,'(A)')FNDIFF
      WRITE(NCREPT,'(A)')FNDIFF
      WRITE(NCREPT,*)
      OPEN(UNIT=NCDIFF,FILE=FNDIFF,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1020)
      CLOSE(NCDIFF)

      WRITE(NCREPT,*)
     +      'Please enter the name of the collision integral data file:'
      READ(NCINPT,'(A)')FNCOLL
      WRITE(NCREPT,'(A)')FNCOLL
      WRITE(NCREPT,*)
      OPEN(UNIT=NCCOLL,FILE=FNCOLL,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1030)
      CLOSE(NCCOLL)

      WRITE(NCREPT,*)'Please enter the name of the output file:'
      READ(NCINPT,'(A)')FNOUTP
      WRITE(NCREPT,'(A)')FNOUTP
      WRITE(NCREPT,*)
      OPEN(UNIT=NCOUTP,FILE=FNOUTP,STATUS='UNKNOWN',FORM='FORMATTED',
     +     ERR=1040)

C     FILENAME ERROR HANDLER
      GOTO 1100
1000  WRITE(NCREPT,*)'Error opening species data file:'
      WRITE(NCREPT,*)FNSPEC
      WRITE(NCREPT,*)'- exiting from PPDIFF'
      WRITE(NCREPT,*)
      STOP
1010  WRITE(NCREPT,*)'Error opening thermochemical data file:'
      WRITE(NCREPT,*)FNTHMO
      WRITE(NCREPT,*)'- exiting from PPDIFF'
      WRITE(NCREPT,*)
      STOP
1020  WRITE(NCREPT,*)'Error opening diffusion data file:'
      WRITE(NCREPT,*)FNDIFF
      WRITE(NCREPT,*)'- exiting from PPDIFF'
      WRITE(NCREPT,*)
      STOP
1030  WRITE(NCREPT,*)'Error opening collision integral data file:'
      WRITE(NCREPT,*)FNCOLL
      WRITE(NCREPT,*)'- exiting from PPDIFF'
      WRITE(NCREPT,*)
      STOP
1040  WRITE(NCREPT,*)'Error opening output file:'
      WRITE(NCREPT,*)FNOUTP
      WRITE(NCREPT,*)'- exiting from PPDIFF'
      WRITE(NCREPT,*)
      STOP
1100  CONTINUE
C     END OF FILENAME ERROR HANDLER

C     =========================================================================


      RETURN
      END
