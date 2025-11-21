      SUBROUTINE FINISH
 
C     *************************************************************************
C
C     FINISH
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     09-MAR-2003:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     TERMINATES THE PROGRAM
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
#ifdef HDF5
      USE hdf5io
#endif
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------

C     BEGIN
C     =====

C     =========================================================================

C     FINAL REPORT
C     ------------
      IF(IPROC.EQ.0)THEN

        OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='OLD',FORM='FORMATTED')

C       GO TO EOF
1000    CONTINUE
          READ(NCREPT,'(A)',END=1010)
          GOTO 1000
C       END OF LOOP 1000

1010    BACKSPACE(NCREPT)
        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Normal termination of program'
        WRITE(NCREPT,*)'-----------------------------'
        CLOSE(NCREPT)

      ENDIF

C     =========================================================================
#ifdef HDF5
      CALL HDF5_CLOSE
#endif
C     TERMINATE PARALLEL MESSAGE PASSING
      CALL P_STOP

C     =========================================================================

      RETURN
      END
