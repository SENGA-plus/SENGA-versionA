      SUBROUTINE ERHAND(ERRSTR,IERR)
 
C     *************************************************************************
C
C     ERHAND
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     21-SEP-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     HANDLES ERRORS IN PPDIFF
C
C     *************************************************************************


C     PARAMETER
C     =========
      INTEGER NCHERR
      PARAMETER(NCHERR=6)


C     ARGUMENTS
C     =========
      INTEGER IERR
      CHARACTER*50 ERRSTR


C     BEGIN
C     =====
      IF(IERR.EQ.0)THEN
        CONTINUE
      ELSE IF(IERR.EQ.1)THEN
        WRITE(NCHERR,9100)ERRSTR
        WRITE(NCHERR,*)
      ELSE IF(IERR.EQ.2)THEN
        WRITE(NCHERR,9200)ERRSTR
        WRITE(NCHERR,*)
      ELSE IF(IERR.EQ.3)THEN
        WRITE(NCHERR,9300)ERRSTR
        WRITE(NCHERR,*)
      ELSE IF(IERR.EQ.4)THEN
        WRITE(NCHERR,9400)ERRSTR
        WRITE(NCHERR,*)'- stopping PPDIFF...'
        WRITE(NCHERR,*)
        STOP
      ELSE IF(IERR.EQ.5)THEN
        WRITE(NCHERR,9500)ERRSTR
        WRITE(NCHERR,*)'- exiting from PPDIFF'
        WRITE(NCHERR,*)
        STOP
      ENDIF


      RETURN

9100  FORMAT('Information: ',A50)
9200  FORMAT('Warning: ',A50)
9300  FORMAT('Error: ',A50)
9400  FORMAT('Serious error: ',A50)
9500  FORMAT('Fatal error: ',A50)

      END
