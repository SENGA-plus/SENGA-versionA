      SUBROUTINE P_SEND(PLOCAL,NPARAY,NCOUNT,IDPROC,ITAG)

C     *************************************************************************
C
C     P_SEND
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     13-APR-1997:  CREATED
C     11-MAY-2003:  RSC MODIFIED FOR SENGA2
C     11-AUG-2009:  RSC STUB VERSION
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     PARALLEL SEND ROUTINE
C     STUB FOR NON-PARALLEL SYSTEMS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'mpif.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NPARAY,NCOUNT,IDPROC,ITAG
      DOUBLE PRECISION PLOCAL(NPARAY)


C     BEGIN
C     =====

C     =========================================================================

C     NO ACTION REQUIRED

C     =========================================================================


      RETURN
      END
