      SUBROUTINE P_BCST(PLOCAL,NPARAY,NCOUNT,IDPROC)

C     *************************************************************************
C
C     P_BCST
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-DEC-2003:  CREATED
C     11-AUG-2009: RSC STUB VERSION
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     PARALLEL BROADCAST ROUTINE
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
      INTEGER NPARAY,NCOUNT,IDPROC
      DOUBLE PRECISION PLOCAL(NPARAY)


C     BEGIN
C     =====

C     =========================================================================

C     NO ACTION REQUIRED

C     =========================================================================


      RETURN
      END
