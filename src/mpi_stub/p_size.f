      SUBROUTINE P_SIZE(NPROC,IPROC)
 
C     *************************************************************************
C
C     P_SIZE
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     24-MAR-2003:  CREATED
C     11-AUG-2009:  RSC STUB VERSION
C
C     DESCRIPTION
C     -----------
C     GETS THE NUMBER OF PROCESSORS NPROC
C     AND THE LOCAL PROCESSOR ID IPROC, IPROC = 0,NPROC-1
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
      INTEGER IPROC,NPROC


C     BEGIN
C     =====

C     =========================================================================

      NPROC = 1
      IPROC = 0

C     =========================================================================


      RETURN
      END
