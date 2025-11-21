      SUBROUTINE P_STOP
 
C     *************************************************************************
C
C     P_STOP
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     24-MAR-2003:  CREATED

C     DESCRIPTION
C     -----------
C     CARRIES OUT TERMINATION OF PARALLEL TASKS USING MPI
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'mpif.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER IERROR


C     BEGIN
C     =====

C     =========================================================================

C     SYNCHRONISE
      CALL MPI_FINALIZE(IERROR)

C     =========================================================================


      RETURN
      END
