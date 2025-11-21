      SUBROUTINE P_SYNC
 
C     *************************************************************************
C
C     P_SYNC
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     24-MAR-2003:  CREATED
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT SYNCHRONISATION OF PARALLEL PROCESSORS USING MPI
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
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

C     =========================================================================


      RETURN
      END
