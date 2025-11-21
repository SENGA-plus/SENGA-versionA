      SUBROUTINE P_INIT
 
C     *************************************************************************
C
C     P_INIT
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
C     CARRIES OUT INITIALISATION OF PARALLEL DATA USING MPI
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

C     INITIALISE MPI
      CALL MPI_INIT(IERROR)

C     =========================================================================


      RETURN
      END
