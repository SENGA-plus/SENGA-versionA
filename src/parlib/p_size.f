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
C
C     DESCRIPTION
C     -----------
C     GETS THE NUMBER OF PROCESSORS NPROC
C     AND THE LOCAL PROCESSOR ID IPROC, IPROC = 0,NPROC-1
C     USING MPI
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


C     LOCAL DATA
C     ==========
      INTEGER IERROR


C     BEGIN
C     =====

C     =========================================================================

C     GET NUMBER OF PROCESSORS
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERROR)

C     =========================================================================

C     GET ID OF THIS PROCESSOR
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IPROC,IERROR)

C     =========================================================================


      RETURN
      END
