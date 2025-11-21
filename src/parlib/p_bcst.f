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
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     PARALLEL BROADCAST ROUTINE
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


C     LOCAL DATA
C     ==========
      INTEGER IERROR


C     BEGIN
C     =====

C     =========================================================================

      CALL MPI_BCAST(PLOCAL,NCOUNT,MPI_DOUBLE_PRECISION,
     +               IDPROC,MPI_COMM_WORLD,IERROR)

C     =========================================================================


      RETURN
      END
