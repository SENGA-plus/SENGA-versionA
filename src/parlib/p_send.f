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
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     PARALLEL SEND ROUTINE
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


C     LOCAL DATA
C     ==========
      INTEGER IERROR


C     BEGIN
C     =====

C     =========================================================================

      CALL MPI_SSEND(PLOCAL,NCOUNT,MPI_DOUBLE_PRECISION,
     +               IDPROC,ITAG,MPI_COMM_WORLD,IERROR)

C     =========================================================================


      RETURN
      END
