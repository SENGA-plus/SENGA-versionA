      SUBROUTINE P_SUMM(ELEMNT,ANSWER)
 
C     *************************************************************************
C
C     P_SUMM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     12-MAR-2003:  CREATED
C
C     DESCRIPTION
C     -----------
C     SUMS DATA FROM ALL PROCESSORS
C     AND DISTRIBUTES THE ANSWER TO ALL PROCESSORS
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
      DOUBLE PRECISION ELEMNT,ANSWER


C     LOCAL DATA
C     ==========
      INTEGER IERROR


C     BEGIN
C     =====

C     =========================================================================

      CALL MPI_ALLREDUCE(ELEMNT,ANSWER,1,MPI_DOUBLE_PRECISION,
     +                   MPI_SUM,MPI_COMM_WORLD,IERROR)

C     =========================================================================


      RETURN
      END
