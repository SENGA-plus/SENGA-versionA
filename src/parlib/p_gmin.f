      SUBROUTINE P_GMIN(ELEMNT,ANSWER)
 
C     *************************************************************************
C
C     P_GMIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     05-JUN-2004:  CREATED
C
C     DESCRIPTION
C     -----------
C     FINDS THE MIN OF DATA FROM ALL PROCESSORS
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
     +                   MPI_MIN,MPI_COMM_WORLD,IERROR)

C     =========================================================================


      RETURN
      END
