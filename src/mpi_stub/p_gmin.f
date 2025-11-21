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
C     11-AUG-2009:  RSC STUB VERSION
C
C     DESCRIPTION
C     -----------
C     FINDS THE MIN OF DATA FROM ALL PROCESSORS
C     AND DISTRIBUTES THE ANSWER TO ALL PROCESSORS
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
      DOUBLE PRECISION ELEMNT,ANSWER


C     BEGIN
C     =====

C     =========================================================================

      ANSWER = ELEMNT

C     =========================================================================


      RETURN
      END
