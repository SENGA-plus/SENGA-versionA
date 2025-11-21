      SUBROUTINE BCDTZR
 
C     *************************************************************************
C
C     BCDTZR
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     26-OCT-2013:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR DENSITY
C     AND ITS TIME DERIVATIVE
C
C     Z-DIRECTION RIGHT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER IC,JC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRDZR,DDDTZR

      DO JC = JSTAL,JSTOL
        DO IC = ISTAL,ISTOL

C         SET DENSITY TO CONSTANT (INITIAL) VALUE
          STRDZR(IC,JC) = DRIN

C         SET DENSITY TIME DERIVATIVE TO ZERO
          DDDTZR(IC,JC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
