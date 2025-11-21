      SUBROUTINE BCTTZL
 
C     *************************************************************************
C
C     BCTTZL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     26-OCT-2013:  CREATED
C     09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
C     AND ITS TIME DERIVATIVE
C
C     Z-DIRECTION LEFT-HAND END
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

C     EVALUATE AND RETURN STRTZL,DTDTZL
      DO JC = JSTAL,JSTOL
        DO IC = ISTAL,ISTOL

          STRTZL(IC,JC) = TRIN

          DTDTZL(IC,JC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================

C     ISOTHERMAL WALL
      IF(NSBCZL.EQ.NSBCW2)THEN

        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STRTZL(IC,JC) = RZLPRM(1)

            DTDTZL(IC,JC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
