      SUBROUTINE BCTTXL
 
C     *************************************************************************
C
C     BCTTXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-DEC-2003:  CREATED
C     09-MAY-2015:  RSC MODIFIED FOR ISOTHERMAL WALL
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR TEMPERATURE
C     AND ITS TIME DERIVATIVE
C
C     X-DIRECTION LEFT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER JC,KC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRTXL,DTDTXL
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL

          STRTXL(JC,KC) = TRIN

          DTDTXL(JC,KC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================

C     ISOTHERMAL WALL
      IF(NSBCXL.EQ.NSBCW2)THEN

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRTXL(JC,KC) = RXLPRM(1)

            DTDTXL(JC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
