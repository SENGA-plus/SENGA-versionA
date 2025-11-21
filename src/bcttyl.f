      SUBROUTINE BCTTYL
 
C     *************************************************************************
C
C     BCTTYL
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
C     Y-DIRECTION LEFT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER IC,KC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRTYL,DTDTYL
      DO KC = KSTAL,KSTOL
        DO IC = ISTAL,ISTOL

          STRTYL(IC,KC) = TRIN

          DTDTYL(IC,KC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================

C     ISOTHERMAL WALL
      IF(NSBCYL.EQ.NSBCW2)THEN

        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRTYL(IC,KC) = RYLPRM(1)

            DTDTYL(IC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
