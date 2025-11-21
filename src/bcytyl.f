      SUBROUTINE BCYTYL
 
C     *************************************************************************
C
C     BCYTYL
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
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
C     AND THEIR TIME DERIVATIVES
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
      INTEGER ISPEC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRYYL,DYDTYL
      DO ISPEC = 1,NSPEC

        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C           SET MASS FRACTIONS TO CONSTANT (INITIAL) VALUES
            STRYYL(IC,KC,ISPEC) = YRIN(ISPEC)

C           SET MASS FRACTION TIME DERIVATIVES TO ZERO
            DYDTYL(IC,KC,ISPEC) = ZERO

          ENDDO
        ENDDO

      ENDDO

C     =========================================================================


      RETURN
      END
