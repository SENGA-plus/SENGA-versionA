      SUBROUTINE BCYTXR
 
C     *************************************************************************
C
C     BCYTXR
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
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR MASS FRACTIONS
C     AND THEIR TIME DERIVATIVES
C
C     X-DIRECTION RIGHT-HAND END
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
      INTEGER ISPEC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     EVALUATE AND RETURN STRYXR,DYDTXR
      DO ISPEC = 1,NSPEC

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C           SET MASS FRACTIONS TO CONSTANT (INITIAL) VALUES
            STRYXR(JC,KC,ISPEC) = YRIN(ISPEC)

C           SET MASS FRACTION TIME DERIVATIVES TO ZERO
            DYDTXR(JC,KC,ISPEC) = ZERO

          ENDDO
        ENDDO

      ENDDO

C     =========================================================================


      RETURN
      END
