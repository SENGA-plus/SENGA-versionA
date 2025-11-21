      SUBROUTINE BCUTXR
 
C     *************************************************************************
C
C     BCUTXR
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
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
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
CKA   FIX INFLOW BC
CKA      DOUBLE PRECISION BTIME
      DOUBLE PRECISION FORNOW,ARGMNT
      INTEGER JC,KC


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
CKA   FIX INFLOW BC
CKA      BTIME = ETIME + RKTIM(IRKSTP)

C     =========================================================================

C     CONSTANT U-VELOCITY
C     PARAMETER I1=1, R1=U-VELOCITY
      IF(NXRPRM(1).EQ.1)THEN

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXR(JC,KC) = RXRPRM(1)
            STRVXR(JC,KC) = ZERO
            STRWXR(JC,KC) = ZERO

            DUDTXR(JC,KC) = ZERO
            DVDTXR(JC,KC) = ZERO
            DWDTXR(JC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     SINUSOIDAL U-VELOCITY
C     PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD
      IF(NXRPRM(1).EQ.2)THEN

        FORNOW = TWO*PI/RXRPRM(2)
        ARGMNT = FORNOW*BTIME

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXR(JC,KC) = RXRPRM(1)*SIN(ARGMNT)
            STRVXR(JC,KC) = ZERO
            STRWXR(JC,KC) = ZERO

            DUDTXR(JC,KC) = FORNOW*RXRPRM(1)*COS(ARGMNT)
            DVDTXR(JC,KC) = ZERO
            DWDTXR(JC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
