      SUBROUTINE BCUTYL
 
C     *************************************************************************
C
C     BCUTYL
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
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
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


C     BEGIN
C     =====

C     =========================================================================

C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)

C     =========================================================================

C     CONSTANT V-VELOCITY
C     PARAMETER I1=1, R1=V-VELOCITY
      IF(NYLPRM(1).EQ.1)THEN

        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRUYL(IC,KC) = ZERO
            STRVYL(IC,KC) = RYLPRM(1)
            STRWYL(IC,KC) = ZERO

            DUDTYL(IC,KC) = ZERO
            DVDTYL(IC,KC) = ZERO
            DWDTYL(IC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
