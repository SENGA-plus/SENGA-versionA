      SUBROUTINE BCUTYR
 
C     *************************************************************************
C
C     BCUTYR
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
C     Y-DIRECTION RIGHT-HAND END
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
      IF(NYRPRM(1).EQ.1)THEN

        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

            STRUYR(IC,KC) = ZERO
            STRVYR(IC,KC) = RYRPRM(1)
            STRWYR(IC,KC) = ZERO

            DUDTYR(IC,KC) = ZERO
            DVDTYR(IC,KC) = ZERO
            DWDTYR(IC,KC) = ZERO

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
