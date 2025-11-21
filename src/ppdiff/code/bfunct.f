      SUBROUTINE BFUNCT(XPOINT,ABFUNC,MCOEFF)

C     *************************************************************************
C
C     BFUNCT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     13-OCT-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     EVALUATES THE FITTING BASIS FUNCTIONS AT A POINT 
C
C     *************************************************************************


C     PARAMETER
C     =========
      DOUBLE PRECISION ONE
      PARAMETER(ONE = 1.0D0)


C     ARGUMENTS
C     =========
      DOUBLE PRECISION ABFUNC(MCOEFF)
      DOUBLE PRECISION XPOINT
      INTEGER MCOEFF


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FORNOW
      INTEGER IC


C     BEGIN
C     =====
C     BASIS FUNCTIONS ARE MONOMIALS
      FORNOW = ONE
      ABFUNC(1) = FORNOW
      DO IC = 2, MCOEFF
        FORNOW = FORNOW*XPOINT
        ABFUNC(IC) = FORNOW
      ENDDO


      RETURN
      END
