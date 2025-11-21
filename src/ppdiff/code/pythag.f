      FUNCTION PYTHAG(AVALUE,BVALUE)

C     *************************************************************************
C
C     PYTHAG
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     29-OCT-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     COMPUTES (A^2 + B^2)^1/2 WITHOUT DESTRUCTIVE UNDERFLOW OR OVERFLOW
C
C     BASED ON NUMERICAL RECIPES FUNCTION PYTHAG
C
C     *************************************************************************

C     PARAMETERS
C     ==========
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)


C     FUNCTION VALUE
C     ==============
      DOUBLE PRECISION PYTHAG


C     ARGUMENTS
C     =========
      DOUBLE PRECISION AVALUE,BVALUE


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ABSA,ABSB


C     BEGIN
C     =====
      ABSA = ABS(AVALUE)
      ABSB = ABS(BVALUE)
      IF(ABSA.GT.ABSB)THEN
        PYTHAG = ABSA*SQRT(ONE+(ABSB/ABSA)**2)
      ELSE
        IF(ABSB.EQ.ZERO)THEN
          PYTHAG = ZERO
        ELSE
          PYTHAG = ABSB*SQRT(ONE+(ABSA/ABSB)**2)
        ENDIF
      ENDIF


      RETURN
      END
