      FUNCTION ERFUNC(ARGMNT)
 
C     *************************************************************************
C
C     ERFUNC
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     09-JAN-2005:  CREATED
C
C     DESCRIPTION
C     -----------
C     COMPUTES THE ERROR FUNCTION
C     REF: ABRAMOWITZ AND STEGUN p299 sec 7.1.26
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ZERO,ONE,HALF
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, HALF=5.0D-1)

      INTEGER NCOEFF,NCOFM1
      PARAMETER(NCOEFF = 5, NCOFM1 = NCOEFF-1)


C     FUNCTION
C     ========
      DOUBLE PRECISION ERFUNC


C     ARGUMENT
C     ========
      DOUBLE PRECISION ARGMNT


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ECOEFF(NCOEFF)
      DOUBLE PRECISION PCOEFF
      DOUBLE PRECISION ETOTAL,ZVALUE,TVALUE
      INTEGER ICOEFF


C     BEGIN
C     =====

C     =========================================================================

C     SET THE COEFFICIENTS
      PCOEFF = 0.3275911D0
      ECOEFF(1) = 0.254829592D0
      ECOEFF(2) =-0.284496736D0
      ECOEFF(3) = 1.421413741D0
      ECOEFF(4) =-1.453152027D0
      ECOEFF(5) = 1.061405429D0

C     EVALUATE ERROR FUNCTION
      ZVALUE = ABS(ARGMNT)
      TVALUE = ONE/(ONE+PCOEFF*ZVALUE)

      ETOTAL = ECOEFF(NCOEFF)
      DO ICOEFF = NCOFM1,1,-1
        ETOTAL = ECOEFF(ICOEFF) + ETOTAL*TVALUE
      ENDDO
      ETOTAL = ETOTAL*TVALUE

      ERFUNC = ONE - ETOTAL*EXP(-ZVALUE*ZVALUE)
      IF(ARGMNT.LT.ZERO)ERFUNC = -ERFUNC

C     =========================================================================


      RETURN
      END
