      FUNCTION ESPOVK(WAVENO)
 
C     *************************************************************************
C
C     ESPOVK
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     20-APR-1997:  CREATED
C     24-MAY-2003:  RSC UPDATED FOR SENGA2
C
C     DESCRIPTION
C     -----------
C     EVALUATES THE TURBULENT ENERGY SPECTRUM DIVIDED BY WAVENUMBER
C     AT THE GIVEN WAVENUMBER MAGNITUDE
C
C     *************************************************************************


C     EXTERNAL FUNCTION
C     =================
      DOUBLE PRECISION ESPECT
      EXTERNAL ESPECT


C     PARAMETER
C     =========
      DOUBLE PRECISION ZERO
      PARAMETER(ZERO = 0.0D0)


C     FUNCTION
C     ========
      DOUBLE PRECISION ESPOVK


C     ARGUMENT
C     ========
      DOUBLE PRECISION WAVENO


C     BEGIN
C     =====

C     =========================================================================

      ESPOVK = ZERO
      IF(WAVENO.GT.ZERO)ESPOVK = ESPECT(WAVENO)/WAVENO

C     =========================================================================


      RETURN
      END
