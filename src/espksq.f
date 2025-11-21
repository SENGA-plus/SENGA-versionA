      FUNCTION ESPKSQ(WAVENO)
 
C     *************************************************************************
C
C     ESPKSQ
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
C     EVALUATES THE TURBULENT ENERGY SPECTRUM TIMES WAVENUMBER SQUARED
C     AT THE GIVEN WAVENUMBER MAGNITUDE
C
C     *************************************************************************


C     EXTERNAL FUNCTION
C     =================
      DOUBLE PRECISION ESPECT
      EXTERNAL ESPECT


C     FUNCTION
C     ========
      DOUBLE PRECISION ESPKSQ


C     ARGUMENT
C     ========
      DOUBLE PRECISION WAVENO


C     BEGIN
C     =====

C     =========================================================================

      ESPKSQ = WAVENO*WAVENO*ESPECT(WAVENO) 

C     =========================================================================


      RETURN
      END
