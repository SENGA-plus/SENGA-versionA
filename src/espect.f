      FUNCTION ESPECT(WAVENO)
 
C     *************************************************************************
C
C     ESPECT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  DEPARTMENT OF MECHANICAL ENGINEERING, UMIST
C
C     CHANGE RECORD
C     -------------
C     28-JUN-1994:  CREATED
C     13-JUL-2001:  RSC NEW VERSION FOR BATCHELOR-TOWNSEND SPECTRUM
C     24-MAY-2003:  RSC UPDATED FOR SENGA2
C
C     DESCRIPTION
C     -----------
C     EVALUATES THE TURBULENT ENERGY SPECTRUM
C     AT THE GIVEN WAVENUMBER MAGNITUDE
C
C     REFERENCES
C     ----------
C     1) BATCHELOR, G.K., TOWNSEND, A.A.: JFM 88(4), 685-709, 1948.
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_espect.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
      DOUBLE PRECISION TWO
      PARAMETER(TWO = 2.0D0)


C     FUNCTION
C     ========
      DOUBLE PRECISION ESPECT


C     ARGUMENT
C     ========
      DOUBLE PRECISION WAVENO


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION WAVRAT


C     BEGIN
C     =====

C     =========================================================================

      WAVRAT = WAVENO*OVK0
      WAVRAT = WAVRAT*WAVRAT
      ESPECT = COVK0*WAVRAT*WAVRAT*EXP(-TWO*WAVRAT)

C     =========================================================================

      
      RETURN
      END
