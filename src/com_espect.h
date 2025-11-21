C     ESPCOM-------------------------------------------------------------------

C     COMMON DATA FOR ENERGY SPECTRUM FUNCTION
C     ----------------------------------------

C     GENERAL SPECTRUM PARAMETERS
C     ---------------------------
C     NUMBER OF SPECTRUM PARAMETERS
      INTEGER NSPARM
      PARAMETER(NSPARM = 4)

C     SPECTRUM PARAMETERS
      DOUBLE PRECISION SPARAM(NSPARM)


C     DATA SPECIFIC TO THE SPECTRUM FUNCTION
C     --------------------------------------
C     BATCHELOR-TOWNSEND SPECTRUM
      DOUBLE PRECISION CONST0,CK0,OVK0,COVK0


      COMMON/ESPCOM/SPARAM,
     +              CONST0,CK0,OVK0,COVK0

C     ESPCOM-------------------------------------------------------------------
