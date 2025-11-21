      PROGRAM PPDIFF
 
C     *************************************************************************
C
C     PPDIFF
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-SEP-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     PREPROCESSES MOLECULAR TRANSPORT DATA FOR SENGA
C
C     *************************************************************************


C     GLOBAL DATA
C     +++++++++++
C     -------------------------------------------------------------------------
      INCLUDE 'com_diffin.h'
      INCLUDE 'com_ppdcom.h'
C     -------------------------------------------------------------------------


C     BEGIN
C     =====

C     =========================================================================

C     ANNOUNCEMENT
C     ------------
      WRITE(NCREPT,*)'========================================'
      WRITE(NCREPT,*)
      WRITE(NCREPT,*)'PPDIFF: molecular transport preprocessor'
      WRITE(NCREPT,*)
      WRITE(NCREPT,*)'========================================'
      WRITE(NCREPT,*)

C     =========================================================================

C     GET THE FILENAMES AND CHECK THE FILES
C     -------------------------------------
      CALL FNDCHK

C     =========================================================================

C     PROCESS THE SPECIES DATA
C     ------------------------
      CALL PDSPEC

C     =========================================================================

C     PROCESS THE THERMOCHEMICAL DATA
C     -------------------------------
      CALL PDTHMO

C     =========================================================================

C     PROCESS THE DIFFUSION DATA
C     --------------------------
      CALL PDDIFF

C     =========================================================================

C     PROCESS THE COLLISION INTEGRAL DATA
C     -----------------------------------
      CALL PDCOLL

C     =========================================================================

C     INITIALISE THE LISTS AND TABLES
C     -------------------------------
      CALL MOLTAB

C     =========================================================================

C     EVALUATE THE VISCOSITY DATA
C     ---------------------------
      CALL VISDAT

C     =========================================================================

C     EVALUATE THE BINARY DIFFUSION COEFFICIENT DATA
C     ----------------------------------------------
      CALL DIFDAT

C     =========================================================================

C     EVALUATE THE THERMAL CONDUCTIVITY DATA
C     --------------------------------------
      CALL CONDAT

C     =========================================================================

C     EVALUATE THE THERMAL DIFFUSION RATIO DATA
C     -----------------------------------------
      CALL TDRDAT

C     =========================================================================

C     WRITE OUT THE LISTS
C     -------------------
      CALL OUTMOL

C     =========================================================================


      STOP
      END
