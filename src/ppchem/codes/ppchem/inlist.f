      SUBROUTINE INLIST
 
C     *************************************************************************
C
C     INLIST
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     INITIALISES LISTS AND TABLES FOR PPCHEM
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_chemin.h'
      INCLUDE 'com_ppccom.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER ISPEC
      INTEGER ISTEP
      INTEGER IBODY


C     BEGIN
C     =====

C     =========================================================================

C     INITIALISE THE STEP-SPECIES TABLES
C     ----------------------------------
      DO ISTEP = 1,NSTPMX
        DO ISPEC = 1,NSPCMX
          NRTABL(ISPEC,ISTEP) = 0
          NPTABL(ISPEC,ISTEP) = 0
        ENDDO
      ENDDO

C     INITIALISE THE STEP-COEFFICIENT TABLES
C     --------------------------------------
      DO ISTEP = 1,NSTPMX
        DO ISPEC = 1,NSPCMX
          CRTABL(ISPEC,ISTEP) = ZERO
          CPTABL(ISPEC,ISTEP) = ZERO
        ENDDO
      ENDDO

C     INITIALISE THE THIRD BODY STEP LIST
C     -----------------------------------
      DO ISTEP = 1,NSTPMX
        MBLIST(ISTEP) = 0
      ENDDO

C     INITIALISE THE GIBBS FUNCTION STEP LIST
C     ---------------------------------------
      NGIBB = 0
      DO ISTEP = 1,NSTPMX
        MGLIST(ISTEP) = 0
      ENDDO

C     INITIALISE THE LINDEMANN FORM STEP LIST
C     ---------------------------------------
      NLIND = 0
      DO ISTEP = 1,NSTPMX
        MLLIST(ISTEP) = 0
      ENDDO

C     INITIALISE THE TROE FORM STEP LIST
C     ----------------------------------
      NTROE = 0
      DO ISTEP = 1,NSTPMX
        MTLIST(ISTEP) = 0
      ENDDO

C     INITIALISE THE SRI FORM STEP LIST
C     ---------------------------------
      NSRIF = 0
      DO ISTEP = 1,NSTPMX
        MSLIST(ISTEP) = 0
      ENDDO

C     INITIALISE THE THIRD-BODY EFFICIENCY TABLE
C     ------------------------------------------
      DO IBODY = 1,NBDYMX
        DO ISPEC = 1,NSPCMX
          EFFY3B(ISPEC,IBODY) = ONE
        ENDDO
      ENDDO

C     =========================================================================
      

      RETURN
      END
