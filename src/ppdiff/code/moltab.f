      SUBROUTINE MOLTAB
 
C     *************************************************************************
C
C     MOLTAB
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     07-SEP-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     INITIALISES LISTS AND TABLES FOR PPDIFF
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_diffin.h'
      INCLUDE 'com_ppdcom.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER ISPEC,JSPEC
      INTEGER ICOEF


C     BEGIN
C     =====

C     =========================================================================

C     CONSTANTS
C     ---------
      PI = TWO*TWO*ATAN(ONE)
      TWOPI = TWO*PI
      ROOTPI = SQRT(PI)
      ROOTWO = SQRT(TWO)

C     =========================================================================

C     REFERENCE TEMPERATURE
C     ---------------------
      TREFGB = 3.0D2

C     =========================================================================

C     INITIALISE THE BINARY DIFFUSION COEFFICIENT TABLE
C     -------------------------------------------------
      DO JSPEC = 1,NSPCMX
        DO ISPEC = 1,NSPCMX
          DO ICOEF = 1,NDCFMX

            DIFFCO(ICOEF,ISPEC,JSPEC) = ZERO

          ENDDO
        ENDDO
      ENDDO

C     INITIALISE THE VISCOSITY COEFFICIENT TABLE
C     ------------------------------------------
      DO ISPEC = 1,NSPCMX
        DO ICOEF = 1,NVCFMX

          VISCCO(ICOEF,ISPEC) = ZERO

        ENDDO
      ENDDO

C     INITIALISE THE CONDUCTIVITY COEFFICIENT TABLE
C     ---------------------------------------------
      DO ISPEC = 1,NSPCMX
        DO ICOEF = 1,NCCFMX

          CONDCO(ICOEF,ISPEC) = ZERO

        ENDDO
      ENDDO

C     =========================================================================
      

      RETURN
      END
