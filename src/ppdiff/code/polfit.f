      SUBROUTINE POLFIT(XPOINT,YPOINT,NPOINT,NFITMX,
     +                  ACOEFS,NCOEFS,NFCOMX,CHISQR)
 
C     *************************************************************************
C
C     POLFIT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     07-OCT-2012:  CREATED
C      
C     DESCRIPTION
C     -----------
C     FITS A POLYNOMIAL TO THE MOLECULAR TRANSPORT DATA FOR PPDIFF
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ZERO, ONE
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)


C     ARGUMENTS
C     =========
      DOUBLE PRECISION XPOINT(NFITMX),YPOINT(NFITMX)
      DOUBLE PRECISION ACOEFS(NFCOMX)
      DOUBLE PRECISION CHISQR
      INTEGER NPOINT,NFITMX,NCOEFS,NFCOMX


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION UMATRX(NFITMX,NFCOMX)
      DOUBLE PRECISION VMATRX(NFCOMX,NFCOMX)
      DOUBLE PRECISION WVECTR(NFCOMX)
      DOUBLE PRECISION SIGDEV(NFITMX)
      INTEGER ICOEFF


C     BEGIN
C     =====

C     =========================================================================

C     INITIALISE
C     -----------

C     PRESET THE "MEASUREMENT ERRORS"
      DO IPOINT = 1, NPOINT
        SIGDEV(IPOINT) = ONE
      ENDDO

C     INITIALISE THE COEFFICIENTS
      DO ICOEFF = 1, NCOEFS
        ACOEFS(ICOEFF) = ZERO
      ENDDO

C     INITIALISE THE GOODNESS-OF-FIT PARAMETER
      CHISQR = ZERO

C     =========================================================================

C     FIT USING SINGULAR VALUE DECOMPOSITION
C     --------------------------------------
      CALL SVDFIT(XPOINT,YPOINT,SIGDEV,NPOINT,ACOEFS,NCOEFS,
     +            UMATRX,VMATRX,WVECTR,NFITMX,NFCOMX,CHISQR)

C     =========================================================================

CC     DIAGNOSTICS
C      WRITE(6,*)'POLFIT:',NCOEFS
C      DO ICOEFF = 1, NCOEFS
C        WRITE(6,'(I5,1PE12.4)')ICOEFF,ACOEFS(ICOEFF)
C      ENDDO


      RETURN
      END
