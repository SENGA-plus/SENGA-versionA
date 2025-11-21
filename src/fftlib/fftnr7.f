      SUBROUTINE FFTNR7(NX,NI,IFORW)

C     *************************************************************************
C
C     FFTNR7
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-MAY-1999:  CREATED
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT AN FFT IN 1D
C     USING TEMPERTON SELF-SORTING PRIME FACTOR ALGORITHM
C     RADIX 7
C     WITH ROTATIONS
C     INITIALISATION ONLY
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 13,3,676-686, 1992
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ONE,SEVEN,EIGHT
      PARAMETER(ONE=1.0D0, SEVEN=7.0D0, EIGHT=8.0D0)


C     GLOBAL DATA
C     ===========
C     FF7COM-------------------------------------------------------------------

      INTEGER NFTMAX
      PARAMETER(NFTMAX=1024)

      DOUBLE PRECISION COSTM7(NFTMAX),SINTM7(NFTMAX)
      DOUBLE PRECISION COEFR7(6),COEFI7(6)

      INTEGER NXSVTH,MPOWR7,MPHLF7,NI7,INCRM7

      COMMON/FF7COM/COSTM7,SINTM7,
     +              COEFR7,COEFI7,
     +              NXSVTH,MPOWR7,MPHLF7,NI7,INCRM7

C     FF7COM-------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NX,NI,IFORW


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ARGMNT,TWOPIN,TWOPI7,REALNX,RLOG7
      INTEGER IADD,IX,II,IROTN,IROT7


C     BEGIN
C     =====

C     USEFUL CONSTANTS
      REALNX = REAL(NI)
      RLOG7  = LOG(SEVEN)
      TWOPIN = REAL(IFORW)*EIGHT*ATAN(ONE)
      TWOPI7 = TWOPIN/SEVEN
      TWOPIN = TWOPIN/REALNX
      NXSVTH = NX/7
      MPOWR7 = NINT(LOG(REALNX)/RLOG7)
      IADD = MOD(MPOWR7,2)
      MPHLF7 = MPOWR7/2+IADD
      INCRM7 = NX/NI
      IROTN = MOD(INCRM7,NI)
      IROT7 = MOD(IROTN,7)
      TWOPI7 = TWOPI7*REAL(IROT7)
      NI7 = NI

C     PRECOMPUTE THE RADIX-7 COEFFICIENTS
      DO IX = 1,6
        ARGMNT = TWOPI7*REAL(IX)
        COEFR7(IX) = COS(ARGMNT)
        COEFI7(IX) = SIN(ARGMNT)
      ENDDO

C     PRECOMPUTE THE TRIG FACTORS
      II = 0
      DO IX = 1,6*NXSVTH-5
        ARGMNT = TWOPIN*REAL(II)
        COSTM7(IX) = COS(ARGMNT)
        SINTM7(IX) = SIN(ARGMNT)
        II = II+IROTN
        IF(II.GT.NI)II = II-NI
      ENDDO


      RETURN
      END

