      SUBROUTINE FFTNR5(NX,NI,IFORW)

C     *************************************************************************
C
C     FFTNR5
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
C     RADIX 5
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
      DOUBLE PRECISION ONE,FIVE,EIGHT
      PARAMETER(ONE=1.0D0, FIVE=5.0D0, EIGHT=8.0D0)


C     GLOBAL DATA
C     ===========
C     FF5COM-------------------------------------------------------------------

      INTEGER NFTMAX
      PARAMETER(NFTMAX=1024)

      DOUBLE PRECISION COSTM5(NFTMAX),SINTM5(NFTMAX)
      DOUBLE PRECISION COEFR5(4),COEFI5(4)

      INTEGER NXFFTH,MPOWR5,MPHLF5,NI5,INCRM5

      COMMON/FF5COM/COSTM5,SINTM5,
     +              COEFR5,COEFI5,
     +              NXFFTH,MPOWR5,MPHLF5,NI5,INCRM5

C     FF5COM-------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NX,NI,IFORW


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ARGMNT,TWOPIN,TWOPI5,REALNX,RLOG5
      INTEGER IADD,IX,II,IROTN,IROT5


C     BEGIN
C     =====

C     USEFUL CONSTANTS
      REALNX = REAL(NI)
      RLOG5  = LOG(FIVE)
      TWOPIN = REAL(IFORW)*EIGHT*ATAN(ONE)
      TWOPI5 = TWOPIN/FIVE
      TWOPIN = TWOPIN/REALNX
      NXFFTH = NX/5
      MPOWR5 = NINT(LOG(REALNX)/RLOG5)
      IADD = MOD(MPOWR5,2)
      MPHLF5 = MPOWR5/2+IADD
      INCRM5 = NX/NI
      IROTN = MOD(INCRM5,NI)
      IROT5 = MOD(IROTN,5)
      TWOPI5 = TWOPI5*REAL(IROT5)
      NI5 = NI


C     PRECOMPUTE THE RADIX-5 COEFFICIENTS
      DO IX = 1,4
        ARGMNT = TWOPI5*REAL(IX)
        COEFR5(IX) = COS(ARGMNT)
        COEFI5(IX) = SIN(ARGMNT)
      ENDDO

C     PRECOMPUTE THE TRIG FACTORS
      II = 0
      DO IX = 1,4*NXFFTH-3
        ARGMNT = TWOPIN*REAL(II)
        COSTM5(IX) = COS(ARGMNT)
        SINTM5(IX) = SIN(ARGMNT)
        II = II+IROTN
        IF(II.GT.NI)II = II-NI
      ENDDO


      RETURN
      END

