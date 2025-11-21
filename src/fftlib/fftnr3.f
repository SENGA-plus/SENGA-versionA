      SUBROUTINE FFTNR3(NX,NI,IFORW)

C     *************************************************************************
C
C     FFTNR3
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
C     RADIX 3
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
      DOUBLE PRECISION ONE,THREE,EIGHT
      PARAMETER(ONE=1.0D0, THREE=3.0D0, EIGHT=8.0D0)


C     GLOBAL DATA
C     ===========
C     FF3COM-------------------------------------------------------------------

      INTEGER NFTMAX
      PARAMETER(NFTMAX=1024)

      DOUBLE PRECISION COSTM3(NFTMAX),SINTM3(NFTMAX)
      DOUBLE PRECISION COEFR3(2),COEFI3(2)

      INTEGER NXTHRD,MPOWR3,MPHLF3,NI3,INCRM3

      COMMON/FF3COM/COSTM3,SINTM3,
     +              COEFR3,COEFI3,
     +              NXTHRD,MPOWR3,MPHLF3,NI3,INCRM3

C     FF3COM-------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NX,NI,IFORW


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ARGMNT,TWOPIN,TWOPI3,REALNX,RLOG3
      INTEGER IX,IROTN,IROT3,II,IADD


C     BEGIN
C     =====

C     USEFUL CONSTANTS
      REALNX = REAL(NI)
      RLOG3  = LOG(THREE)
      TWOPIN = REAL(IFORW)*EIGHT*ATAN(ONE)
      TWOPI3 = TWOPIN/THREE
      TWOPIN = TWOPIN/REALNX
      NXTHRD = NX/3
      MPOWR3 = NINT(LOG(REALNX)/RLOG3)
      IADD = MOD(MPOWR3,2)
      MPHLF3 = MPOWR3/2+IADD
      INCRM3 = NX/NI
      IROTN = MOD(INCRM3,NI)
      IROT3 = MOD(IROTN,3)
      TWOPI3 = TWOPI3*REAL(IROT3)
      NI3 = NI

C     PRECOMPUTE THE RADIX-3 COEFFICIENTS
      DO IX = 1,2
        ARGMNT = TWOPI3*REAL(IX)
        COEFR3(IX) = COS(ARGMNT)
        COEFI3(IX) = SIN(ARGMNT)
      ENDDO

C     PRECOMPUTE THE TRIG FACTORS
      II = 0
      DO IX = 1,2*NXTHRD-1
        ARGMNT = TWOPIN*REAL(II)
        COSTM3(IX) = COS(ARGMNT)
        SINTM3(IX) = SIN(ARGMNT)
        II = II+IROTN
        IF(II.GT.NI)II = II-NI
      ENDDO


      RETURN
      END


