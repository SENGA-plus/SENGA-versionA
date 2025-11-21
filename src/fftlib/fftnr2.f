      SUBROUTINE FFTNR2(NX,NI,IFORW)

C     *************************************************************************
C
C     FFTNR2
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     21-MAY-1999:  CREATED
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT AN FFT IN 1D
C     USING TEMPERTON SELF-SORTING PRIME FACTOR ALGORITHM
C     RADIX 2
C     WITH ROTATIONS
C     INITIALISATION ONLY
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 12,4,808-823, 1991
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      DOUBLE PRECISION ONE,TWO,EIGHT
      PARAMETER(ONE=1.0D0, TWO=2.0D0, EIGHT=8.0D0)


C     GLOBAL DATA
C     ===========
C     FF2COM-------------------------------------------------------------------

      INTEGER NFTMAX
      PARAMETER(NFTMAX=1024)

      DOUBLE PRECISION COSTM2(NFTMAX),SINTM2(NFTMAX)

      INTEGER NXHALF,MPOWR2,MPHLF2,NI2,INCRM2

      COMMON/FF2COM/COSTM2,SINTM2,
     +              NXHALF,MPOWR2,MPHLF2,NI2,INCRM2

C     FF2COM-------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NX,NI,IFORW


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION ARGMNT,TWOPIN,REALNX,RLOG2
      INTEGER IADD,IX,KK,IROTN


C     BEGIN
C     =====

C     USEFUL CONSTANTS
      REALNX = REAL(NI)
      RLOG2  = LOG(TWO)
      TWOPIN = REAL(IFORW)*EIGHT*ATAN(ONE)/REALNX
      NXHALF = NX/2
      MPOWR2 = NINT(LOG(REALNX)/RLOG2)
      IADD = MOD(MPOWR2,2)
      MPHLF2 = MPOWR2/2+IADD
      INCRM2 = NX/NI
      NI2 = NI

C     PRECOMPUTE THE TRIG FACTORS
      IROTN = MOD(INCRM2,NI)
      KK = 0
      DO IX = 1,NXHALF
        ARGMNT = TWOPIN*REAL(KK)
        COSTM2(IX) = COS(ARGMNT)
        SINTM2(IX) = SIN(ARGMNT)
        KK = KK+IROTN
        IF(KK.GT.NI)KK = KK-NI
      ENDDO


      RETURN
      END

