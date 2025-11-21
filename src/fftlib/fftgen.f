      SUBROUTINE FFTGEN(CARRAY,NX)

C     *************************************************************************
C
C     FFTGEN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-MAY-2001:  CREATED
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT AN FFT IN 1D
C     USING TEMPERTON SELF-SORTING PRIME FACTOR ALGORITHM
C     FOR RADIX POWERS OF 2,3,5,7,11
C     BLUESTEIN CONVOLUTIVE ALGORITHM (CHIRP-Z) FOR ANY OTHER SIZE
C     SEPARATE INITIALISATION USING SUBROUTINE FFTGIN
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 13,3,676-686, 1992
C
C     *************************************************************************



C     GLOBAL DATA
C     ===========
C     FGGCOM-------------------------------------------------------------------

      INTEGER NFACTR
      PARAMETER(NFACTR=6)

      INTEGER NIG(NFACTR)

      COMMON/FGGCOM/NIG

C     FGGCOM-------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NX
      DOUBLE PRECISION CARRAY(2*NX)


C     BEGIN
C     =====

C     CALL TEMPERTON RADIX-n FFTS, DEFAULTING TO BLUESTEIN (CHIRP-Z)
      IF(NIG(1).GT.1)CALL FFTFRA(CARRAY,NX)
      IF(NIG(2).GT.1)CALL FFTFR7(CARRAY,NX)
      IF(NIG(3).GT.1)CALL FFTFR5(CARRAY,NX)
      IF(NIG(4).GT.1)CALL FFTFR3(CARRAY,NX)
      IF(NIG(5).GT.1)CALL FFTFR2(CARRAY,NX)
      IF(NIG(6).GT.1)CALL FFTP1D(CARRAY,NX)
      

      RETURN
      END
