      SUBROUTINE FFTGIN(NX,IFORW)

C     *************************************************************************
C
C     FFTGIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-MAY-2001:  CREATED
C     31-DEC-2004:  RSC BUG FIX AND MINOR TIDYING-UP
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT AN FFT IN 1D
C     USING TEMPERTON SELF-SORTING PRIME FACTOR ALGORITHM
C     RADIX POWERS OF 2,3,5,7,11
C     BLUESTEIN CONVOLUTIVE ALGORITHM (CHIRP-Z) FOR ANY OTHER SIZE
C     INITIALISATION ONLY
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
      INTEGER NX,IFORW


C     LOCAL DATA
C     ==========
      INTEGER NSEQ(NFACTR)
      INTEGER NUMBER,IREM,IFACTR,NFACM1


C     BEGIN
C     =====

C     SET FACTORS
      NSEQ(1)=11
      NSEQ(2)=7
      NSEQ(3)=5
      NSEQ(4)=3
      NSEQ(5)=2
 
C     INITIALISE TWIDDLE FACTORS
      DO IFACTR = 1,NFACTR
        NIG(IFACTR) = 1
      ENDDO
      NFACM1 = NFACTR-1

C     FACTORISE
      NUMBER = NX
      DO IFACTR = 1,NFACM1
1000    CONTINUE      
        IREM = MOD(NUMBER,NSEQ(IFACTR))
        IF(IREM.EQ.0)THEN
          NUMBER = NUMBER/NSEQ(IFACTR)
          NIG(IFACTR) = NSEQ(IFACTR)*NIG(IFACTR)
          GOTO 1000
        ENDIF
C       END OF LOOP 1000
      ENDDO

C     CHECK FACTORISATION
C     RSC 31-DEC-2004 BUG FIX DEFAULT TO BLUESTEIN
      IF(NUMBER.GT.1)THEN
        DO IFACTR = 1,NFACM1
          NIG(IFACTR) = 1
        ENDDO
        NIG(NFACTR) = 2
      ENDIF

C     INITIALISE RADIX-n FFTS
      IF(NIG(1).GT.1)CALL FFTNRA(NX,NIG(1),IFORW,NSEQ(1))
      IF(NIG(2).GT.1)CALL FFTNR7(NX,NIG(2),IFORW)
      IF(NIG(3).GT.1)CALL FFTNR5(NX,NIG(3),IFORW)
      IF(NIG(4).GT.1)CALL FFTNR3(NX,NIG(4),IFORW)
      IF(NIG(5).GT.1)CALL FFTNR2(NX,NIG(5),IFORW)
      IF(NIG(6).GT.1)CALL FFTPIN(NX,IFORW)


      RETURN
      END
