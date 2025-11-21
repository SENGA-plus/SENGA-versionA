      SUBROUTINE FFTFR7(CARRAY,NX)

C     *************************************************************************
C
C     FFTFR7
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
C     SEPARATE INITIALISATION
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 13,3,676-686, 1992
C
C     *************************************************************************


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
      INTEGER NX
      DOUBLE PRECISION CARRAY(2*NX)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TWDLRE(6),TWDLIM(6)
      DOUBLE PRECISION TEMPRE(15),TEMPIM(15)
      DOUBLE PRECISION DIFFRE,DIFFIM
      INTEGER IBASE(7,7),JBASE(7,7),IRE(7,7),IIM(7,7)
      INTEGER LA,NXONLA,LABY7
      INTEGER K1,K2,K3,K4,K5,K6
      INTEGER IX,JX,KX,LX
      INTEGER IFLY,INPT
      INTEGER II,LAINC,LAINC7


C     BEGIN
C     =====


C     FIRST HALF
C     ----------
      LA = 1

C     LOOP OVER STAGES
      DO LX = 1,MPHLF7
        NXONLA = NX/LA
        JBASE(1,1) = 0
        JBASE(1,2) = NXSVTH/LA
        JBASE(1,3) = 2*NXSVTH/LA
        JBASE(1,4) = 3*NXSVTH/LA
        JBASE(1,5) = 4*NXSVTH/LA
        JBASE(1,6) = 5*NXSVTH/LA
        JBASE(1,7) = 6*NXSVTH/LA
        K1 = 1
        K2 = 1
        K3 = 1
        K4 = 1
        K5 = 1
        K6 = 1

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRM7

          TWDLRE(1) = COSTM7(K1)
          TWDLIM(1) = SINTM7(K1)
          TWDLRE(2) = COSTM7(K2)
          TWDLIM(2) = SINTM7(K2)
          TWDLRE(3) = COSTM7(K3)
          TWDLIM(3) = SINTM7(K3)
          TWDLRE(4) = COSTM7(K4)
          TWDLIM(4) = SINTM7(K4)
          TWDLRE(5) = COSTM7(K5)
          TWDLIM(5) = SINTM7(K5)
          TWDLRE(6) = COSTM7(K6)
          TWDLIM(6) = SINTM7(K6)

C         LOOP OVER BLOCKS
          DO IX = KX+1,NX,NXONLA

            IBASE(1,1) = JBASE(1,1)+IX
            IBASE(1,2) = JBASE(1,2)+IX
            IBASE(1,3) = JBASE(1,3)+IX
            IBASE(1,4) = JBASE(1,4)+IX
            IBASE(1,5) = JBASE(1,5)+IX
            IBASE(1,6) = JBASE(1,6)+IX
            IBASE(1,7) = JBASE(1,7)+IX

C           LOOP OVER INCREMENTS
            DO II = 1, INCRM7

C             LOCAL INDEXING
              DO INPT = 1,7
                IIM(1,INPT) = 2*IBASE(1,INPT)
                IRE(1,INPT) = IIM(1,INPT)-1
              ENDDO

C             SINGLE BUTTERFLY
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(1)*CARRAY(IRE(1,2))-COEFI7(1)*CARRAY(IIM(1,2))
     +          + COEFR7(2)*CARRAY(IRE(1,3))-COEFI7(2)*CARRAY(IIM(1,3))
     +          + COEFR7(3)*CARRAY(IRE(1,4))-COEFI7(3)*CARRAY(IIM(1,4))
     +          + COEFR7(4)*CARRAY(IRE(1,5))-COEFI7(4)*CARRAY(IIM(1,5))
     +          + COEFR7(5)*CARRAY(IRE(1,6))-COEFI7(5)*CARRAY(IIM(1,6))
     +          + COEFR7(6)*CARRAY(IRE(1,7))-COEFI7(6)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(1)*CARRAY(IRE(1,2))+COEFR7(1)*CARRAY(IIM(1,2))
     +          + COEFI7(2)*CARRAY(IRE(1,3))+COEFR7(2)*CARRAY(IIM(1,3))
     +          + COEFI7(3)*CARRAY(IRE(1,4))+COEFR7(3)*CARRAY(IIM(1,4))
     +          + COEFI7(4)*CARRAY(IRE(1,5))+COEFR7(4)*CARRAY(IIM(1,5))
     +          + COEFI7(5)*CARRAY(IRE(1,6))+COEFR7(5)*CARRAY(IIM(1,6))
     +          + COEFI7(6)*CARRAY(IRE(1,7))+COEFR7(6)*CARRAY(IIM(1,7))
              TEMPRE(1) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              TEMPIM(1) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(2)*CARRAY(IRE(1,2))-COEFI7(2)*CARRAY(IIM(1,2))
     +          + COEFR7(4)*CARRAY(IRE(1,3))-COEFI7(4)*CARRAY(IIM(1,3))
     +          + COEFR7(6)*CARRAY(IRE(1,4))-COEFI7(6)*CARRAY(IIM(1,4))
     +          + COEFR7(1)*CARRAY(IRE(1,5))-COEFI7(1)*CARRAY(IIM(1,5))
     +          + COEFR7(3)*CARRAY(IRE(1,6))-COEFI7(3)*CARRAY(IIM(1,6))
     +          + COEFR7(5)*CARRAY(IRE(1,7))-COEFI7(5)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(2)*CARRAY(IRE(1,2))+COEFR7(2)*CARRAY(IIM(1,2))
     +          + COEFI7(4)*CARRAY(IRE(1,3))+COEFR7(4)*CARRAY(IIM(1,3))
     +          + COEFI7(6)*CARRAY(IRE(1,4))+COEFR7(6)*CARRAY(IIM(1,4))
     +          + COEFI7(1)*CARRAY(IRE(1,5))+COEFR7(1)*CARRAY(IIM(1,5))
     +          + COEFI7(3)*CARRAY(IRE(1,6))+COEFR7(3)*CARRAY(IIM(1,6))
     +          + COEFI7(5)*CARRAY(IRE(1,7))+COEFR7(5)*CARRAY(IIM(1,7))
              TEMPRE(2) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              TEMPIM(2) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(3)*CARRAY(IRE(1,2))-COEFI7(3)*CARRAY(IIM(1,2))
     +          + COEFR7(6)*CARRAY(IRE(1,3))-COEFI7(6)*CARRAY(IIM(1,3))
     +          + COEFR7(2)*CARRAY(IRE(1,4))-COEFI7(2)*CARRAY(IIM(1,4))
     +          + COEFR7(5)*CARRAY(IRE(1,5))-COEFI7(5)*CARRAY(IIM(1,5))
     +          + COEFR7(1)*CARRAY(IRE(1,6))-COEFI7(1)*CARRAY(IIM(1,6))
     +          + COEFR7(4)*CARRAY(IRE(1,7))-COEFI7(4)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(3)*CARRAY(IRE(1,2))+COEFR7(3)*CARRAY(IIM(1,2))
     +          + COEFI7(6)*CARRAY(IRE(1,3))+COEFR7(6)*CARRAY(IIM(1,3))
     +          + COEFI7(2)*CARRAY(IRE(1,4))+COEFR7(2)*CARRAY(IIM(1,4))
     +          + COEFI7(5)*CARRAY(IRE(1,5))+COEFR7(5)*CARRAY(IIM(1,5))
     +          + COEFI7(1)*CARRAY(IRE(1,6))+COEFR7(1)*CARRAY(IIM(1,6))
     +          + COEFI7(4)*CARRAY(IRE(1,7))+COEFR7(4)*CARRAY(IIM(1,7))
              TEMPRE(3) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(3) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(4)*CARRAY(IRE(1,2))-COEFI7(4)*CARRAY(IIM(1,2))
     +          + COEFR7(1)*CARRAY(IRE(1,3))-COEFI7(1)*CARRAY(IIM(1,3))
     +          + COEFR7(5)*CARRAY(IRE(1,4))-COEFI7(5)*CARRAY(IIM(1,4))
     +          + COEFR7(2)*CARRAY(IRE(1,5))-COEFI7(2)*CARRAY(IIM(1,5))
     +          + COEFR7(6)*CARRAY(IRE(1,6))-COEFI7(6)*CARRAY(IIM(1,6))
     +          + COEFR7(3)*CARRAY(IRE(1,7))-COEFI7(3)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(4)*CARRAY(IRE(1,2))+COEFR7(4)*CARRAY(IIM(1,2))
     +          + COEFI7(1)*CARRAY(IRE(1,3))+COEFR7(1)*CARRAY(IIM(1,3))
     +          + COEFI7(5)*CARRAY(IRE(1,4))+COEFR7(5)*CARRAY(IIM(1,4))
     +          + COEFI7(2)*CARRAY(IRE(1,5))+COEFR7(2)*CARRAY(IIM(1,5))
     +          + COEFI7(6)*CARRAY(IRE(1,6))+COEFR7(6)*CARRAY(IIM(1,6))
     +          + COEFI7(3)*CARRAY(IRE(1,7))+COEFR7(3)*CARRAY(IIM(1,7))
              TEMPRE(4) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(4) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(5)*CARRAY(IRE(1,2))-COEFI7(5)*CARRAY(IIM(1,2))
     +          + COEFR7(3)*CARRAY(IRE(1,3))-COEFI7(3)*CARRAY(IIM(1,3))
     +          + COEFR7(1)*CARRAY(IRE(1,4))-COEFI7(1)*CARRAY(IIM(1,4))
     +          + COEFR7(6)*CARRAY(IRE(1,5))-COEFI7(6)*CARRAY(IIM(1,5))
     +          + COEFR7(4)*CARRAY(IRE(1,6))-COEFI7(4)*CARRAY(IIM(1,6))
     +          + COEFR7(2)*CARRAY(IRE(1,7))-COEFI7(2)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(5)*CARRAY(IRE(1,2))+COEFR7(5)*CARRAY(IIM(1,2))
     +          + COEFI7(3)*CARRAY(IRE(1,3))+COEFR7(3)*CARRAY(IIM(1,3))
     +          + COEFI7(1)*CARRAY(IRE(1,4))+COEFR7(1)*CARRAY(IIM(1,4))
     +          + COEFI7(6)*CARRAY(IRE(1,5))+COEFR7(6)*CARRAY(IIM(1,5))
     +          + COEFI7(4)*CARRAY(IRE(1,6))+COEFR7(4)*CARRAY(IIM(1,6))
     +          + COEFI7(2)*CARRAY(IRE(1,7))+COEFR7(2)*CARRAY(IIM(1,7))
              TEMPRE(5) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              TEMPIM(5) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(6)*CARRAY(IRE(1,2))-COEFI7(6)*CARRAY(IIM(1,2))
     +          + COEFR7(5)*CARRAY(IRE(1,3))-COEFI7(5)*CARRAY(IIM(1,3))
     +          + COEFR7(4)*CARRAY(IRE(1,4))-COEFI7(4)*CARRAY(IIM(1,4))
     +          + COEFR7(3)*CARRAY(IRE(1,5))-COEFI7(3)*CARRAY(IIM(1,5))
     +          + COEFR7(2)*CARRAY(IRE(1,6))-COEFI7(2)*CARRAY(IIM(1,6))
     +          + COEFR7(1)*CARRAY(IRE(1,7))-COEFI7(1)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(6)*CARRAY(IRE(1,2))+COEFR7(6)*CARRAY(IIM(1,2))
     +          + COEFI7(5)*CARRAY(IRE(1,3))+COEFR7(5)*CARRAY(IIM(1,3))
     +          + COEFI7(4)*CARRAY(IRE(1,4))+COEFR7(4)*CARRAY(IIM(1,4))
     +          + COEFI7(3)*CARRAY(IRE(1,5))+COEFR7(3)*CARRAY(IIM(1,5))
     +          + COEFI7(2)*CARRAY(IRE(1,6))+COEFR7(2)*CARRAY(IIM(1,6))
     +          + COEFI7(1)*CARRAY(IRE(1,7))+COEFR7(1)*CARRAY(IIM(1,7))
              TEMPRE(6) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              TEMPIM(6) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              CARRAY(IRE(1,1)) = CARRAY(IRE(1,1))
     +                         + CARRAY(IRE(1,2))
     +                         + CARRAY(IRE(1,3))
     +                         + CARRAY(IRE(1,4))
     +                         + CARRAY(IRE(1,5))
     +                         + CARRAY(IRE(1,6))
     +                         + CARRAY(IRE(1,7))
              CARRAY(IIM(1,1)) = CARRAY(IIM(1,1))
     +                         + CARRAY(IIM(1,2))
     +                         + CARRAY(IIM(1,3))
     +                         + CARRAY(IIM(1,4))
     +                         + CARRAY(IIM(1,5))
     +                         + CARRAY(IIM(1,6))
     +                         + CARRAY(IIM(1,7))
              CARRAY(IRE(1,2)) = TEMPRE(1)
              CARRAY(IIM(1,2)) = TEMPIM(1)
              CARRAY(IRE(1,3)) = TEMPRE(2)
              CARRAY(IIM(1,3)) = TEMPIM(2)
              CARRAY(IRE(1,4)) = TEMPRE(3)
              CARRAY(IIM(1,4)) = TEMPIM(3)
              CARRAY(IRE(1,5)) = TEMPRE(4)
              CARRAY(IIM(1,5)) = TEMPIM(4)
              CARRAY(IRE(1,6)) = TEMPRE(5)
              CARRAY(IIM(1,6)) = TEMPIM(5)
              CARRAY(IRE(1,7)) = TEMPRE(6)
              CARRAY(IIM(1,7)) = TEMPIM(6)

              IBASE(1,1) = IBASE(1,1)+NI7
              IBASE(1,2) = IBASE(1,2)+NI7
              IBASE(1,3) = IBASE(1,3)+NI7
              IBASE(1,4) = IBASE(1,4)+NI7
              IBASE(1,5) = IBASE(1,5)+NI7
              IBASE(1,6) = IBASE(1,6)+NI7
              IBASE(1,7) = IBASE(1,7)+NI7
              IF(IBASE(1,1).GT.NX)IBASE(1,1) = IBASE(1,1)-NX
              IF(IBASE(1,2).GT.NX)IBASE(1,2) = IBASE(1,2)-NX
              IF(IBASE(1,3).GT.NX)IBASE(1,3) = IBASE(1,3)-NX
              IF(IBASE(1,4).GT.NX)IBASE(1,4) = IBASE(1,4)-NX
              IF(IBASE(1,5).GT.NX)IBASE(1,5) = IBASE(1,5)-NX
              IF(IBASE(1,6).GT.NX)IBASE(1,6) = IBASE(1,6)-NX
              IF(IBASE(1,7).GT.NX)IBASE(1,7) = IBASE(1,7)-NX

            ENDDO
C           LOOP OVER INCREMENTS

          ENDDO
C         LOOP OVER BLOCKS

          K1 = K1+LA
          K2 = 2*K1-1
          K3 = 3*K1-2
          K4 = 4*K1-3
          K5 = 5*K1-4
          K6 = 6*K1-5

        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = 7*LA

      ENDDO
C     LOOP OVER STAGES


C     SECOND HALF
C     -----------
C     LOOP OVER STAGES
      DO LX = MPHLF7+1,MPOWR7

        NXONLA = NX/LA
        LABY7 = LA*7
        LAINC = LA*INCRM7
        LAINC7 = LABY7*INCRM7
        DO INPT = 1,7
          JBASE(1,INPT) = (INPT-1)*NXSVTH/LA
          DO IFLY = 2,7
            JBASE(IFLY,INPT) = JBASE(IFLY-1,INPT)+LAINC
          ENDDO
        ENDDO
        K1 = 1
        K2 = 1
        K3 = 1
        K4 = 1
        K5 = 1
        K6 = 1

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRM7

          TWDLRE(1) = COSTM7(K1)
          TWDLIM(1) = SINTM7(K1)
          TWDLRE(2) = COSTM7(K2)
          TWDLIM(2) = SINTM7(K2)
          TWDLRE(3) = COSTM7(K3)
          TWDLIM(3) = SINTM7(K3)
          TWDLRE(4) = COSTM7(K4)
          TWDLIM(4) = SINTM7(K4)
          TWDLRE(5) = COSTM7(K5)
          TWDLIM(5) = SINTM7(K5)
          TWDLRE(6) = COSTM7(K6)
          TWDLIM(6) = SINTM7(K6)

C         LOOP OVER SETS OF COUPLED BUTTERFLIES
          DO JX = KX,LAINC-1,NXONLA

C           LOOP OVER UNCOUPLED BLOCKS
            DO IX = JX+1,NX,LAINC7

C             LOCAL INDEXING
              DO INPT = 1,7
                DO IFLY = 1,7
                  IBASE(IFLY,INPT) = JBASE(IFLY,INPT)+IX
                ENDDO
              ENDDO

C             LOOP OVER INCREMENTS
              DO II = 1, INCRM7

C               LOCAL INDEXING
                DO INPT = 1,7
                  DO IFLY = 1,7
                    IIM(IFLY,INPT) = 2*IBASE(IFLY,INPT)
                    IRE(IFLY,INPT) = IIM(IFLY,INPT)-1
                  ENDDO
                ENDDO

C             FIRST BUTTERFLY
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(1)*CARRAY(IRE(1,2))-COEFI7(1)*CARRAY(IIM(1,2))
     +          + COEFR7(2)*CARRAY(IRE(1,3))-COEFI7(2)*CARRAY(IIM(1,3))
     +          + COEFR7(3)*CARRAY(IRE(1,4))-COEFI7(3)*CARRAY(IIM(1,4))
     +          + COEFR7(4)*CARRAY(IRE(1,5))-COEFI7(4)*CARRAY(IIM(1,5))
     +          + COEFR7(5)*CARRAY(IRE(1,6))-COEFI7(5)*CARRAY(IIM(1,6))
     +          + COEFR7(6)*CARRAY(IRE(1,7))-COEFI7(6)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(1)*CARRAY(IRE(1,2))+COEFR7(1)*CARRAY(IIM(1,2))
     +          + COEFI7(2)*CARRAY(IRE(1,3))+COEFR7(2)*CARRAY(IIM(1,3))
     +          + COEFI7(3)*CARRAY(IRE(1,4))+COEFR7(3)*CARRAY(IIM(1,4))
     +          + COEFI7(4)*CARRAY(IRE(1,5))+COEFR7(4)*CARRAY(IIM(1,5))
     +          + COEFI7(5)*CARRAY(IRE(1,6))+COEFR7(5)*CARRAY(IIM(1,6))
     +          + COEFI7(6)*CARRAY(IRE(1,7))+COEFR7(6)*CARRAY(IIM(1,7))
              TEMPRE(1) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              TEMPIM(1) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(2)*CARRAY(IRE(1,2))-COEFI7(2)*CARRAY(IIM(1,2))
     +          + COEFR7(4)*CARRAY(IRE(1,3))-COEFI7(4)*CARRAY(IIM(1,3))
     +          + COEFR7(6)*CARRAY(IRE(1,4))-COEFI7(6)*CARRAY(IIM(1,4))
     +          + COEFR7(1)*CARRAY(IRE(1,5))-COEFI7(1)*CARRAY(IIM(1,5))
     +          + COEFR7(3)*CARRAY(IRE(1,6))-COEFI7(3)*CARRAY(IIM(1,6))
     +          + COEFR7(5)*CARRAY(IRE(1,7))-COEFI7(5)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(2)*CARRAY(IRE(1,2))+COEFR7(2)*CARRAY(IIM(1,2))
     +          + COEFI7(4)*CARRAY(IRE(1,3))+COEFR7(4)*CARRAY(IIM(1,3))
     +          + COEFI7(6)*CARRAY(IRE(1,4))+COEFR7(6)*CARRAY(IIM(1,4))
     +          + COEFI7(1)*CARRAY(IRE(1,5))+COEFR7(1)*CARRAY(IIM(1,5))
     +          + COEFI7(3)*CARRAY(IRE(1,6))+COEFR7(3)*CARRAY(IIM(1,6))
     +          + COEFI7(5)*CARRAY(IRE(1,7))+COEFR7(5)*CARRAY(IIM(1,7))
              TEMPRE(2) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              TEMPIM(2) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(3)*CARRAY(IRE(1,2))-COEFI7(3)*CARRAY(IIM(1,2))
     +          + COEFR7(6)*CARRAY(IRE(1,3))-COEFI7(6)*CARRAY(IIM(1,3))
     +          + COEFR7(2)*CARRAY(IRE(1,4))-COEFI7(2)*CARRAY(IIM(1,4))
     +          + COEFR7(5)*CARRAY(IRE(1,5))-COEFI7(5)*CARRAY(IIM(1,5))
     +          + COEFR7(1)*CARRAY(IRE(1,6))-COEFI7(1)*CARRAY(IIM(1,6))
     +          + COEFR7(4)*CARRAY(IRE(1,7))-COEFI7(4)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(3)*CARRAY(IRE(1,2))+COEFR7(3)*CARRAY(IIM(1,2))
     +          + COEFI7(6)*CARRAY(IRE(1,3))+COEFR7(6)*CARRAY(IIM(1,3))
     +          + COEFI7(2)*CARRAY(IRE(1,4))+COEFR7(2)*CARRAY(IIM(1,4))
     +          + COEFI7(5)*CARRAY(IRE(1,5))+COEFR7(5)*CARRAY(IIM(1,5))
     +          + COEFI7(1)*CARRAY(IRE(1,6))+COEFR7(1)*CARRAY(IIM(1,6))
     +          + COEFI7(4)*CARRAY(IRE(1,7))+COEFR7(4)*CARRAY(IIM(1,7))
              TEMPRE(3) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(3) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(4)*CARRAY(IRE(1,2))-COEFI7(4)*CARRAY(IIM(1,2))
     +          + COEFR7(1)*CARRAY(IRE(1,3))-COEFI7(1)*CARRAY(IIM(1,3))
     +          + COEFR7(5)*CARRAY(IRE(1,4))-COEFI7(5)*CARRAY(IIM(1,4))
     +          + COEFR7(2)*CARRAY(IRE(1,5))-COEFI7(2)*CARRAY(IIM(1,5))
     +          + COEFR7(6)*CARRAY(IRE(1,6))-COEFI7(6)*CARRAY(IIM(1,6))
     +          + COEFR7(3)*CARRAY(IRE(1,7))-COEFI7(3)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(4)*CARRAY(IRE(1,2))+COEFR7(4)*CARRAY(IIM(1,2))
     +          + COEFI7(1)*CARRAY(IRE(1,3))+COEFR7(1)*CARRAY(IIM(1,3))
     +          + COEFI7(5)*CARRAY(IRE(1,4))+COEFR7(5)*CARRAY(IIM(1,4))
     +          + COEFI7(2)*CARRAY(IRE(1,5))+COEFR7(2)*CARRAY(IIM(1,5))
     +          + COEFI7(6)*CARRAY(IRE(1,6))+COEFR7(6)*CARRAY(IIM(1,6))
     +          + COEFI7(3)*CARRAY(IRE(1,7))+COEFR7(3)*CARRAY(IIM(1,7))
                TEMPRE(4) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
                TEMPIM(4) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
                DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(5)*CARRAY(IRE(1,2))-COEFI7(5)*CARRAY(IIM(1,2))
     +          + COEFR7(3)*CARRAY(IRE(1,3))-COEFI7(3)*CARRAY(IIM(1,3))
     +          + COEFR7(1)*CARRAY(IRE(1,4))-COEFI7(1)*CARRAY(IIM(1,4))
     +          + COEFR7(6)*CARRAY(IRE(1,5))-COEFI7(6)*CARRAY(IIM(1,5))
     +          + COEFR7(4)*CARRAY(IRE(1,6))-COEFI7(4)*CARRAY(IIM(1,6))
     +          + COEFR7(2)*CARRAY(IRE(1,7))-COEFI7(2)*CARRAY(IIM(1,7))
                DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(5)*CARRAY(IRE(1,2))+COEFR7(5)*CARRAY(IIM(1,2))
     +          + COEFI7(3)*CARRAY(IRE(1,3))+COEFR7(3)*CARRAY(IIM(1,3))
     +          + COEFI7(1)*CARRAY(IRE(1,4))+COEFR7(1)*CARRAY(IIM(1,4))
     +          + COEFI7(6)*CARRAY(IRE(1,5))+COEFR7(6)*CARRAY(IIM(1,5))
     +          + COEFI7(4)*CARRAY(IRE(1,6))+COEFR7(4)*CARRAY(IIM(1,6))
     +          + COEFI7(2)*CARRAY(IRE(1,7))+COEFR7(2)*CARRAY(IIM(1,7))
              TEMPRE(5) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              TEMPIM(5) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR7(6)*CARRAY(IRE(1,2))-COEFI7(6)*CARRAY(IIM(1,2))
     +          + COEFR7(5)*CARRAY(IRE(1,3))-COEFI7(5)*CARRAY(IIM(1,3))
     +          + COEFR7(4)*CARRAY(IRE(1,4))-COEFI7(4)*CARRAY(IIM(1,4))
     +          + COEFR7(3)*CARRAY(IRE(1,5))-COEFI7(3)*CARRAY(IIM(1,5))
     +          + COEFR7(2)*CARRAY(IRE(1,6))-COEFI7(2)*CARRAY(IIM(1,6))
     +          + COEFR7(1)*CARRAY(IRE(1,7))-COEFI7(1)*CARRAY(IIM(1,7))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI7(6)*CARRAY(IRE(1,2))+COEFR7(6)*CARRAY(IIM(1,2))
     +          + COEFI7(5)*CARRAY(IRE(1,3))+COEFR7(5)*CARRAY(IIM(1,3))
     +          + COEFI7(4)*CARRAY(IRE(1,4))+COEFR7(4)*CARRAY(IIM(1,4))
     +          + COEFI7(3)*CARRAY(IRE(1,5))+COEFR7(3)*CARRAY(IIM(1,5))
     +          + COEFI7(2)*CARRAY(IRE(1,6))+COEFR7(2)*CARRAY(IIM(1,6))
     +          + COEFI7(1)*CARRAY(IRE(1,7))+COEFR7(1)*CARRAY(IIM(1,7))
              TEMPRE(6) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              TEMPIM(6) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              CARRAY(IRE(1,1)) = CARRAY(IRE(1,1))
     +                         + CARRAY(IRE(1,2))
     +                         + CARRAY(IRE(1,3))
     +                         + CARRAY(IRE(1,4))
     +                         + CARRAY(IRE(1,5))
     +                         + CARRAY(IRE(1,6))
     +                         + CARRAY(IRE(1,7))
              CARRAY(IIM(1,1)) = CARRAY(IIM(1,1))
     +                         + CARRAY(IIM(1,2))
     +                         + CARRAY(IIM(1,3))
     +                         + CARRAY(IIM(1,4))
     +                         + CARRAY(IIM(1,5))
     +                         + CARRAY(IIM(1,6))
     +                         + CARRAY(IIM(1,7))

C             SECOND BUTTERFLY
              CARRAY(IRE(1,2)) = CARRAY(IRE(2,1))
     +                         + CARRAY(IRE(2,2))
     +                         + CARRAY(IRE(2,3))
     +                         + CARRAY(IRE(2,4))
     +                         + CARRAY(IRE(2,5))
     +                         + CARRAY(IRE(2,6))
     +                         + CARRAY(IRE(2,7))
              CARRAY(IIM(1,2)) = CARRAY(IIM(2,1))
     +                         + CARRAY(IIM(2,2))
     +                         + CARRAY(IIM(2,3))
     +                         + CARRAY(IIM(2,4))
     +                         + CARRAY(IIM(2,5))
     +                         + CARRAY(IIM(2,6))
     +                         + CARRAY(IIM(2,7))
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR7(2)*CARRAY(IRE(2,2))-COEFI7(2)*CARRAY(IIM(2,2))
     +          + COEFR7(4)*CARRAY(IRE(2,3))-COEFI7(4)*CARRAY(IIM(2,3))
     +          + COEFR7(6)*CARRAY(IRE(2,4))-COEFI7(6)*CARRAY(IIM(2,4))
     +          + COEFR7(1)*CARRAY(IRE(2,5))-COEFI7(1)*CARRAY(IIM(2,5))
     +          + COEFR7(3)*CARRAY(IRE(2,6))-COEFI7(3)*CARRAY(IIM(2,6))
     +          + COEFR7(5)*CARRAY(IRE(2,7))-COEFI7(5)*CARRAY(IIM(2,7))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI7(2)*CARRAY(IRE(2,2))+COEFR7(2)*CARRAY(IIM(2,2))
     +          + COEFI7(4)*CARRAY(IRE(2,3))+COEFR7(4)*CARRAY(IIM(2,3))
     +          + COEFI7(6)*CARRAY(IRE(2,4))+COEFR7(6)*CARRAY(IIM(2,4))
     +          + COEFI7(1)*CARRAY(IRE(2,5))+COEFR7(1)*CARRAY(IIM(2,5))
     +          + COEFI7(3)*CARRAY(IRE(2,6))+COEFR7(3)*CARRAY(IIM(2,6))
     +          + COEFI7(5)*CARRAY(IRE(2,7))+COEFR7(5)*CARRAY(IIM(2,7))
              TEMPRE(7) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              TEMPIM(7) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR7(3)*CARRAY(IRE(2,2))-COEFI7(3)*CARRAY(IIM(2,2))
     +          + COEFR7(6)*CARRAY(IRE(2,3))-COEFI7(6)*CARRAY(IIM(2,3))
     +          + COEFR7(2)*CARRAY(IRE(2,4))-COEFI7(2)*CARRAY(IIM(2,4))
     +          + COEFR7(5)*CARRAY(IRE(2,5))-COEFI7(5)*CARRAY(IIM(2,5))
     +          + COEFR7(1)*CARRAY(IRE(2,6))-COEFI7(1)*CARRAY(IIM(2,6))
     +          + COEFR7(4)*CARRAY(IRE(2,7))-COEFI7(4)*CARRAY(IIM(2,7))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI7(3)*CARRAY(IRE(2,2))+COEFR7(3)*CARRAY(IIM(2,2))
     +          + COEFI7(6)*CARRAY(IRE(2,3))+COEFR7(6)*CARRAY(IIM(2,3))
     +          + COEFI7(2)*CARRAY(IRE(2,4))+COEFR7(2)*CARRAY(IIM(2,4))
     +          + COEFI7(5)*CARRAY(IRE(2,5))+COEFR7(5)*CARRAY(IIM(2,5))
     +          + COEFI7(1)*CARRAY(IRE(2,6))+COEFR7(1)*CARRAY(IIM(2,6))
     +          + COEFI7(4)*CARRAY(IRE(2,7))+COEFR7(4)*CARRAY(IIM(2,7))
              TEMPRE(8) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(8) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR7(4)*CARRAY(IRE(2,2))-COEFI7(4)*CARRAY(IIM(2,2))
     +          + COEFR7(1)*CARRAY(IRE(2,3))-COEFI7(1)*CARRAY(IIM(2,3))
     +          + COEFR7(5)*CARRAY(IRE(2,4))-COEFI7(5)*CARRAY(IIM(2,4))
     +          + COEFR7(2)*CARRAY(IRE(2,5))-COEFI7(2)*CARRAY(IIM(2,5))
     +          + COEFR7(6)*CARRAY(IRE(2,6))-COEFI7(6)*CARRAY(IIM(2,6))
     +          + COEFR7(3)*CARRAY(IRE(2,7))-COEFI7(3)*CARRAY(IIM(2,7))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI7(4)*CARRAY(IRE(2,2))+COEFR7(4)*CARRAY(IIM(2,2))
     +          + COEFI7(1)*CARRAY(IRE(2,3))+COEFR7(1)*CARRAY(IIM(2,3))
     +          + COEFI7(5)*CARRAY(IRE(2,4))+COEFR7(5)*CARRAY(IIM(2,4))
     +          + COEFI7(2)*CARRAY(IRE(2,5))+COEFR7(2)*CARRAY(IIM(2,5))
     +          + COEFI7(6)*CARRAY(IRE(2,6))+COEFR7(6)*CARRAY(IIM(2,6))
     +          + COEFI7(3)*CARRAY(IRE(2,7))+COEFR7(3)*CARRAY(IIM(2,7))
              TEMPRE(9) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(9) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR7(5)*CARRAY(IRE(2,2))-COEFI7(5)*CARRAY(IIM(2,2))
     +          + COEFR7(3)*CARRAY(IRE(2,3))-COEFI7(3)*CARRAY(IIM(2,3))
     +          + COEFR7(1)*CARRAY(IRE(2,4))-COEFI7(1)*CARRAY(IIM(2,4))
     +          + COEFR7(6)*CARRAY(IRE(2,5))-COEFI7(6)*CARRAY(IIM(2,5))
     +          + COEFR7(4)*CARRAY(IRE(2,6))-COEFI7(4)*CARRAY(IIM(2,6))
     +          + COEFR7(2)*CARRAY(IRE(2,7))-COEFI7(2)*CARRAY(IIM(2,7))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI7(5)*CARRAY(IRE(2,2))+COEFR7(5)*CARRAY(IIM(2,2))
     +          + COEFI7(3)*CARRAY(IRE(2,3))+COEFR7(3)*CARRAY(IIM(2,3))
     +          + COEFI7(1)*CARRAY(IRE(2,4))+COEFR7(1)*CARRAY(IIM(2,4))
     +          + COEFI7(6)*CARRAY(IRE(2,5))+COEFR7(6)*CARRAY(IIM(2,5))
     +          + COEFI7(4)*CARRAY(IRE(2,6))+COEFR7(4)*CARRAY(IIM(2,6))
     +          + COEFI7(2)*CARRAY(IRE(2,7))+COEFR7(2)*CARRAY(IIM(2,7))
              TEMPRE(10) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              TEMPIM(10) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR7(6)*CARRAY(IRE(2,2))-COEFI7(6)*CARRAY(IIM(2,2))
     +          + COEFR7(5)*CARRAY(IRE(2,3))-COEFI7(5)*CARRAY(IIM(2,3))
     +          + COEFR7(4)*CARRAY(IRE(2,4))-COEFI7(4)*CARRAY(IIM(2,4))
     +          + COEFR7(3)*CARRAY(IRE(2,5))-COEFI7(3)*CARRAY(IIM(2,5))
     +          + COEFR7(2)*CARRAY(IRE(2,6))-COEFI7(2)*CARRAY(IIM(2,6))
     +          + COEFR7(1)*CARRAY(IRE(2,7))-COEFI7(1)*CARRAY(IIM(2,7))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI7(6)*CARRAY(IRE(2,2))+COEFR7(6)*CARRAY(IIM(2,2))
     +          + COEFI7(5)*CARRAY(IRE(2,3))+COEFR7(5)*CARRAY(IIM(2,3))
     +          + COEFI7(4)*CARRAY(IRE(2,4))+COEFR7(4)*CARRAY(IIM(2,4))
     +          + COEFI7(3)*CARRAY(IRE(2,5))+COEFR7(3)*CARRAY(IIM(2,5))
     +          + COEFI7(2)*CARRAY(IRE(2,6))+COEFR7(2)*CARRAY(IIM(2,6))
     +          + COEFI7(1)*CARRAY(IRE(2,7))+COEFR7(1)*CARRAY(IIM(2,7))
              TEMPRE(11) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              TEMPIM(11) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR7(1)*CARRAY(IRE(2,2))-COEFI7(1)*CARRAY(IIM(2,2))
     +          + COEFR7(2)*CARRAY(IRE(2,3))-COEFI7(2)*CARRAY(IIM(2,3))
     +          + COEFR7(3)*CARRAY(IRE(2,4))-COEFI7(3)*CARRAY(IIM(2,4))
     +          + COEFR7(4)*CARRAY(IRE(2,5))-COEFI7(4)*CARRAY(IIM(2,5))
     +          + COEFR7(5)*CARRAY(IRE(2,6))-COEFI7(5)*CARRAY(IIM(2,6))
     +          + COEFR7(6)*CARRAY(IRE(2,7))-COEFI7(6)*CARRAY(IIM(2,7))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI7(1)*CARRAY(IRE(2,2))+COEFR7(1)*CARRAY(IIM(2,2))
     +          + COEFI7(2)*CARRAY(IRE(2,3))+COEFR7(2)*CARRAY(IIM(2,3))
     +          + COEFI7(3)*CARRAY(IRE(2,4))+COEFR7(3)*CARRAY(IIM(2,4))
     +          + COEFI7(4)*CARRAY(IRE(2,5))+COEFR7(4)*CARRAY(IIM(2,5))
     +          + COEFI7(5)*CARRAY(IRE(2,6))+COEFR7(5)*CARRAY(IIM(2,6))
     +          + COEFI7(6)*CARRAY(IRE(2,7))+COEFR7(6)*CARRAY(IIM(2,7))
              CARRAY(IRE(2,2)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,2)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              CARRAY(IRE(2,1)) = TEMPRE(1)
              CARRAY(IIM(2,1)) = TEMPIM(1)

C             THIRD BUTTERFLY
              CARRAY(IRE(1,3)) = CARRAY(IRE(3,1))
     +                         + CARRAY(IRE(3,2))
     +                         + CARRAY(IRE(3,3))
     +                         + CARRAY(IRE(3,4))
     +                         + CARRAY(IRE(3,5))
     +                         + CARRAY(IRE(3,6))
     +                         + CARRAY(IRE(3,7))
              CARRAY(IIM(1,3)) = CARRAY(IIM(3,1))
     +                         + CARRAY(IIM(3,2))
     +                         + CARRAY(IIM(3,3))
     +                         + CARRAY(IIM(3,4))
     +                         + CARRAY(IIM(3,5))
     +                         + CARRAY(IIM(3,6))
     +                         + CARRAY(IIM(3,7))
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR7(1)*CARRAY(IRE(3,2))-COEFI7(1)*CARRAY(IIM(3,2))
     +          + COEFR7(2)*CARRAY(IRE(3,3))-COEFI7(2)*CARRAY(IIM(3,3))
     +          + COEFR7(3)*CARRAY(IRE(3,4))-COEFI7(3)*CARRAY(IIM(3,4))
     +          + COEFR7(4)*CARRAY(IRE(3,5))-COEFI7(4)*CARRAY(IIM(3,5))
     +          + COEFR7(5)*CARRAY(IRE(3,6))-COEFI7(5)*CARRAY(IIM(3,6))
     +          + COEFR7(6)*CARRAY(IRE(3,7))-COEFI7(6)*CARRAY(IIM(3,7))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI7(1)*CARRAY(IRE(3,2))+COEFR7(1)*CARRAY(IIM(3,2))
     +          + COEFI7(2)*CARRAY(IRE(3,3))+COEFR7(2)*CARRAY(IIM(3,3))
     +          + COEFI7(3)*CARRAY(IRE(3,4))+COEFR7(3)*CARRAY(IIM(3,4))
     +          + COEFI7(4)*CARRAY(IRE(3,5))+COEFR7(4)*CARRAY(IIM(3,5))
     +          + COEFI7(5)*CARRAY(IRE(3,6))+COEFR7(5)*CARRAY(IIM(3,6))
     +          + COEFI7(6)*CARRAY(IRE(3,7))+COEFR7(6)*CARRAY(IIM(3,7))
              CARRAY(IRE(2,3)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,3)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR7(3)*CARRAY(IRE(3,2))-COEFI7(3)*CARRAY(IIM(3,2))
     +          + COEFR7(6)*CARRAY(IRE(3,3))-COEFI7(6)*CARRAY(IIM(3,3))
     +          + COEFR7(2)*CARRAY(IRE(3,4))-COEFI7(2)*CARRAY(IIM(3,4))
     +          + COEFR7(5)*CARRAY(IRE(3,5))-COEFI7(5)*CARRAY(IIM(3,5))
     +          + COEFR7(1)*CARRAY(IRE(3,6))-COEFI7(1)*CARRAY(IIM(3,6))
     +          + COEFR7(4)*CARRAY(IRE(3,7))-COEFI7(4)*CARRAY(IIM(3,7))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI7(3)*CARRAY(IRE(3,2))+COEFR7(3)*CARRAY(IIM(3,2))
     +          + COEFI7(6)*CARRAY(IRE(3,3))+COEFR7(6)*CARRAY(IIM(3,3))
     +          + COEFI7(2)*CARRAY(IRE(3,4))+COEFR7(2)*CARRAY(IIM(3,4))
     +          + COEFI7(5)*CARRAY(IRE(3,5))+COEFR7(5)*CARRAY(IIM(3,5))
     +          + COEFI7(1)*CARRAY(IRE(3,6))+COEFR7(1)*CARRAY(IIM(3,6))
     +          + COEFI7(4)*CARRAY(IRE(3,7))+COEFR7(4)*CARRAY(IIM(3,7))
              TEMPRE(1) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(1) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR7(4)*CARRAY(IRE(3,2))-COEFI7(4)*CARRAY(IIM(3,2))
     +          + COEFR7(1)*CARRAY(IRE(3,3))-COEFI7(1)*CARRAY(IIM(3,3))
     +          + COEFR7(5)*CARRAY(IRE(3,4))-COEFI7(5)*CARRAY(IIM(3,4))
     +          + COEFR7(2)*CARRAY(IRE(3,5))-COEFI7(2)*CARRAY(IIM(3,5))
     +          + COEFR7(6)*CARRAY(IRE(3,6))-COEFI7(6)*CARRAY(IIM(3,6))
     +          + COEFR7(3)*CARRAY(IRE(3,7))-COEFI7(3)*CARRAY(IIM(3,7))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI7(4)*CARRAY(IRE(3,2))+COEFR7(4)*CARRAY(IIM(3,2))
     +          + COEFI7(1)*CARRAY(IRE(3,3))+COEFR7(1)*CARRAY(IIM(3,3))
     +          + COEFI7(5)*CARRAY(IRE(3,4))+COEFR7(5)*CARRAY(IIM(3,4))
     +          + COEFI7(2)*CARRAY(IRE(3,5))+COEFR7(2)*CARRAY(IIM(3,5))
     +          + COEFI7(6)*CARRAY(IRE(3,6))+COEFR7(6)*CARRAY(IIM(3,6))
     +          + COEFI7(3)*CARRAY(IRE(3,7))+COEFR7(3)*CARRAY(IIM(3,7))
              TEMPRE(12) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(12) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR7(5)*CARRAY(IRE(3,2))-COEFI7(5)*CARRAY(IIM(3,2))
     +          + COEFR7(3)*CARRAY(IRE(3,3))-COEFI7(3)*CARRAY(IIM(3,3))
     +          + COEFR7(1)*CARRAY(IRE(3,4))-COEFI7(1)*CARRAY(IIM(3,4))
     +          + COEFR7(6)*CARRAY(IRE(3,5))-COEFI7(6)*CARRAY(IIM(3,5))
     +          + COEFR7(4)*CARRAY(IRE(3,6))-COEFI7(4)*CARRAY(IIM(3,6))
     +          + COEFR7(2)*CARRAY(IRE(3,7))-COEFI7(2)*CARRAY(IIM(3,7))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI7(5)*CARRAY(IRE(3,2))+COEFR7(5)*CARRAY(IIM(3,2))
     +          + COEFI7(3)*CARRAY(IRE(3,3))+COEFR7(3)*CARRAY(IIM(3,3))
     +          + COEFI7(1)*CARRAY(IRE(3,4))+COEFR7(1)*CARRAY(IIM(3,4))
     +          + COEFI7(6)*CARRAY(IRE(3,5))+COEFR7(6)*CARRAY(IIM(3,5))
     +          + COEFI7(4)*CARRAY(IRE(3,6))+COEFR7(4)*CARRAY(IIM(3,6))
     +          + COEFI7(2)*CARRAY(IRE(3,7))+COEFR7(2)*CARRAY(IIM(3,7))
              TEMPRE(13) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              TEMPIM(13) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR7(6)*CARRAY(IRE(3,2))-COEFI7(6)*CARRAY(IIM(3,2))
     +          + COEFR7(5)*CARRAY(IRE(3,3))-COEFI7(5)*CARRAY(IIM(3,3))
     +          + COEFR7(4)*CARRAY(IRE(3,4))-COEFI7(4)*CARRAY(IIM(3,4))
     +          + COEFR7(3)*CARRAY(IRE(3,5))-COEFI7(3)*CARRAY(IIM(3,5))
     +          + COEFR7(2)*CARRAY(IRE(3,6))-COEFI7(2)*CARRAY(IIM(3,6))
     +          + COEFR7(1)*CARRAY(IRE(3,7))-COEFI7(1)*CARRAY(IIM(3,7))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI7(6)*CARRAY(IRE(3,2))+COEFR7(6)*CARRAY(IIM(3,2))
     +          + COEFI7(5)*CARRAY(IRE(3,3))+COEFR7(5)*CARRAY(IIM(3,3))
     +          + COEFI7(4)*CARRAY(IRE(3,4))+COEFR7(4)*CARRAY(IIM(3,4))
     +          + COEFI7(3)*CARRAY(IRE(3,5))+COEFR7(3)*CARRAY(IIM(3,5))
     +          + COEFI7(2)*CARRAY(IRE(3,6))+COEFR7(2)*CARRAY(IIM(3,6))
     +          + COEFI7(1)*CARRAY(IRE(3,7))+COEFR7(1)*CARRAY(IIM(3,7))
              TEMPRE(14) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              TEMPIM(14) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR7(2)*CARRAY(IRE(3,2))-COEFI7(2)*CARRAY(IIM(3,2))
     +          + COEFR7(4)*CARRAY(IRE(3,3))-COEFI7(4)*CARRAY(IIM(3,3))
     +          + COEFR7(6)*CARRAY(IRE(3,4))-COEFI7(6)*CARRAY(IIM(3,4))
     +          + COEFR7(1)*CARRAY(IRE(3,5))-COEFI7(1)*CARRAY(IIM(3,5))
     +          + COEFR7(3)*CARRAY(IRE(3,6))-COEFI7(3)*CARRAY(IIM(3,6))
     +          + COEFR7(5)*CARRAY(IRE(3,7))-COEFI7(5)*CARRAY(IIM(3,7))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI7(2)*CARRAY(IRE(3,2))+COEFR7(2)*CARRAY(IIM(3,2))
     +          + COEFI7(4)*CARRAY(IRE(3,3))+COEFR7(4)*CARRAY(IIM(3,3))
     +          + COEFI7(6)*CARRAY(IRE(3,4))+COEFR7(6)*CARRAY(IIM(3,4))
     +          + COEFI7(1)*CARRAY(IRE(3,5))+COEFR7(1)*CARRAY(IIM(3,5))
     +          + COEFI7(3)*CARRAY(IRE(3,6))+COEFR7(3)*CARRAY(IIM(3,6))
     +          + COEFI7(5)*CARRAY(IRE(3,7))+COEFR7(5)*CARRAY(IIM(3,7))
              CARRAY(IRE(3,3)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,3)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              CARRAY(IRE(3,1)) = TEMPRE(2)
              CARRAY(IIM(3,1)) = TEMPIM(2)
              CARRAY(IRE(3,2)) = TEMPRE(7)
              CARRAY(IIM(3,2)) = TEMPIM(7)

C             FOURTH BUTTERFLY
              CARRAY(IRE(1,4)) = CARRAY(IRE(4,1))
     +                         + CARRAY(IRE(4,2))
     +                         + CARRAY(IRE(4,3))
     +                         + CARRAY(IRE(4,4))
     +                         + CARRAY(IRE(4,5))
     +                         + CARRAY(IRE(4,6))
     +                         + CARRAY(IRE(4,7))
              CARRAY(IIM(1,4)) = CARRAY(IIM(4,1))
     +                         + CARRAY(IIM(4,2))
     +                         + CARRAY(IIM(4,3))
     +                         + CARRAY(IIM(4,4))
     +                         + CARRAY(IIM(4,5))
     +                         + CARRAY(IIM(4,6))
     +                         + CARRAY(IIM(4,7))
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR7(1)*CARRAY(IRE(4,2))-COEFI7(1)*CARRAY(IIM(4,2))
     +          + COEFR7(2)*CARRAY(IRE(4,3))-COEFI7(2)*CARRAY(IIM(4,3))
     +          + COEFR7(3)*CARRAY(IRE(4,4))-COEFI7(3)*CARRAY(IIM(4,4))
     +          + COEFR7(4)*CARRAY(IRE(4,5))-COEFI7(4)*CARRAY(IIM(4,5))
     +          + COEFR7(5)*CARRAY(IRE(4,6))-COEFI7(5)*CARRAY(IIM(4,6))
     +          + COEFR7(6)*CARRAY(IRE(4,7))-COEFI7(6)*CARRAY(IIM(4,7))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI7(1)*CARRAY(IRE(4,2))+COEFR7(1)*CARRAY(IIM(4,2))
     +          + COEFI7(2)*CARRAY(IRE(4,3))+COEFR7(2)*CARRAY(IIM(4,3))
     +          + COEFI7(3)*CARRAY(IRE(4,4))+COEFR7(3)*CARRAY(IIM(4,4))
     +          + COEFI7(4)*CARRAY(IRE(4,5))+COEFR7(4)*CARRAY(IIM(4,5))
     +          + COEFI7(5)*CARRAY(IRE(4,6))+COEFR7(5)*CARRAY(IIM(4,6))
     +          + COEFI7(6)*CARRAY(IRE(4,7))+COEFR7(6)*CARRAY(IIM(4,7))
              CARRAY(IRE(2,4)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,4)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR7(2)*CARRAY(IRE(4,2))-COEFI7(2)*CARRAY(IIM(4,2))
     +          + COEFR7(4)*CARRAY(IRE(4,3))-COEFI7(4)*CARRAY(IIM(4,3))
     +          + COEFR7(6)*CARRAY(IRE(4,4))-COEFI7(6)*CARRAY(IIM(4,4))
     +          + COEFR7(1)*CARRAY(IRE(4,5))-COEFI7(1)*CARRAY(IIM(4,5))
     +          + COEFR7(3)*CARRAY(IRE(4,6))-COEFI7(3)*CARRAY(IIM(4,6))
     +          + COEFR7(5)*CARRAY(IRE(4,7))-COEFI7(5)*CARRAY(IIM(4,7))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI7(2)*CARRAY(IRE(4,2))+COEFR7(2)*CARRAY(IIM(4,2))
     +          + COEFI7(4)*CARRAY(IRE(4,3))+COEFR7(4)*CARRAY(IIM(4,3))
     +          + COEFI7(6)*CARRAY(IRE(4,4))+COEFR7(6)*CARRAY(IIM(4,4))
     +          + COEFI7(1)*CARRAY(IRE(4,5))+COEFR7(1)*CARRAY(IIM(4,5))
     +          + COEFI7(3)*CARRAY(IRE(4,6))+COEFR7(3)*CARRAY(IIM(4,6))
     +          + COEFI7(5)*CARRAY(IRE(4,7))+COEFR7(5)*CARRAY(IIM(4,7))
              CARRAY(IRE(3,4)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,4)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR7(4)*CARRAY(IRE(4,2))-COEFI7(4)*CARRAY(IIM(4,2))
     +          + COEFR7(1)*CARRAY(IRE(4,3))-COEFI7(1)*CARRAY(IIM(4,3))
     +          + COEFR7(5)*CARRAY(IRE(4,4))-COEFI7(5)*CARRAY(IIM(4,4))
     +          + COEFR7(2)*CARRAY(IRE(4,5))-COEFI7(2)*CARRAY(IIM(4,5))
     +          + COEFR7(6)*CARRAY(IRE(4,6))-COEFI7(6)*CARRAY(IIM(4,6))
     +          + COEFR7(3)*CARRAY(IRE(4,7))-COEFI7(3)*CARRAY(IIM(4,7))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI7(4)*CARRAY(IRE(4,2))+COEFR7(4)*CARRAY(IIM(4,2))
     +          + COEFI7(1)*CARRAY(IRE(4,3))+COEFR7(1)*CARRAY(IIM(4,3))
     +          + COEFI7(5)*CARRAY(IRE(4,4))+COEFR7(5)*CARRAY(IIM(4,4))
     +          + COEFI7(2)*CARRAY(IRE(4,5))+COEFR7(2)*CARRAY(IIM(4,5))
     +          + COEFI7(6)*CARRAY(IRE(4,6))+COEFR7(6)*CARRAY(IIM(4,6))
     +          + COEFI7(3)*CARRAY(IRE(4,7))+COEFR7(3)*CARRAY(IIM(4,7))
              TEMPRE(2) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(2) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR7(5)*CARRAY(IRE(4,2))-COEFI7(5)*CARRAY(IIM(4,2))
     +          + COEFR7(3)*CARRAY(IRE(4,3))-COEFI7(3)*CARRAY(IIM(4,3))
     +          + COEFR7(1)*CARRAY(IRE(4,4))-COEFI7(1)*CARRAY(IIM(4,4))
     +          + COEFR7(6)*CARRAY(IRE(4,5))-COEFI7(6)*CARRAY(IIM(4,5))
     +          + COEFR7(4)*CARRAY(IRE(4,6))-COEFI7(4)*CARRAY(IIM(4,6))
     +          + COEFR7(2)*CARRAY(IRE(4,7))-COEFI7(2)*CARRAY(IIM(4,7))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI7(5)*CARRAY(IRE(4,2))+COEFR7(5)*CARRAY(IIM(4,2))
     +          + COEFI7(3)*CARRAY(IRE(4,3))+COEFR7(3)*CARRAY(IIM(4,3))
     +          + COEFI7(1)*CARRAY(IRE(4,4))+COEFR7(1)*CARRAY(IIM(4,4))
     +          + COEFI7(6)*CARRAY(IRE(4,5))+COEFR7(6)*CARRAY(IIM(4,5))
     +          + COEFI7(4)*CARRAY(IRE(4,6))+COEFR7(4)*CARRAY(IIM(4,6))
     +          + COEFI7(2)*CARRAY(IRE(4,7))+COEFR7(2)*CARRAY(IIM(4,7))
              TEMPRE(7) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              TEMPIM(7) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR7(6)*CARRAY(IRE(4,2))-COEFI7(6)*CARRAY(IIM(4,2))
     +          + COEFR7(5)*CARRAY(IRE(4,3))-COEFI7(5)*CARRAY(IIM(4,3))
     +          + COEFR7(4)*CARRAY(IRE(4,4))-COEFI7(4)*CARRAY(IIM(4,4))
     +          + COEFR7(3)*CARRAY(IRE(4,5))-COEFI7(3)*CARRAY(IIM(4,5))
     +          + COEFR7(2)*CARRAY(IRE(4,6))-COEFI7(2)*CARRAY(IIM(4,6))
     +          + COEFR7(1)*CARRAY(IRE(4,7))-COEFI7(1)*CARRAY(IIM(4,7))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI7(6)*CARRAY(IRE(4,2))+COEFR7(6)*CARRAY(IIM(4,2))
     +          + COEFI7(5)*CARRAY(IRE(4,3))+COEFR7(5)*CARRAY(IIM(4,3))
     +          + COEFI7(4)*CARRAY(IRE(4,4))+COEFR7(4)*CARRAY(IIM(4,4))
     +          + COEFI7(3)*CARRAY(IRE(4,5))+COEFR7(3)*CARRAY(IIM(4,5))
     +          + COEFI7(2)*CARRAY(IRE(4,6))+COEFR7(2)*CARRAY(IIM(4,6))
     +          + COEFI7(1)*CARRAY(IRE(4,7))+COEFR7(1)*CARRAY(IIM(4,7))
              TEMPRE(15) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              TEMPIM(15) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR7(3)*CARRAY(IRE(4,2))-COEFI7(3)*CARRAY(IIM(4,2))
     +          + COEFR7(6)*CARRAY(IRE(4,3))-COEFI7(6)*CARRAY(IIM(4,3))
     +          + COEFR7(2)*CARRAY(IRE(4,4))-COEFI7(2)*CARRAY(IIM(4,4))
     +          + COEFR7(5)*CARRAY(IRE(4,5))-COEFI7(5)*CARRAY(IIM(4,5))
     +          + COEFR7(1)*CARRAY(IRE(4,6))-COEFI7(1)*CARRAY(IIM(4,6))
     +          + COEFR7(4)*CARRAY(IRE(4,7))-COEFI7(4)*CARRAY(IIM(4,7))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI7(3)*CARRAY(IRE(4,2))+COEFR7(3)*CARRAY(IIM(4,2))
     +          + COEFI7(6)*CARRAY(IRE(4,3))+COEFR7(6)*CARRAY(IIM(4,3))
     +          + COEFI7(2)*CARRAY(IRE(4,4))+COEFR7(2)*CARRAY(IIM(4,4))
     +          + COEFI7(5)*CARRAY(IRE(4,5))+COEFR7(5)*CARRAY(IIM(4,5))
     +          + COEFI7(1)*CARRAY(IRE(4,6))+COEFR7(1)*CARRAY(IIM(4,6))
     +          + COEFI7(4)*CARRAY(IRE(4,7))+COEFR7(4)*CARRAY(IIM(4,7))
              CARRAY(IRE(4,4)) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              CARRAY(IIM(4,4)) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              CARRAY(IRE(4,1)) = TEMPRE(3)
              CARRAY(IIM(4,1)) = TEMPIM(3)
              CARRAY(IRE(4,2)) = TEMPRE(8)
              CARRAY(IIM(4,2)) = TEMPIM(8)
              CARRAY(IRE(4,3)) = TEMPRE(1)
              CARRAY(IIM(4,3)) = TEMPIM(1)

C             FIFTH BUTTERFLY
              CARRAY(IRE(1,5)) = CARRAY(IRE(5,1))
     +                         + CARRAY(IRE(5,2))
     +                         + CARRAY(IRE(5,3))
     +                         + CARRAY(IRE(5,4))
     +                         + CARRAY(IRE(5,5))
     +                         + CARRAY(IRE(5,6))
     +                         + CARRAY(IRE(5,7))
              CARRAY(IIM(1,5)) = CARRAY(IIM(5,1))
     +                         + CARRAY(IIM(5,2))
     +                         + CARRAY(IIM(5,3))
     +                         + CARRAY(IIM(5,4))
     +                         + CARRAY(IIM(5,5))
     +                         + CARRAY(IIM(5,6))
     +                         + CARRAY(IIM(5,7))
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR7(1)*CARRAY(IRE(5,2))-COEFI7(1)*CARRAY(IIM(5,2))
     +          + COEFR7(2)*CARRAY(IRE(5,3))-COEFI7(2)*CARRAY(IIM(5,3))
     +          + COEFR7(3)*CARRAY(IRE(5,4))-COEFI7(3)*CARRAY(IIM(5,4))
     +          + COEFR7(4)*CARRAY(IRE(5,5))-COEFI7(4)*CARRAY(IIM(5,5))
     +          + COEFR7(5)*CARRAY(IRE(5,6))-COEFI7(5)*CARRAY(IIM(5,6))
     +          + COEFR7(6)*CARRAY(IRE(5,7))-COEFI7(6)*CARRAY(IIM(5,7))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI7(1)*CARRAY(IRE(5,2))+COEFR7(1)*CARRAY(IIM(5,2))
     +          + COEFI7(2)*CARRAY(IRE(5,3))+COEFR7(2)*CARRAY(IIM(5,3))
     +          + COEFI7(3)*CARRAY(IRE(5,4))+COEFR7(3)*CARRAY(IIM(5,4))
     +          + COEFI7(4)*CARRAY(IRE(5,5))+COEFR7(4)*CARRAY(IIM(5,5))
     +          + COEFI7(5)*CARRAY(IRE(5,6))+COEFR7(5)*CARRAY(IIM(5,6))
     +          + COEFI7(6)*CARRAY(IRE(5,7))+COEFR7(6)*CARRAY(IIM(5,7))
              CARRAY(IRE(2,5)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,5)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR7(2)*CARRAY(IRE(5,2))-COEFI7(2)*CARRAY(IIM(5,2))
     +          + COEFR7(4)*CARRAY(IRE(5,3))-COEFI7(4)*CARRAY(IIM(5,3))
     +          + COEFR7(6)*CARRAY(IRE(5,4))-COEFI7(6)*CARRAY(IIM(5,4))
     +          + COEFR7(1)*CARRAY(IRE(5,5))-COEFI7(1)*CARRAY(IIM(5,5))
     +          + COEFR7(3)*CARRAY(IRE(5,6))-COEFI7(3)*CARRAY(IIM(5,6))
     +          + COEFR7(5)*CARRAY(IRE(5,7))-COEFI7(5)*CARRAY(IIM(5,7))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI7(2)*CARRAY(IRE(5,2))+COEFR7(2)*CARRAY(IIM(5,2))
     +          + COEFI7(4)*CARRAY(IRE(5,3))+COEFR7(4)*CARRAY(IIM(5,3))
     +          + COEFI7(6)*CARRAY(IRE(5,4))+COEFR7(6)*CARRAY(IIM(5,4))
     +          + COEFI7(1)*CARRAY(IRE(5,5))+COEFR7(1)*CARRAY(IIM(5,5))
     +          + COEFI7(3)*CARRAY(IRE(5,6))+COEFR7(3)*CARRAY(IIM(5,6))
     +          + COEFI7(5)*CARRAY(IRE(5,7))+COEFR7(5)*CARRAY(IIM(5,7))
              CARRAY(IRE(3,5)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,5)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR7(3)*CARRAY(IRE(5,2))-COEFI7(3)*CARRAY(IIM(5,2))
     +          + COEFR7(6)*CARRAY(IRE(5,3))-COEFI7(6)*CARRAY(IIM(5,3))
     +          + COEFR7(2)*CARRAY(IRE(5,4))-COEFI7(2)*CARRAY(IIM(5,4))
     +          + COEFR7(5)*CARRAY(IRE(5,5))-COEFI7(5)*CARRAY(IIM(5,5))
     +          + COEFR7(1)*CARRAY(IRE(5,6))-COEFI7(1)*CARRAY(IIM(5,6))
     +          + COEFR7(4)*CARRAY(IRE(5,7))-COEFI7(4)*CARRAY(IIM(5,7))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI7(3)*CARRAY(IRE(5,2))+COEFR7(3)*CARRAY(IIM(5,2))
     +          + COEFI7(6)*CARRAY(IRE(5,3))+COEFR7(6)*CARRAY(IIM(5,3))
     +          + COEFI7(2)*CARRAY(IRE(5,4))+COEFR7(2)*CARRAY(IIM(5,4))
     +          + COEFI7(5)*CARRAY(IRE(5,5))+COEFR7(5)*CARRAY(IIM(5,5))
     +          + COEFI7(1)*CARRAY(IRE(5,6))+COEFR7(1)*CARRAY(IIM(5,6))
     +          + COEFI7(4)*CARRAY(IRE(5,7))+COEFR7(4)*CARRAY(IIM(5,7))
              CARRAY(IRE(4,5)) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              CARRAY(IIM(4,5)) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR7(5)*CARRAY(IRE(5,2))-COEFI7(5)*CARRAY(IIM(5,2))
     +          + COEFR7(3)*CARRAY(IRE(5,3))-COEFI7(3)*CARRAY(IIM(5,3))
     +          + COEFR7(1)*CARRAY(IRE(5,4))-COEFI7(1)*CARRAY(IIM(5,4))
     +          + COEFR7(6)*CARRAY(IRE(5,5))-COEFI7(6)*CARRAY(IIM(5,5))
     +          + COEFR7(4)*CARRAY(IRE(5,6))-COEFI7(4)*CARRAY(IIM(5,6))
     +          + COEFR7(2)*CARRAY(IRE(5,7))-COEFI7(2)*CARRAY(IIM(5,7))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI7(5)*CARRAY(IRE(5,2))+COEFR7(5)*CARRAY(IIM(5,2))
     +          + COEFI7(3)*CARRAY(IRE(5,3))+COEFR7(3)*CARRAY(IIM(5,3))
     +          + COEFI7(1)*CARRAY(IRE(5,4))+COEFR7(1)*CARRAY(IIM(5,4))
     +          + COEFI7(6)*CARRAY(IRE(5,5))+COEFR7(6)*CARRAY(IIM(5,5))
     +          + COEFI7(4)*CARRAY(IRE(5,6))+COEFR7(4)*CARRAY(IIM(5,6))
     +          + COEFI7(2)*CARRAY(IRE(5,7))+COEFR7(2)*CARRAY(IIM(5,7))
              TEMPRE(1) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              TEMPIM(1) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR7(6)*CARRAY(IRE(5,2))-COEFI7(6)*CARRAY(IIM(5,2))
     +          + COEFR7(5)*CARRAY(IRE(5,3))-COEFI7(5)*CARRAY(IIM(5,3))
     +          + COEFR7(4)*CARRAY(IRE(5,4))-COEFI7(4)*CARRAY(IIM(5,4))
     +          + COEFR7(3)*CARRAY(IRE(5,5))-COEFI7(3)*CARRAY(IIM(5,5))
     +          + COEFR7(2)*CARRAY(IRE(5,6))-COEFI7(2)*CARRAY(IIM(5,6))
     +          + COEFR7(1)*CARRAY(IRE(5,7))-COEFI7(1)*CARRAY(IIM(5,7))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI7(6)*CARRAY(IRE(5,2))+COEFR7(6)*CARRAY(IIM(5,2))
     +          + COEFI7(5)*CARRAY(IRE(5,3))+COEFR7(5)*CARRAY(IIM(5,3))
     +          + COEFI7(4)*CARRAY(IRE(5,4))+COEFR7(4)*CARRAY(IIM(5,4))
     +          + COEFI7(3)*CARRAY(IRE(5,5))+COEFR7(3)*CARRAY(IIM(5,5))
     +          + COEFI7(2)*CARRAY(IRE(5,6))+COEFR7(2)*CARRAY(IIM(5,6))
     +          + COEFI7(1)*CARRAY(IRE(5,7))+COEFR7(1)*CARRAY(IIM(5,7))
              TEMPRE(3) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              TEMPIM(3) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR7(4)*CARRAY(IRE(5,2))-COEFI7(4)*CARRAY(IIM(5,2))
     +          + COEFR7(1)*CARRAY(IRE(5,3))-COEFI7(1)*CARRAY(IIM(5,3))
     +          + COEFR7(5)*CARRAY(IRE(5,4))-COEFI7(5)*CARRAY(IIM(5,4))
     +          + COEFR7(2)*CARRAY(IRE(5,5))-COEFI7(2)*CARRAY(IIM(5,5))
     +          + COEFR7(6)*CARRAY(IRE(5,6))-COEFI7(6)*CARRAY(IIM(5,6))
     +          + COEFR7(3)*CARRAY(IRE(5,7))-COEFI7(3)*CARRAY(IIM(5,7))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI7(4)*CARRAY(IRE(5,2))+COEFR7(4)*CARRAY(IIM(5,2))
     +          + COEFI7(1)*CARRAY(IRE(5,3))+COEFR7(1)*CARRAY(IIM(5,3))
     +          + COEFI7(5)*CARRAY(IRE(5,4))+COEFR7(5)*CARRAY(IIM(5,4))
     +          + COEFI7(2)*CARRAY(IRE(5,5))+COEFR7(2)*CARRAY(IIM(5,5))
     +          + COEFI7(6)*CARRAY(IRE(5,6))+COEFR7(6)*CARRAY(IIM(5,6))
     +          + COEFI7(3)*CARRAY(IRE(5,7))+COEFR7(3)*CARRAY(IIM(5,7))
              CARRAY(IRE(5,5)) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              CARRAY(IIM(5,5)) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              CARRAY(IRE(5,1)) = TEMPRE(4)
              CARRAY(IIM(5,1)) = TEMPIM(4)
              CARRAY(IRE(5,2)) = TEMPRE(9)
              CARRAY(IIM(5,2)) = TEMPIM(9)
              CARRAY(IRE(5,3)) = TEMPRE(12)
              CARRAY(IIM(5,3)) = TEMPIM(12)
              CARRAY(IRE(5,4)) = TEMPRE(2)
              CARRAY(IIM(5,4)) = TEMPIM(2)

C             SIXTH BUTTERFLY
              CARRAY(IRE(1,6)) = CARRAY(IRE(6,1))
     +                         + CARRAY(IRE(6,2))
     +                         + CARRAY(IRE(6,3))
     +                         + CARRAY(IRE(6,4))
     +                         + CARRAY(IRE(6,5))
     +                         + CARRAY(IRE(6,6))
     +                         + CARRAY(IRE(6,7))
              CARRAY(IIM(1,6)) = CARRAY(IIM(6,1))
     +                         + CARRAY(IIM(6,2))
     +                         + CARRAY(IIM(6,3))
     +                         + CARRAY(IIM(6,4))
     +                         + CARRAY(IIM(6,5))
     +                         + CARRAY(IIM(6,6))
     +                         + CARRAY(IIM(6,7))
              DIFFRE = CARRAY(IRE(6,1))
     +          + COEFR7(1)*CARRAY(IRE(6,2))-COEFI7(1)*CARRAY(IIM(6,2))
     +          + COEFR7(2)*CARRAY(IRE(6,3))-COEFI7(2)*CARRAY(IIM(6,3))
     +          + COEFR7(3)*CARRAY(IRE(6,4))-COEFI7(3)*CARRAY(IIM(6,4))
     +          + COEFR7(4)*CARRAY(IRE(6,5))-COEFI7(4)*CARRAY(IIM(6,5))
     +          + COEFR7(5)*CARRAY(IRE(6,6))-COEFI7(5)*CARRAY(IIM(6,6))
     +          + COEFR7(6)*CARRAY(IRE(6,7))-COEFI7(6)*CARRAY(IIM(6,7))
              DIFFIM = CARRAY(IIM(6,1))
     +          + COEFI7(1)*CARRAY(IRE(6,2))+COEFR7(1)*CARRAY(IIM(6,2))
     +          + COEFI7(2)*CARRAY(IRE(6,3))+COEFR7(2)*CARRAY(IIM(6,3))
     +          + COEFI7(3)*CARRAY(IRE(6,4))+COEFR7(3)*CARRAY(IIM(6,4))
     +          + COEFI7(4)*CARRAY(IRE(6,5))+COEFR7(4)*CARRAY(IIM(6,5))
     +          + COEFI7(5)*CARRAY(IRE(6,6))+COEFR7(5)*CARRAY(IIM(6,6))
     +          + COEFI7(6)*CARRAY(IRE(6,7))+COEFR7(6)*CARRAY(IIM(6,7))
              CARRAY(IRE(2,6)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,6)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(6,1))
     +          + COEFR7(2)*CARRAY(IRE(6,2))-COEFI7(2)*CARRAY(IIM(6,2))
     +          + COEFR7(4)*CARRAY(IRE(6,3))-COEFI7(4)*CARRAY(IIM(6,3))
     +          + COEFR7(6)*CARRAY(IRE(6,4))-COEFI7(6)*CARRAY(IIM(6,4))
     +          + COEFR7(1)*CARRAY(IRE(6,5))-COEFI7(1)*CARRAY(IIM(6,5))
     +          + COEFR7(3)*CARRAY(IRE(6,6))-COEFI7(3)*CARRAY(IIM(6,6))
     +          + COEFR7(5)*CARRAY(IRE(6,7))-COEFI7(5)*CARRAY(IIM(6,7))
              DIFFIM = CARRAY(IIM(6,1))
     +          + COEFI7(2)*CARRAY(IRE(6,2))+COEFR7(2)*CARRAY(IIM(6,2))
     +          + COEFI7(4)*CARRAY(IRE(6,3))+COEFR7(4)*CARRAY(IIM(6,3))
     +          + COEFI7(6)*CARRAY(IRE(6,4))+COEFR7(6)*CARRAY(IIM(6,4))
     +          + COEFI7(1)*CARRAY(IRE(6,5))+COEFR7(1)*CARRAY(IIM(6,5))
     +          + COEFI7(3)*CARRAY(IRE(6,6))+COEFR7(3)*CARRAY(IIM(6,6))
     +          + COEFI7(5)*CARRAY(IRE(6,7))+COEFR7(5)*CARRAY(IIM(6,7))
              CARRAY(IRE(3,6)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,6)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(6,1))
     +          + COEFR7(3)*CARRAY(IRE(6,2))-COEFI7(3)*CARRAY(IIM(6,2))
     +          + COEFR7(6)*CARRAY(IRE(6,3))-COEFI7(6)*CARRAY(IIM(6,3))
     +          + COEFR7(2)*CARRAY(IRE(6,4))-COEFI7(2)*CARRAY(IIM(6,4))
     +          + COEFR7(5)*CARRAY(IRE(6,5))-COEFI7(5)*CARRAY(IIM(6,5))
     +          + COEFR7(1)*CARRAY(IRE(6,6))-COEFI7(1)*CARRAY(IIM(6,6))
     +          + COEFR7(4)*CARRAY(IRE(6,7))-COEFI7(4)*CARRAY(IIM(6,7))
              DIFFIM = CARRAY(IIM(6,1))
     +          + COEFI7(3)*CARRAY(IRE(6,2))+COEFR7(3)*CARRAY(IIM(6,2))
     +          + COEFI7(6)*CARRAY(IRE(6,3))+COEFR7(6)*CARRAY(IIM(6,3))
     +          + COEFI7(2)*CARRAY(IRE(6,4))+COEFR7(2)*CARRAY(IIM(6,4))
     +          + COEFI7(5)*CARRAY(IRE(6,5))+COEFR7(5)*CARRAY(IIM(6,5))
     +          + COEFI7(1)*CARRAY(IRE(6,6))+COEFR7(1)*CARRAY(IIM(6,6))
     +          + COEFI7(4)*CARRAY(IRE(6,7))+COEFR7(4)*CARRAY(IIM(6,7))
              CARRAY(IRE(4,6)) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              CARRAY(IIM(4,6)) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(6,1))
     +          + COEFR7(4)*CARRAY(IRE(6,2))-COEFI7(4)*CARRAY(IIM(6,2))
     +          + COEFR7(1)*CARRAY(IRE(6,3))-COEFI7(1)*CARRAY(IIM(6,3))
     +          + COEFR7(5)*CARRAY(IRE(6,4))-COEFI7(5)*CARRAY(IIM(6,4))
     +          + COEFR7(2)*CARRAY(IRE(6,5))-COEFI7(2)*CARRAY(IIM(6,5))
     +          + COEFR7(6)*CARRAY(IRE(6,6))-COEFI7(6)*CARRAY(IIM(6,6))
     +          + COEFR7(3)*CARRAY(IRE(6,7))-COEFI7(3)*CARRAY(IIM(6,7))
              DIFFIM = CARRAY(IIM(6,1))
     +          + COEFI7(4)*CARRAY(IRE(6,2))+COEFR7(4)*CARRAY(IIM(6,2))
     +          + COEFI7(1)*CARRAY(IRE(6,3))+COEFR7(1)*CARRAY(IIM(6,3))
     +          + COEFI7(5)*CARRAY(IRE(6,4))+COEFR7(5)*CARRAY(IIM(6,4))
     +          + COEFI7(2)*CARRAY(IRE(6,5))+COEFR7(2)*CARRAY(IIM(6,5))
     +          + COEFI7(6)*CARRAY(IRE(6,6))+COEFR7(6)*CARRAY(IIM(6,6))
     +          + COEFI7(3)*CARRAY(IRE(6,7))+COEFR7(3)*CARRAY(IIM(6,7))
              CARRAY(IRE(5,6)) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              CARRAY(IIM(5,6)) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(6,1))
     +          + COEFR7(6)*CARRAY(IRE(6,2))-COEFI7(6)*CARRAY(IIM(6,2))
     +          + COEFR7(5)*CARRAY(IRE(6,3))-COEFI7(5)*CARRAY(IIM(6,3))
     +          + COEFR7(4)*CARRAY(IRE(6,4))-COEFI7(4)*CARRAY(IIM(6,4))
     +          + COEFR7(3)*CARRAY(IRE(6,5))-COEFI7(3)*CARRAY(IIM(6,5))
     +          + COEFR7(2)*CARRAY(IRE(6,6))-COEFI7(2)*CARRAY(IIM(6,6))
     +          + COEFR7(1)*CARRAY(IRE(6,7))-COEFI7(1)*CARRAY(IIM(6,7))
              DIFFIM = CARRAY(IIM(6,1))
     +          + COEFI7(6)*CARRAY(IRE(6,2))+COEFR7(6)*CARRAY(IIM(6,2))
     +          + COEFI7(5)*CARRAY(IRE(6,3))+COEFR7(5)*CARRAY(IIM(6,3))
     +          + COEFI7(4)*CARRAY(IRE(6,4))+COEFR7(4)*CARRAY(IIM(6,4))
     +          + COEFI7(3)*CARRAY(IRE(6,5))+COEFR7(3)*CARRAY(IIM(6,5))
     +          + COEFI7(2)*CARRAY(IRE(6,6))+COEFR7(2)*CARRAY(IIM(6,6))
     +          + COEFI7(1)*CARRAY(IRE(6,7))+COEFR7(1)*CARRAY(IIM(6,7))
              TEMPRE(2) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              TEMPIM(2) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              DIFFRE = CARRAY(IRE(6,1))
     +          + COEFR7(5)*CARRAY(IRE(6,2))-COEFI7(5)*CARRAY(IIM(6,2))
     +          + COEFR7(3)*CARRAY(IRE(6,3))-COEFI7(3)*CARRAY(IIM(6,3))
     +          + COEFR7(1)*CARRAY(IRE(6,4))-COEFI7(1)*CARRAY(IIM(6,4))
     +          + COEFR7(6)*CARRAY(IRE(6,5))-COEFI7(6)*CARRAY(IIM(6,5))
     +          + COEFR7(4)*CARRAY(IRE(6,6))-COEFI7(4)*CARRAY(IIM(6,6))
     +          + COEFR7(2)*CARRAY(IRE(6,7))-COEFI7(2)*CARRAY(IIM(6,7))
              DIFFIM = CARRAY(IIM(6,1))
     +          + COEFI7(5)*CARRAY(IRE(6,2))+COEFR7(5)*CARRAY(IIM(6,2))
     +          + COEFI7(3)*CARRAY(IRE(6,3))+COEFR7(3)*CARRAY(IIM(6,3))
     +          + COEFI7(1)*CARRAY(IRE(6,4))+COEFR7(1)*CARRAY(IIM(6,4))
     +          + COEFI7(6)*CARRAY(IRE(6,5))+COEFR7(6)*CARRAY(IIM(6,5))
     +          + COEFI7(4)*CARRAY(IRE(6,6))+COEFR7(4)*CARRAY(IIM(6,6))
     +          + COEFI7(2)*CARRAY(IRE(6,7))+COEFR7(2)*CARRAY(IIM(6,7))
              CARRAY(IRE(6,6)) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              CARRAY(IIM(6,6)) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              CARRAY(IRE(6,1)) = TEMPRE(5)
              CARRAY(IIM(6,1)) = TEMPIM(5)
              CARRAY(IRE(6,2)) = TEMPRE(10)
              CARRAY(IIM(6,2)) = TEMPIM(10)
              CARRAY(IRE(6,3)) = TEMPRE(13)
              CARRAY(IIM(6,3)) = TEMPIM(13)
              CARRAY(IRE(6,4)) = TEMPRE(7)
              CARRAY(IIM(6,4)) = TEMPIM(7)
              CARRAY(IRE(6,5)) = TEMPRE(1)
              CARRAY(IIM(6,5)) = TEMPIM(1)

C             SEVENTH BUTTERFLY
              CARRAY(IRE(1,7)) = CARRAY(IRE(7,1))
     +                         + CARRAY(IRE(7,2))
     +                         + CARRAY(IRE(7,3))
     +                         + CARRAY(IRE(7,4))
     +                         + CARRAY(IRE(7,5))
     +                         + CARRAY(IRE(7,6))
     +                         + CARRAY(IRE(7,7))
              CARRAY(IIM(1,7)) = CARRAY(IIM(7,1))
     +                         + CARRAY(IIM(7,2))
     +                         + CARRAY(IIM(7,3))
     +                         + CARRAY(IIM(7,4))
     +                         + CARRAY(IIM(7,5))
     +                         + CARRAY(IIM(7,6))
     +                         + CARRAY(IIM(7,7))
              DIFFRE = CARRAY(IRE(7,1))
     +          + COEFR7(1)*CARRAY(IRE(7,2))-COEFI7(1)*CARRAY(IIM(7,2))
     +          + COEFR7(2)*CARRAY(IRE(7,3))-COEFI7(2)*CARRAY(IIM(7,3))
     +          + COEFR7(3)*CARRAY(IRE(7,4))-COEFI7(3)*CARRAY(IIM(7,4))
     +          + COEFR7(4)*CARRAY(IRE(7,5))-COEFI7(4)*CARRAY(IIM(7,5))
     +          + COEFR7(5)*CARRAY(IRE(7,6))-COEFI7(5)*CARRAY(IIM(7,6))
     +          + COEFR7(6)*CARRAY(IRE(7,7))-COEFI7(6)*CARRAY(IIM(7,7))
              DIFFIM = CARRAY(IIM(7,1))
     +          + COEFI7(1)*CARRAY(IRE(7,2))+COEFR7(1)*CARRAY(IIM(7,2))
     +          + COEFI7(2)*CARRAY(IRE(7,3))+COEFR7(2)*CARRAY(IIM(7,3))
     +          + COEFI7(3)*CARRAY(IRE(7,4))+COEFR7(3)*CARRAY(IIM(7,4))
     +          + COEFI7(4)*CARRAY(IRE(7,5))+COEFR7(4)*CARRAY(IIM(7,5))
     +          + COEFI7(5)*CARRAY(IRE(7,6))+COEFR7(5)*CARRAY(IIM(7,6))
     +          + COEFI7(6)*CARRAY(IRE(7,7))+COEFR7(6)*CARRAY(IIM(7,7))
              CARRAY(IRE(2,7)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,7)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(7,1))
     +          + COEFR7(2)*CARRAY(IRE(7,2))-COEFI7(2)*CARRAY(IIM(7,2))
     +          + COEFR7(4)*CARRAY(IRE(7,3))-COEFI7(4)*CARRAY(IIM(7,3))
     +          + COEFR7(6)*CARRAY(IRE(7,4))-COEFI7(6)*CARRAY(IIM(7,4))
     +          + COEFR7(1)*CARRAY(IRE(7,5))-COEFI7(1)*CARRAY(IIM(7,5))
     +          + COEFR7(3)*CARRAY(IRE(7,6))-COEFI7(3)*CARRAY(IIM(7,6))
     +          + COEFR7(5)*CARRAY(IRE(7,7))-COEFI7(5)*CARRAY(IIM(7,7))
              DIFFIM = CARRAY(IIM(7,1))
     +          + COEFI7(2)*CARRAY(IRE(7,2))+COEFR7(2)*CARRAY(IIM(7,2))
     +          + COEFI7(4)*CARRAY(IRE(7,3))+COEFR7(4)*CARRAY(IIM(7,3))
     +          + COEFI7(6)*CARRAY(IRE(7,4))+COEFR7(6)*CARRAY(IIM(7,4))
     +          + COEFI7(1)*CARRAY(IRE(7,5))+COEFR7(1)*CARRAY(IIM(7,5))
     +          + COEFI7(3)*CARRAY(IRE(7,6))+COEFR7(3)*CARRAY(IIM(7,6))
     +          + COEFI7(5)*CARRAY(IRE(7,7))+COEFR7(5)*CARRAY(IIM(7,7))
              CARRAY(IRE(3,7)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,7)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(7,1))
     +          + COEFR7(3)*CARRAY(IRE(7,2))-COEFI7(3)*CARRAY(IIM(7,2))
     +          + COEFR7(6)*CARRAY(IRE(7,3))-COEFI7(6)*CARRAY(IIM(7,3))
     +          + COEFR7(2)*CARRAY(IRE(7,4))-COEFI7(2)*CARRAY(IIM(7,4))
     +          + COEFR7(5)*CARRAY(IRE(7,5))-COEFI7(5)*CARRAY(IIM(7,5))
     +          + COEFR7(1)*CARRAY(IRE(7,6))-COEFI7(1)*CARRAY(IIM(7,6))
     +          + COEFR7(4)*CARRAY(IRE(7,7))-COEFI7(4)*CARRAY(IIM(7,7))
              DIFFIM = CARRAY(IIM(7,1))
     +          + COEFI7(3)*CARRAY(IRE(7,2))+COEFR7(3)*CARRAY(IIM(7,2))
     +          + COEFI7(6)*CARRAY(IRE(7,3))+COEFR7(6)*CARRAY(IIM(7,3))
     +          + COEFI7(2)*CARRAY(IRE(7,4))+COEFR7(2)*CARRAY(IIM(7,4))
     +          + COEFI7(5)*CARRAY(IRE(7,5))+COEFR7(5)*CARRAY(IIM(7,5))
     +          + COEFI7(1)*CARRAY(IRE(7,6))+COEFR7(1)*CARRAY(IIM(7,6))
     +          + COEFI7(4)*CARRAY(IRE(7,7))+COEFR7(4)*CARRAY(IIM(7,7))
              CARRAY(IRE(4,7)) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              CARRAY(IIM(4,7)) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(7,1))
     +          + COEFR7(4)*CARRAY(IRE(7,2))-COEFI7(4)*CARRAY(IIM(7,2))
     +          + COEFR7(1)*CARRAY(IRE(7,3))-COEFI7(1)*CARRAY(IIM(7,3))
     +          + COEFR7(5)*CARRAY(IRE(7,4))-COEFI7(5)*CARRAY(IIM(7,4))
     +          + COEFR7(2)*CARRAY(IRE(7,5))-COEFI7(2)*CARRAY(IIM(7,5))
     +          + COEFR7(6)*CARRAY(IRE(7,6))-COEFI7(6)*CARRAY(IIM(7,6))
     +          + COEFR7(3)*CARRAY(IRE(7,7))-COEFI7(3)*CARRAY(IIM(7,7))
              DIFFIM = CARRAY(IIM(7,1))
     +          + COEFI7(4)*CARRAY(IRE(7,2))+COEFR7(4)*CARRAY(IIM(7,2))
     +          + COEFI7(1)*CARRAY(IRE(7,3))+COEFR7(1)*CARRAY(IIM(7,3))
     +          + COEFI7(5)*CARRAY(IRE(7,4))+COEFR7(5)*CARRAY(IIM(7,4))
     +          + COEFI7(2)*CARRAY(IRE(7,5))+COEFR7(2)*CARRAY(IIM(7,5))
     +          + COEFI7(6)*CARRAY(IRE(7,6))+COEFR7(6)*CARRAY(IIM(7,6))
     +          + COEFI7(3)*CARRAY(IRE(7,7))+COEFR7(3)*CARRAY(IIM(7,7))
              CARRAY(IRE(5,7)) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              CARRAY(IIM(5,7)) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(7,1))
     +          + COEFR7(5)*CARRAY(IRE(7,2))-COEFI7(5)*CARRAY(IIM(7,2))
     +          + COEFR7(3)*CARRAY(IRE(7,3))-COEFI7(3)*CARRAY(IIM(7,3))
     +          + COEFR7(1)*CARRAY(IRE(7,4))-COEFI7(1)*CARRAY(IIM(7,4))
     +          + COEFR7(6)*CARRAY(IRE(7,5))-COEFI7(6)*CARRAY(IIM(7,5))
     +          + COEFR7(4)*CARRAY(IRE(7,6))-COEFI7(4)*CARRAY(IIM(7,6))
     +          + COEFR7(2)*CARRAY(IRE(7,7))-COEFI7(2)*CARRAY(IIM(7,7))
              DIFFIM = CARRAY(IIM(7,1))
     +          + COEFI7(5)*CARRAY(IRE(7,2))+COEFR7(5)*CARRAY(IIM(7,2))
     +          + COEFI7(3)*CARRAY(IRE(7,3))+COEFR7(3)*CARRAY(IIM(7,3))
     +          + COEFI7(1)*CARRAY(IRE(7,4))+COEFR7(1)*CARRAY(IIM(7,4))
     +          + COEFI7(6)*CARRAY(IRE(7,5))+COEFR7(6)*CARRAY(IIM(7,5))
     +          + COEFI7(4)*CARRAY(IRE(7,6))+COEFR7(4)*CARRAY(IIM(7,6))
     +          + COEFI7(2)*CARRAY(IRE(7,7))+COEFR7(2)*CARRAY(IIM(7,7))
              CARRAY(IRE(6,7)) = TWDLRE(5)*DIFFRE-TWDLIM(5)*DIFFIM
              CARRAY(IIM(6,7)) = TWDLIM(5)*DIFFRE+TWDLRE(5)*DIFFIM
              DIFFRE = CARRAY(IRE(7,1))
     +          + COEFR7(6)*CARRAY(IRE(7,2))-COEFI7(6)*CARRAY(IIM(7,2))
     +          + COEFR7(5)*CARRAY(IRE(7,3))-COEFI7(5)*CARRAY(IIM(7,3))
     +          + COEFR7(4)*CARRAY(IRE(7,4))-COEFI7(4)*CARRAY(IIM(7,4))
     +          + COEFR7(3)*CARRAY(IRE(7,5))-COEFI7(3)*CARRAY(IIM(7,5))
     +          + COEFR7(2)*CARRAY(IRE(7,6))-COEFI7(2)*CARRAY(IIM(7,6))
     +          + COEFR7(1)*CARRAY(IRE(7,7))-COEFI7(1)*CARRAY(IIM(7,7))
              DIFFIM = CARRAY(IIM(7,1))
     +          + COEFI7(6)*CARRAY(IRE(7,2))+COEFR7(6)*CARRAY(IIM(7,2))
     +          + COEFI7(5)*CARRAY(IRE(7,3))+COEFR7(5)*CARRAY(IIM(7,3))
     +          + COEFI7(4)*CARRAY(IRE(7,4))+COEFR7(4)*CARRAY(IIM(7,4))
     +          + COEFI7(3)*CARRAY(IRE(7,5))+COEFR7(3)*CARRAY(IIM(7,5))
     +          + COEFI7(2)*CARRAY(IRE(7,6))+COEFR7(2)*CARRAY(IIM(7,6))
     +          + COEFI7(1)*CARRAY(IRE(7,7))+COEFR7(1)*CARRAY(IIM(7,7))
              CARRAY(IRE(7,7)) = TWDLRE(6)*DIFFRE-TWDLIM(6)*DIFFIM
              CARRAY(IIM(7,7)) = TWDLIM(6)*DIFFRE+TWDLRE(6)*DIFFIM
              CARRAY(IRE(7,1)) = TEMPRE(6)
              CARRAY(IIM(7,1)) = TEMPIM(6)
              CARRAY(IRE(7,2)) = TEMPRE(11)
              CARRAY(IIM(7,2)) = TEMPIM(11)
              CARRAY(IRE(7,3)) = TEMPRE(14)
              CARRAY(IIM(7,3)) = TEMPIM(14)
              CARRAY(IRE(7,4)) = TEMPRE(15)
              CARRAY(IIM(7,4)) = TEMPIM(15)
              CARRAY(IRE(7,5)) = TEMPRE(3)
              CARRAY(IIM(7,5)) = TEMPIM(3)
              CARRAY(IRE(7,6)) = TEMPRE(2)
              CARRAY(IIM(7,6)) = TEMPIM(2)

                DO INPT = 1,7
                  DO IFLY = 1,7
                    IBASE(IFLY,INPT) = IBASE(IFLY,INPT)+NI7
                    IF(IBASE(IFLY,INPT).GT.NX)
     +                 IBASE(IFLY,INPT) = IBASE(IFLY,INPT)-NX
                  ENDDO
                ENDDO

              ENDDO
C             LOOP OVER INCREMENTS

            ENDDO
C           LOOP OVER UNCOUPLED BLOCKS

          ENDDO
C         LOOP OVER SETS OF COUPLED BUTTERFLIES

          K1 = K1+LA
          K2 = 2*K1-1
          K3 = 3*K1-2
          K4 = 4*K1-3
          K5 = 5*K1-4
          K6 = 6*K1-5

        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = LABY7

      ENDDO
C     LOOP OVER STAGES


      RETURN
      END

