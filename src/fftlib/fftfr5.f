      SUBROUTINE FFTFR5(CARRAY,NX)

C     *************************************************************************
C
C     FFTFR5
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
C     SEPARATE INITIALISATION
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 13,3,676-686, 1992
C
C     *************************************************************************


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
      INTEGER NX
      DOUBLE PRECISION CARRAY(2*NX)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TWDLRE(4),TWDLIM(4)
      DOUBLE PRECISION TEMPRE(8),TEMPIM(8)
      DOUBLE PRECISION DIFFRE,DIFFIM
      INTEGER IBASE(5,5),JBASE(5,5),IRE(5,5),IIM(5,5)
      INTEGER LA,NXONLA,LABY5
      INTEGER K1,K2,K3,K4
      INTEGER IX,JX,KX,LX
      INTEGER IFLY,INPT
      INTEGER II,LAINC,LAINC5


C     BEGIN
C     =====


C     FIRST HALF
C     ----------
      LA = 1

C     LOOP OVER STAGES
      DO LX = 1,MPHLF5
        NXONLA = NX/LA
        JBASE(1,1) = 0
        JBASE(1,2) = NXFFTH/LA
        JBASE(1,3) = 2*NXFFTH/LA
        JBASE(1,4) = 3*NXFFTH/LA
        JBASE(1,5) = 4*NXFFTH/LA
        K1 = 1
        K2 = 1
        K3 = 1
        K4 = 1

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRM5

          TWDLRE(1) = COSTM5(K1)
          TWDLIM(1) = SINTM5(K1)
          TWDLRE(2) = COSTM5(K2)
          TWDLIM(2) = SINTM5(K2)
          TWDLRE(3) = COSTM5(K3)
          TWDLIM(3) = SINTM5(K3)
          TWDLRE(4) = COSTM5(K4)
          TWDLIM(4) = SINTM5(K4)

C         LOOP OVER BLOCKS
          DO IX = KX+1,NX,NXONLA

            IBASE(1,1) = JBASE(1,1)+IX
            IBASE(1,2) = JBASE(1,2)+IX
            IBASE(1,3) = JBASE(1,3)+IX
            IBASE(1,4) = JBASE(1,4)+IX
            IBASE(1,5) = JBASE(1,5)+IX

C           LOOP OVER INCREMENTS
            DO II = 1, INCRM5

C             LOCAL INDEXING
              IIM(1,1) = 2*IBASE(1,1)
              IRE(1,1) = IIM(1,1)-1
              IIM(1,2) = 2*IBASE(1,2)
              IRE(1,2) = IIM(1,2)-1
              IIM(1,3) = 2*IBASE(1,3)
              IRE(1,3) = IIM(1,3)-1
              IIM(1,4) = 2*IBASE(1,4)
              IRE(1,4) = IIM(1,4)-1
              IIM(1,5) = 2*IBASE(1,5)
              IRE(1,5) = IIM(1,5)-1

C             SINGLE BUTTERFLY
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(1)*CARRAY(IRE(1,2))-COEFI5(1)*CARRAY(IIM(1,2))
     +          + COEFR5(2)*CARRAY(IRE(1,3))-COEFI5(2)*CARRAY(IIM(1,3))
     +          + COEFR5(3)*CARRAY(IRE(1,4))-COEFI5(3)*CARRAY(IIM(1,4))
     +          + COEFR5(4)*CARRAY(IRE(1,5))-COEFI5(4)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(1)*CARRAY(IRE(1,2))+COEFR5(1)*CARRAY(IIM(1,2))
     +          + COEFI5(2)*CARRAY(IRE(1,3))+COEFR5(2)*CARRAY(IIM(1,3))
     +          + COEFI5(3)*CARRAY(IRE(1,4))+COEFR5(3)*CARRAY(IIM(1,4))
     +          + COEFI5(4)*CARRAY(IRE(1,5))+COEFR5(4)*CARRAY(IIM(1,5))
              TEMPRE(1) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              TEMPIM(1) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(2)*CARRAY(IRE(1,2))-COEFI5(2)*CARRAY(IIM(1,2))
     +          + COEFR5(4)*CARRAY(IRE(1,3))-COEFI5(4)*CARRAY(IIM(1,3))
     +          + COEFR5(1)*CARRAY(IRE(1,4))-COEFI5(1)*CARRAY(IIM(1,4))
     +          + COEFR5(3)*CARRAY(IRE(1,5))-COEFI5(3)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(2)*CARRAY(IRE(1,2))+COEFR5(2)*CARRAY(IIM(1,2))
     +          + COEFI5(4)*CARRAY(IRE(1,3))+COEFR5(4)*CARRAY(IIM(1,3))
     +          + COEFI5(1)*CARRAY(IRE(1,4))+COEFR5(1)*CARRAY(IIM(1,4))
     +          + COEFI5(3)*CARRAY(IRE(1,5))+COEFR5(3)*CARRAY(IIM(1,5))
              TEMPRE(2) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              TEMPIM(2) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(3)*CARRAY(IRE(1,2))-COEFI5(3)*CARRAY(IIM(1,2))
     +          + COEFR5(1)*CARRAY(IRE(1,3))-COEFI5(1)*CARRAY(IIM(1,3))
     +          + COEFR5(4)*CARRAY(IRE(1,4))-COEFI5(4)*CARRAY(IIM(1,4))
     +          + COEFR5(2)*CARRAY(IRE(1,5))-COEFI5(2)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(3)*CARRAY(IRE(1,2))+COEFR5(3)*CARRAY(IIM(1,2))
     +          + COEFI5(1)*CARRAY(IRE(1,3))+COEFR5(1)*CARRAY(IIM(1,3))
     +          + COEFI5(4)*CARRAY(IRE(1,4))+COEFR5(4)*CARRAY(IIM(1,4))
     +          + COEFI5(2)*CARRAY(IRE(1,5))+COEFR5(2)*CARRAY(IIM(1,5))
              TEMPRE(3) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(3) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(4)*CARRAY(IRE(1,2))-COEFI5(4)*CARRAY(IIM(1,2))
     +          + COEFR5(3)*CARRAY(IRE(1,3))-COEFI5(3)*CARRAY(IIM(1,3))
     +          + COEFR5(2)*CARRAY(IRE(1,4))-COEFI5(2)*CARRAY(IIM(1,4))
     +          + COEFR5(1)*CARRAY(IRE(1,5))-COEFI5(1)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(4)*CARRAY(IRE(1,2))+COEFR5(4)*CARRAY(IIM(1,2))
     +          + COEFI5(3)*CARRAY(IRE(1,3))+COEFR5(3)*CARRAY(IIM(1,3))
     +          + COEFI5(2)*CARRAY(IRE(1,4))+COEFR5(2)*CARRAY(IIM(1,4))
     +        + COEFI5(1)*CARRAY(IRE(1,5))+COEFR5(1)*CARRAY(IIM(1,5))
              TEMPRE(4) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(4) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              CARRAY(IRE(1,1)) = CARRAY(IRE(1,1))
     +                         + CARRAY(IRE(1,2))
     +                         + CARRAY(IRE(1,3))
     +                         + CARRAY(IRE(1,4))
     +                         + CARRAY(IRE(1,5))
              CARRAY(IIM(1,1)) = CARRAY(IIM(1,1))
     +                         + CARRAY(IIM(1,2))
     +                         + CARRAY(IIM(1,3))
     +                         + CARRAY(IIM(1,4))
     +                         + CARRAY(IIM(1,5))
              CARRAY(IRE(1,2)) = TEMPRE(1)
              CARRAY(IIM(1,2)) = TEMPIM(1)
              CARRAY(IRE(1,3)) = TEMPRE(2)
              CARRAY(IIM(1,3)) = TEMPIM(2)
              CARRAY(IRE(1,4)) = TEMPRE(3)
              CARRAY(IIM(1,4)) = TEMPIM(3)
              CARRAY(IRE(1,5)) = TEMPRE(4)
              CARRAY(IIM(1,5)) = TEMPIM(4)

              IBASE(1,1) = IBASE(1,1)+NI5
              IBASE(1,2) = IBASE(1,2)+NI5
              IBASE(1,3) = IBASE(1,3)+NI5
              IBASE(1,4) = IBASE(1,4)+NI5
              IBASE(1,5) = IBASE(1,5)+NI5
              IF(IBASE(1,1).GT.NX)IBASE(1,1) = IBASE(1,1)-NX
              IF(IBASE(1,2).GT.NX)IBASE(1,2) = IBASE(1,2)-NX
              IF(IBASE(1,3).GT.NX)IBASE(1,3) = IBASE(1,3)-NX
              IF(IBASE(1,4).GT.NX)IBASE(1,4) = IBASE(1,4)-NX
              IF(IBASE(1,5).GT.NX)IBASE(1,5) = IBASE(1,5)-NX

            ENDDO
C           LOOP OVER INCREMENTS

          ENDDO
C         LOOP OVER BLOCKS

          K1 = K1+LA
          K2 = 2*K1-1
          K3 = 3*K1-2
          K4 = 4*K1-3

        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = 5*LA

      ENDDO
C     LOOP OVER STAGES


C     SECOND HALF
C     -----------
C     LOOP OVER STAGES
      DO LX = MPHLF5+1,MPOWR5

        NXONLA = NX/LA
        LABY5 = LA*5
        LAINC = LA*INCRM5
        LAINC5 = LABY5*INCRM5
        DO INPT = 1,5
          JBASE(1,INPT) = (INPT-1)*NXFFTH/LA
          DO IFLY = 2,5
            JBASE(IFLY,INPT) = JBASE(IFLY-1,INPT)+LAINC
          ENDDO
        ENDDO
        K1 = 1
        K2 = 1
        K3 = 1
        K4 = 1

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRM5

          TWDLRE(1) = COSTM5(K1)
          TWDLIM(1) = SINTM5(K1)
          TWDLRE(2) = COSTM5(K2)
          TWDLIM(2) = SINTM5(K2)
          TWDLRE(3) = COSTM5(K3)
          TWDLIM(3) = SINTM5(K3)
          TWDLRE(4) = COSTM5(K4)
          TWDLIM(4) = SINTM5(K4)

C         LOOP OVER SETS OF COUPLED BUTTERFLIES
          DO JX = KX,LAINC-1,NXONLA

C           LOOP OVER UNCOUPLED BLOCKS
            DO IX = JX+1,NX,LAINC5

C             LOCAL INDEXING
              DO INPT = 1,5
                DO IFLY = 1,5
                  IBASE(IFLY,INPT) = JBASE(IFLY,INPT)+IX
                ENDDO
              ENDDO

C             LOOP OVER INCREMENTS
              DO II = 1, INCRM5

C               LOCAL INDEXING
                DO INPT = 1,5
                  DO IFLY = 1,5
                    IIM(IFLY,INPT) = 2*IBASE(IFLY,INPT)
                    IRE(IFLY,INPT) = IIM(IFLY,INPT)-1
                  ENDDO
                ENDDO

C             FIRST BUTTERFLY
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(1)*CARRAY(IRE(1,2))-COEFI5(1)*CARRAY(IIM(1,2))
     +          + COEFR5(2)*CARRAY(IRE(1,3))-COEFI5(2)*CARRAY(IIM(1,3))
     +          + COEFR5(3)*CARRAY(IRE(1,4))-COEFI5(3)*CARRAY(IIM(1,4))
     +          + COEFR5(4)*CARRAY(IRE(1,5))-COEFI5(4)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(1)*CARRAY(IRE(1,2))+COEFR5(1)*CARRAY(IIM(1,2))
     +          + COEFI5(2)*CARRAY(IRE(1,3))+COEFR5(2)*CARRAY(IIM(1,3))
     +          + COEFI5(3)*CARRAY(IRE(1,4))+COEFR5(3)*CARRAY(IIM(1,4))
     +          + COEFI5(4)*CARRAY(IRE(1,5))+COEFR5(4)*CARRAY(IIM(1,5))
              TEMPRE(1) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              TEMPIM(1) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(2)*CARRAY(IRE(1,2))-COEFI5(2)*CARRAY(IIM(1,2))
     +          + COEFR5(4)*CARRAY(IRE(1,3))-COEFI5(4)*CARRAY(IIM(1,3))
     +          + COEFR5(1)*CARRAY(IRE(1,4))-COEFI5(1)*CARRAY(IIM(1,4))
     +          + COEFR5(3)*CARRAY(IRE(1,5))-COEFI5(3)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(2)*CARRAY(IRE(1,2))+COEFR5(2)*CARRAY(IIM(1,2))
     +          + COEFI5(4)*CARRAY(IRE(1,3))+COEFR5(4)*CARRAY(IIM(1,3))
     +          + COEFI5(1)*CARRAY(IRE(1,4))+COEFR5(1)*CARRAY(IIM(1,4))
     +          + COEFI5(3)*CARRAY(IRE(1,5))+COEFR5(3)*CARRAY(IIM(1,5))
              TEMPRE(2) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              TEMPIM(2) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(3)*CARRAY(IRE(1,2))-COEFI5(3)*CARRAY(IIM(1,2))
     +          + COEFR5(1)*CARRAY(IRE(1,3))-COEFI5(1)*CARRAY(IIM(1,3))
     +          + COEFR5(4)*CARRAY(IRE(1,4))-COEFI5(4)*CARRAY(IIM(1,4))
     +          + COEFR5(2)*CARRAY(IRE(1,5))-COEFI5(2)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(3)*CARRAY(IRE(1,2))+COEFR5(3)*CARRAY(IIM(1,2))
     +          + COEFI5(1)*CARRAY(IRE(1,3))+COEFR5(1)*CARRAY(IIM(1,3))
     +          + COEFI5(4)*CARRAY(IRE(1,4))+COEFR5(4)*CARRAY(IIM(1,4))
     +          + COEFI5(2)*CARRAY(IRE(1,5))+COEFR5(2)*CARRAY(IIM(1,5))
              TEMPRE(3) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(3) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR5(4)*CARRAY(IRE(1,2))-COEFI5(4)*CARRAY(IIM(1,2))
     +          + COEFR5(3)*CARRAY(IRE(1,3))-COEFI5(3)*CARRAY(IIM(1,3))
     +          + COEFR5(2)*CARRAY(IRE(1,4))-COEFI5(2)*CARRAY(IIM(1,4))
     +          + COEFR5(1)*CARRAY(IRE(1,5))-COEFI5(1)*CARRAY(IIM(1,5))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI5(4)*CARRAY(IRE(1,2))+COEFR5(4)*CARRAY(IIM(1,2))
     +          + COEFI5(3)*CARRAY(IRE(1,3))+COEFR5(3)*CARRAY(IIM(1,3))
     +          + COEFI5(2)*CARRAY(IRE(1,4))+COEFR5(2)*CARRAY(IIM(1,4))
     +          + COEFI5(1)*CARRAY(IRE(1,5))+COEFR5(1)*CARRAY(IIM(1,5))
              TEMPRE(4) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(4) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              CARRAY(IRE(1,1)) = CARRAY(IRE(1,1))
     +                         + CARRAY(IRE(1,2))
     +                         + CARRAY(IRE(1,3))
     +                         + CARRAY(IRE(1,4))
     +                         + CARRAY(IRE(1,5))
              CARRAY(IIM(1,1)) = CARRAY(IIM(1,1))
     +                         + CARRAY(IIM(1,2))
     +                         + CARRAY(IIM(1,3))
     +                         + CARRAY(IIM(1,4))
     +                         + CARRAY(IIM(1,5))

C             SECOND BUTTERFLY
              CARRAY(IRE(1,2)) = CARRAY(IRE(2,1))
     +                         + CARRAY(IRE(2,2))
     +                         + CARRAY(IRE(2,3))
     +                         + CARRAY(IRE(2,4))
     +                         + CARRAY(IRE(2,5))
              CARRAY(IIM(1,2)) = CARRAY(IIM(2,1))
     +                         + CARRAY(IIM(2,2))
     +                         + CARRAY(IIM(2,3))
     +                         + CARRAY(IIM(2,4))
     +                         + CARRAY(IIM(2,5))
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR5(2)*CARRAY(IRE(2,2))-COEFI5(2)*CARRAY(IIM(2,2))
     +          + COEFR5(4)*CARRAY(IRE(2,3))-COEFI5(4)*CARRAY(IIM(2,3))
     +          + COEFR5(1)*CARRAY(IRE(2,4))-COEFI5(1)*CARRAY(IIM(2,4))
     +          + COEFR5(3)*CARRAY(IRE(2,5))-COEFI5(3)*CARRAY(IIM(2,5))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI5(2)*CARRAY(IRE(2,2))+COEFR5(2)*CARRAY(IIM(2,2))
     +          + COEFI5(4)*CARRAY(IRE(2,3))+COEFR5(4)*CARRAY(IIM(2,3))
     +          + COEFI5(1)*CARRAY(IRE(2,4))+COEFR5(1)*CARRAY(IIM(2,4))
     +          + COEFI5(3)*CARRAY(IRE(2,5))+COEFR5(3)*CARRAY(IIM(2,5))
              TEMPRE(5) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              TEMPIM(5) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR5(3)*CARRAY(IRE(2,2))-COEFI5(3)*CARRAY(IIM(2,2))
     +          + COEFR5(1)*CARRAY(IRE(2,3))-COEFI5(1)*CARRAY(IIM(2,3))
     +          + COEFR5(4)*CARRAY(IRE(2,4))-COEFI5(4)*CARRAY(IIM(2,4))
     +          + COEFR5(2)*CARRAY(IRE(2,5))-COEFI5(2)*CARRAY(IIM(2,5))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI5(3)*CARRAY(IRE(2,2))+COEFR5(3)*CARRAY(IIM(2,2))
     +          + COEFI5(1)*CARRAY(IRE(2,3))+COEFR5(1)*CARRAY(IIM(2,3))
     +          + COEFI5(4)*CARRAY(IRE(2,4))+COEFR5(4)*CARRAY(IIM(2,4))
     +          + COEFI5(2)*CARRAY(IRE(2,5))+COEFR5(2)*CARRAY(IIM(2,5))
              TEMPRE(6) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(6) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR5(4)*CARRAY(IRE(2,2))-COEFI5(4)*CARRAY(IIM(2,2))
     +          + COEFR5(3)*CARRAY(IRE(2,3))-COEFI5(3)*CARRAY(IIM(2,3))
     +          + COEFR5(2)*CARRAY(IRE(2,4))-COEFI5(2)*CARRAY(IIM(2,4))
     +          + COEFR5(1)*CARRAY(IRE(2,5))-COEFI5(1)*CARRAY(IIM(2,5))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI5(4)*CARRAY(IRE(2,2))+COEFR5(4)*CARRAY(IIM(2,2))
     +          + COEFI5(3)*CARRAY(IRE(2,3))+COEFR5(3)*CARRAY(IIM(2,3))
     +          + COEFI5(2)*CARRAY(IRE(2,4))+COEFR5(2)*CARRAY(IIM(2,4))
     +          + COEFI5(1)*CARRAY(IRE(2,5))+COEFR5(1)*CARRAY(IIM(2,5))
              TEMPRE(7) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(7) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR5(1)*CARRAY(IRE(2,2))-COEFI5(1)*CARRAY(IIM(2,2))
     +          + COEFR5(2)*CARRAY(IRE(2,3))-COEFI5(2)*CARRAY(IIM(2,3))
     +          + COEFR5(3)*CARRAY(IRE(2,4))-COEFI5(3)*CARRAY(IIM(2,4))
     +          + COEFR5(4)*CARRAY(IRE(2,5))-COEFI5(4)*CARRAY(IIM(2,5))
              DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI5(1)*CARRAY(IRE(2,2))+COEFR5(1)*CARRAY(IIM(2,2))
     +          + COEFI5(2)*CARRAY(IRE(2,3))+COEFR5(2)*CARRAY(IIM(2,3))
     +          + COEFI5(3)*CARRAY(IRE(2,4))+COEFR5(3)*CARRAY(IIM(2,4))
     +          + COEFI5(4)*CARRAY(IRE(2,5))+COEFR5(4)*CARRAY(IIM(2,5))
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
              CARRAY(IIM(1,3)) = CARRAY(IIM(3,1))
     +                         + CARRAY(IIM(3,2))
     +                         + CARRAY(IIM(3,3))
     +                         + CARRAY(IIM(3,4))
     +                         + CARRAY(IIM(3,5))
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR5(1)*CARRAY(IRE(3,2))-COEFI5(1)*CARRAY(IIM(3,2))
     +          + COEFR5(2)*CARRAY(IRE(3,3))-COEFI5(2)*CARRAY(IIM(3,3))
     +          + COEFR5(3)*CARRAY(IRE(3,4))-COEFI5(3)*CARRAY(IIM(3,4))
     +          + COEFR5(4)*CARRAY(IRE(3,5))-COEFI5(4)*CARRAY(IIM(3,5))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI5(1)*CARRAY(IRE(3,2))+COEFR5(1)*CARRAY(IIM(3,2))
     +          + COEFI5(2)*CARRAY(IRE(3,3))+COEFR5(2)*CARRAY(IIM(3,3))
     +          + COEFI5(3)*CARRAY(IRE(3,4))+COEFR5(3)*CARRAY(IIM(3,4))
     +          + COEFI5(4)*CARRAY(IRE(3,5))+COEFR5(4)*CARRAY(IIM(3,5))
              CARRAY(IRE(2,3)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,3)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR5(3)*CARRAY(IRE(3,2))-COEFI5(3)*CARRAY(IIM(3,2))
     +          + COEFR5(1)*CARRAY(IRE(3,3))-COEFI5(1)*CARRAY(IIM(3,3))
     +          + COEFR5(4)*CARRAY(IRE(3,4))-COEFI5(4)*CARRAY(IIM(3,4))
     +          + COEFR5(2)*CARRAY(IRE(3,5))-COEFI5(2)*CARRAY(IIM(3,5))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI5(3)*CARRAY(IRE(3,2))+COEFR5(3)*CARRAY(IIM(3,2))
     +          + COEFI5(1)*CARRAY(IRE(3,3))+COEFR5(1)*CARRAY(IIM(3,3))
     +          + COEFI5(4)*CARRAY(IRE(3,4))+COEFR5(4)*CARRAY(IIM(3,4))
     +          + COEFI5(2)*CARRAY(IRE(3,5))+COEFR5(2)*CARRAY(IIM(3,5))
              TEMPRE(1) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              TEMPIM(1) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR5(4)*CARRAY(IRE(3,2))-COEFI5(4)*CARRAY(IIM(3,2))
     +          + COEFR5(3)*CARRAY(IRE(3,3))-COEFI5(3)*CARRAY(IIM(3,3))
     +          + COEFR5(2)*CARRAY(IRE(3,4))-COEFI5(2)*CARRAY(IIM(3,4))
     +          + COEFR5(1)*CARRAY(IRE(3,5))-COEFI5(1)*CARRAY(IIM(3,5))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI5(4)*CARRAY(IRE(3,2))+COEFR5(4)*CARRAY(IIM(3,2))
     +          + COEFI5(3)*CARRAY(IRE(3,3))+COEFR5(3)*CARRAY(IIM(3,3))
     +          + COEFI5(2)*CARRAY(IRE(3,4))+COEFR5(2)*CARRAY(IIM(3,4))
     +          + COEFI5(1)*CARRAY(IRE(3,5))+COEFR5(1)*CARRAY(IIM(3,5))
              TEMPRE(8) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(8) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR5(2)*CARRAY(IRE(3,2))-COEFI5(2)*CARRAY(IIM(3,2))
     +          + COEFR5(4)*CARRAY(IRE(3,3))-COEFI5(4)*CARRAY(IIM(3,3))
     +          + COEFR5(1)*CARRAY(IRE(3,4))-COEFI5(1)*CARRAY(IIM(3,4))
     +          + COEFR5(3)*CARRAY(IRE(3,5))-COEFI5(3)*CARRAY(IIM(3,5))
              DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI5(2)*CARRAY(IRE(3,2))+COEFR5(2)*CARRAY(IIM(3,2))
     +          + COEFI5(4)*CARRAY(IRE(3,3))+COEFR5(4)*CARRAY(IIM(3,3))
     +          + COEFI5(1)*CARRAY(IRE(3,4))+COEFR5(1)*CARRAY(IIM(3,4))
     +          + COEFI5(3)*CARRAY(IRE(3,5))+COEFR5(3)*CARRAY(IIM(3,5))
              CARRAY(IRE(3,3)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,3)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              CARRAY(IRE(3,1)) = TEMPRE(2)
              CARRAY(IIM(3,1)) = TEMPIM(2)
              CARRAY(IRE(3,2)) = TEMPRE(5)
              CARRAY(IIM(3,2)) = TEMPIM(5)

C             FOURTH BUTTERFLY
              CARRAY(IRE(1,4)) = CARRAY(IRE(4,1))
     +                         + CARRAY(IRE(4,2))
     +                         + CARRAY(IRE(4,3))
     +                         + CARRAY(IRE(4,4))
     +                         + CARRAY(IRE(4,5))
              CARRAY(IIM(1,4)) = CARRAY(IIM(4,1))
     +                         + CARRAY(IIM(4,2))
     +                         + CARRAY(IIM(4,3))
     +                         + CARRAY(IIM(4,4))
     +                         + CARRAY(IIM(4,5))
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR5(1)*CARRAY(IRE(4,2))-COEFI5(1)*CARRAY(IIM(4,2))
     +          + COEFR5(2)*CARRAY(IRE(4,3))-COEFI5(2)*CARRAY(IIM(4,3))
     +          + COEFR5(3)*CARRAY(IRE(4,4))-COEFI5(3)*CARRAY(IIM(4,4))
     +          + COEFR5(4)*CARRAY(IRE(4,5))-COEFI5(4)*CARRAY(IIM(4,5))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI5(1)*CARRAY(IRE(4,2))+COEFR5(1)*CARRAY(IIM(4,2))
     +          + COEFI5(2)*CARRAY(IRE(4,3))+COEFR5(2)*CARRAY(IIM(4,3))
     +          + COEFI5(3)*CARRAY(IRE(4,4))+COEFR5(3)*CARRAY(IIM(4,4))
     +          + COEFI5(4)*CARRAY(IRE(4,5))+COEFR5(4)*CARRAY(IIM(4,5))
              CARRAY(IRE(2,4)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,4)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR5(2)*CARRAY(IRE(4,2))-COEFI5(2)*CARRAY(IIM(4,2))
     +          + COEFR5(4)*CARRAY(IRE(4,3))-COEFI5(4)*CARRAY(IIM(4,3))
     +          + COEFR5(1)*CARRAY(IRE(4,4))-COEFI5(1)*CARRAY(IIM(4,4))
     +          + COEFR5(3)*CARRAY(IRE(4,5))-COEFI5(3)*CARRAY(IIM(4,5))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI5(2)*CARRAY(IRE(4,2))+COEFR5(2)*CARRAY(IIM(4,2))
     +          + COEFI5(4)*CARRAY(IRE(4,3))+COEFR5(4)*CARRAY(IIM(4,3))
     +          + COEFI5(1)*CARRAY(IRE(4,4))+COEFR5(1)*CARRAY(IIM(4,4))
     +          + COEFI5(3)*CARRAY(IRE(4,5))+COEFR5(3)*CARRAY(IIM(4,5))
              CARRAY(IRE(3,4)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,4)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR5(4)*CARRAY(IRE(4,2))-COEFI5(4)*CARRAY(IIM(4,2))
     +          + COEFR5(3)*CARRAY(IRE(4,3))-COEFI5(3)*CARRAY(IIM(4,3))
     +          + COEFR5(2)*CARRAY(IRE(4,4))-COEFI5(2)*CARRAY(IIM(4,4))
     +          + COEFR5(1)*CARRAY(IRE(4,5))-COEFI5(1)*CARRAY(IIM(4,5))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI5(4)*CARRAY(IRE(4,2))+COEFR5(4)*CARRAY(IIM(4,2))
     +          + COEFI5(3)*CARRAY(IRE(4,3))+COEFR5(3)*CARRAY(IIM(4,3))
     +          + COEFI5(2)*CARRAY(IRE(4,4))+COEFR5(2)*CARRAY(IIM(4,4))
     +          + COEFI5(1)*CARRAY(IRE(4,5))+COEFR5(1)*CARRAY(IIM(4,5))
              TEMPRE(2) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              TEMPIM(2) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              DIFFRE = CARRAY(IRE(4,1))
     +          + COEFR5(3)*CARRAY(IRE(4,2))-COEFI5(3)*CARRAY(IIM(4,2))
     +          + COEFR5(1)*CARRAY(IRE(4,3))-COEFI5(1)*CARRAY(IIM(4,3))
     +          + COEFR5(4)*CARRAY(IRE(4,4))-COEFI5(4)*CARRAY(IIM(4,4))
     +          + COEFR5(2)*CARRAY(IRE(4,5))-COEFI5(2)*CARRAY(IIM(4,5))
              DIFFIM = CARRAY(IIM(4,1))
     +          + COEFI5(3)*CARRAY(IRE(4,2))+COEFR5(3)*CARRAY(IIM(4,2))
     +          + COEFI5(1)*CARRAY(IRE(4,3))+COEFR5(1)*CARRAY(IIM(4,3))
     +          + COEFI5(4)*CARRAY(IRE(4,4))+COEFR5(4)*CARRAY(IIM(4,4))
     +          + COEFI5(2)*CARRAY(IRE(4,5))+COEFR5(2)*CARRAY(IIM(4,5))
              CARRAY(IRE(4,4)) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              CARRAY(IIM(4,4)) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              CARRAY(IRE(4,1)) = TEMPRE(3)
              CARRAY(IIM(4,1)) = TEMPIM(3)
              CARRAY(IRE(4,2)) = TEMPRE(6)
              CARRAY(IIM(4,2)) = TEMPIM(6)
              CARRAY(IRE(4,3)) = TEMPRE(1)
              CARRAY(IIM(4,3)) = TEMPIM(1)

C             FIFTH BUTTERFLY
              CARRAY(IRE(1,5)) = CARRAY(IRE(5,1))
     +                         + CARRAY(IRE(5,2))
     +                         + CARRAY(IRE(5,3))
     +                         + CARRAY(IRE(5,4))
     +                         + CARRAY(IRE(5,5))
              CARRAY(IIM(1,5)) = CARRAY(IIM(5,1))
     +                         + CARRAY(IIM(5,2))
     +                         + CARRAY(IIM(5,3))
     +                         + CARRAY(IIM(5,4))
     +                         + CARRAY(IIM(5,5))
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR5(1)*CARRAY(IRE(5,2))-COEFI5(1)*CARRAY(IIM(5,2))
     +          + COEFR5(2)*CARRAY(IRE(5,3))-COEFI5(2)*CARRAY(IIM(5,3))
     +          + COEFR5(3)*CARRAY(IRE(5,4))-COEFI5(3)*CARRAY(IIM(5,4))
     +          + COEFR5(4)*CARRAY(IRE(5,5))-COEFI5(4)*CARRAY(IIM(5,5))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI5(1)*CARRAY(IRE(5,2))+COEFR5(1)*CARRAY(IIM(5,2))
     +          + COEFI5(2)*CARRAY(IRE(5,3))+COEFR5(2)*CARRAY(IIM(5,3))
     +          + COEFI5(3)*CARRAY(IRE(5,4))+COEFR5(3)*CARRAY(IIM(5,4))
     +          + COEFI5(4)*CARRAY(IRE(5,5))+COEFR5(4)*CARRAY(IIM(5,5))
              CARRAY(IRE(2,5)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              CARRAY(IIM(2,5)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR5(2)*CARRAY(IRE(5,2))-COEFI5(2)*CARRAY(IIM(5,2))
     +          + COEFR5(4)*CARRAY(IRE(5,3))-COEFI5(4)*CARRAY(IIM(5,3))
     +          + COEFR5(1)*CARRAY(IRE(5,4))-COEFI5(1)*CARRAY(IIM(5,4))
     +          + COEFR5(3)*CARRAY(IRE(5,5))-COEFI5(3)*CARRAY(IIM(5,5))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI5(2)*CARRAY(IRE(5,2))+COEFR5(2)*CARRAY(IIM(5,2))
     +          + COEFI5(4)*CARRAY(IRE(5,3))+COEFR5(4)*CARRAY(IIM(5,3))
     +          + COEFI5(1)*CARRAY(IRE(5,4))+COEFR5(1)*CARRAY(IIM(5,4))
     +          + COEFI5(3)*CARRAY(IRE(5,5))+COEFR5(3)*CARRAY(IIM(5,5))
              CARRAY(IRE(3,5)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              CARRAY(IIM(3,5)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR5(3)*CARRAY(IRE(5,2))-COEFI5(3)*CARRAY(IIM(5,2))
     +          + COEFR5(1)*CARRAY(IRE(5,3))-COEFI5(1)*CARRAY(IIM(5,3))
     +          + COEFR5(4)*CARRAY(IRE(5,4))-COEFI5(4)*CARRAY(IIM(5,4))
     +          + COEFR5(2)*CARRAY(IRE(5,5))-COEFI5(2)*CARRAY(IIM(5,5))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI5(3)*CARRAY(IRE(5,2))+COEFR5(3)*CARRAY(IIM(5,2))
     +          + COEFI5(1)*CARRAY(IRE(5,3))+COEFR5(1)*CARRAY(IIM(5,3))
     +          + COEFI5(4)*CARRAY(IRE(5,4))+COEFR5(4)*CARRAY(IIM(5,4))
     +          + COEFI5(2)*CARRAY(IRE(5,5))+COEFR5(2)*CARRAY(IIM(5,5))
              CARRAY(IRE(4,5)) = TWDLRE(3)*DIFFRE-TWDLIM(3)*DIFFIM
              CARRAY(IIM(4,5)) = TWDLIM(3)*DIFFRE+TWDLRE(3)*DIFFIM
              DIFFRE = CARRAY(IRE(5,1))
     +          + COEFR5(4)*CARRAY(IRE(5,2))-COEFI5(4)*CARRAY(IIM(5,2))
     +          + COEFR5(3)*CARRAY(IRE(5,3))-COEFI5(3)*CARRAY(IIM(5,3))
     +          + COEFR5(2)*CARRAY(IRE(5,4))-COEFI5(2)*CARRAY(IIM(5,4))
     +          + COEFR5(1)*CARRAY(IRE(5,5))-COEFI5(1)*CARRAY(IIM(5,5))
              DIFFIM = CARRAY(IIM(5,1))
     +          + COEFI5(4)*CARRAY(IRE(5,2))+COEFR5(4)*CARRAY(IIM(5,2))
     +          + COEFI5(3)*CARRAY(IRE(5,3))+COEFR5(3)*CARRAY(IIM(5,3))
     +          + COEFI5(2)*CARRAY(IRE(5,4))+COEFR5(2)*CARRAY(IIM(5,4))
     +          + COEFI5(1)*CARRAY(IRE(5,5))+COEFR5(1)*CARRAY(IIM(5,5))
              CARRAY(IRE(5,5)) = TWDLRE(4)*DIFFRE-TWDLIM(4)*DIFFIM
              CARRAY(IIM(5,5)) = TWDLIM(4)*DIFFRE+TWDLRE(4)*DIFFIM
              CARRAY(IRE(5,1)) = TEMPRE(4)
              CARRAY(IIM(5,1)) = TEMPIM(4)
              CARRAY(IRE(5,2)) = TEMPRE(7)
              CARRAY(IIM(5,2)) = TEMPIM(7)
              CARRAY(IRE(5,3)) = TEMPRE(8)
              CARRAY(IIM(5,3)) = TEMPIM(8)
              CARRAY(IRE(5,4)) = TEMPRE(2)
              CARRAY(IIM(5,4)) = TEMPIM(2)

                DO INPT = 1,5
                  DO IFLY = 1,5
                    IBASE(IFLY,INPT) = IBASE(IFLY,INPT)+NI5
                    IF(IBASE(IFLY,INPT).GT.NX)
     +                 IBASE(IFLY,INPT)=IBASE(IFLY,INPT)-NX
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

        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = LABY5

      ENDDO
C     LOOP OVER STAGES


      RETURN
      END

