      SUBROUTINE FFTFR3(CARRAY,NX)

C     *************************************************************************
C
C     FFTFR3
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
C     SEPARATE INITIALISATION
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 13,3,676-686, 1992
C
C     *************************************************************************


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
      INTEGER NX
      DOUBLE PRECISION CARRAY(2*NX)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TWDLRE(2),TWDLIM(2)
      DOUBLE PRECISION TEMPRE(3),TEMPIM(3)
      DOUBLE PRECISION DIFFRE,DIFFIM
      INTEGER IBASE(3,3),JBASE(3,3),IRE(3,3),IIM(3,3)
      INTEGER LA,K1,K2,NXONLA,LABY3
      INTEGER IX,JX,KX,LX
      INTEGER IFLY,INPT
      INTEGER II,LAINC,LAINC3


C     BEGIN
C     =====


C     FIRST HALF
C     ----------
      LA = 1

C     LOOP OVER STAGES
      DO LX = 1,MPHLF3
        NXONLA = NX/LA
        JBASE(1,1) = 0
        JBASE(1,2) = NXTHRD/LA
        JBASE(1,3) = 2*NXTHRD/LA
        K1 = 1
        K2 = 1

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRM3

          TWDLRE(1) = COSTM3(K1)
          TWDLIM(1) = SINTM3(K1)
          TWDLRE(2) = COSTM3(K2)
          TWDLIM(2) = SINTM3(K2)

C         LOOP OVER BLOCKS
          DO IX = KX+1,NX,NXONLA

            IBASE(1,1) = JBASE(1,1)+IX
            IBASE(1,2) = JBASE(1,2)+IX
            IBASE(1,3) = JBASE(1,3)+IX

C           LOOP OVER INCREMENTS
            DO II = 1, INCRM3

C             LOCAL INDEXING
              IIM(1,1) = 2*IBASE(1,1)
              IRE(1,1) = IIM(1,1)-1
              IIM(1,2) = 2*IBASE(1,2)
              IRE(1,2) = IIM(1,2)-1
              IIM(1,3) = 2*IBASE(1,3)
              IRE(1,3) = IIM(1,3)-1

C             SINGLE BUTTERFLY
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR3(1)*CARRAY(IRE(1,2))-COEFI3(1)*CARRAY(IIM(1,2))
     +          + COEFR3(2)*CARRAY(IRE(1,3))-COEFI3(2)*CARRAY(IIM(1,3))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI3(1)*CARRAY(IRE(1,2))+COEFR3(1)*CARRAY(IIM(1,2))
     +          + COEFI3(2)*CARRAY(IRE(1,3))+COEFR3(2)*CARRAY(IIM(1,3))
              TEMPRE(1) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
              TEMPIM(1) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
              DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR3(2)*CARRAY(IRE(1,2))-COEFI3(2)*CARRAY(IIM(1,2))
     +          + COEFR3(1)*CARRAY(IRE(1,3))-COEFI3(1)*CARRAY(IIM(1,3))
              DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI3(2)*CARRAY(IRE(1,2))+COEFR3(2)*CARRAY(IIM(1,2))
     +          + COEFI3(1)*CARRAY(IRE(1,3))+COEFR3(1)*CARRAY(IIM(1,3))
              TEMPRE(2) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
              TEMPIM(2) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
              CARRAY(IRE(1,1)) = CARRAY(IRE(1,1))
     +                         + CARRAY(IRE(1,2))
     +                         + CARRAY(IRE(1,3))
              CARRAY(IIM(1,1)) = CARRAY(IIM(1,1))
     +                         + CARRAY(IIM(1,2))
     +                         + CARRAY(IIM(1,3))
              CARRAY(IRE(1,2)) = TEMPRE(1)
              CARRAY(IIM(1,2)) = TEMPIM(1)
              CARRAY(IRE(1,3)) = TEMPRE(2)
              CARRAY(IIM(1,3)) = TEMPIM(2)

              IBASE(1,1) = IBASE(1,1)+NI3
              IBASE(1,2) = IBASE(1,2)+NI3
              IBASE(1,3) = IBASE(1,3)+NI3
              IF(IBASE(1,1).GT.NX)IBASE(1,1) = IBASE(1,1)-NX
              IF(IBASE(1,2).GT.NX)IBASE(1,2) = IBASE(1,2)-NX
              IF(IBASE(1,3).GT.NX)IBASE(1,3) = IBASE(1,3)-NX

            ENDDO
C           LOOP OVER INCREMENTS

          ENDDO
C         LOOP OVER BLOCKS

          K1 = K1+LA
          K2 = 2*K1-1

        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = 3*LA

      ENDDO
C     LOOP OVER STAGES


C     SECOND HALF
C     -----------

C     LOOP OVER STAGES
      DO LX = MPHLF3+1,MPOWR3

        NXONLA = NX/LA
        LABY3 = LA*3
        LAINC = LA*INCRM3
        LAINC3 = LABY3*INCRM3
        JBASE(1,1) = 0
        JBASE(1,2) = NXTHRD/LA
        JBASE(1,3) = 2*NXTHRD/LA
        JBASE(2,1) = LAINC
        JBASE(2,2) = JBASE(1,2)+LAINC
        JBASE(2,3) = JBASE(1,3)+LAINC
        JBASE(3,1) = JBASE(2,1)+LAINC
        JBASE(3,2) = JBASE(2,2)+LAINC
        JBASE(3,3) = JBASE(2,3)+LAINC
        K1 = 1
        K2 = 1

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRM3

          TWDLRE(1) = COSTM3(K1)
          TWDLIM(1) = SINTM3(K1)
          TWDLRE(2) = COSTM3(K2)
          TWDLIM(2) = SINTM3(K2)

C         LOOP OVER SETS OF COUPLED BUTTERFLIES
          DO JX = KX,LAINC-1,NXONLA

C           LOOP OVER UNCOUPLED BLOCKS
            DO IX = JX+1,NX,LAINC3

              DO INPT = 1,3
                DO IFLY = 1,3
                  IBASE(IFLY,INPT) = JBASE(IFLY,INPT)+IX
                ENDDO
              ENDDO

C             LOOP OVER INCREMENTS
              DO II = 1, INCRM3

C               LOCAL INDEXING
                DO INPT = 1,3
                  DO IFLY = 1,3
                    IIM(IFLY,INPT) = 2*IBASE(IFLY,INPT)
                    IRE(IFLY,INPT) = IIM(IFLY,INPT)-1
                  ENDDO
                ENDDO

C               FIRST BUTTERFLY
                DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR3(1)*CARRAY(IRE(1,2))-COEFI3(1)*CARRAY(IIM(1,2))
     +          + COEFR3(2)*CARRAY(IRE(1,3))-COEFI3(2)*CARRAY(IIM(1,3))
                DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI3(1)*CARRAY(IRE(1,2))+COEFR3(1)*CARRAY(IIM(1,2))
     +          + COEFI3(2)*CARRAY(IRE(1,3))+COEFR3(2)*CARRAY(IIM(1,3))
                TEMPRE(1) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
                TEMPIM(1) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
                DIFFRE = CARRAY(IRE(1,1))
     +          + COEFR3(2)*CARRAY(IRE(1,2))-COEFI3(2)*CARRAY(IIM(1,2))
     +          + COEFR3(1)*CARRAY(IRE(1,3))-COEFI3(1)*CARRAY(IIM(1,3))
                DIFFIM = CARRAY(IIM(1,1))
     +          + COEFI3(2)*CARRAY(IRE(1,2))+COEFR3(2)*CARRAY(IIM(1,2))
     +          + COEFI3(1)*CARRAY(IRE(1,3))+COEFR3(1)*CARRAY(IIM(1,3))
                TEMPRE(2) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
                TEMPIM(2) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
                CARRAY(IRE(1,1)) = CARRAY(IRE(1,1))
     +                           + CARRAY(IRE(1,2))
     +                           + CARRAY(IRE(1,3))
                CARRAY(IIM(1,1)) = CARRAY(IIM(1,1))
     +                           + CARRAY(IIM(1,2))
     +                           + CARRAY(IIM(1,3))

C               SECOND BUTTERFLY
                CARRAY(IRE(1,2)) = CARRAY(IRE(2,1))
     +                           + CARRAY(IRE(2,2))
     +                           + CARRAY(IRE(2,3))
                CARRAY(IIM(1,2)) = CARRAY(IIM(2,1))
     +                           + CARRAY(IIM(2,2))
     +                           + CARRAY(IIM(2,3))
                DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR3(2)*CARRAY(IRE(2,2))-COEFI3(2)*CARRAY(IIM(2,2))
     +          + COEFR3(1)*CARRAY(IRE(2,3))-COEFI3(1)*CARRAY(IIM(2,3))
                DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI3(2)*CARRAY(IRE(2,2))+COEFR3(2)*CARRAY(IIM(2,2))
     +          + COEFI3(1)*CARRAY(IRE(2,3))+COEFR3(1)*CARRAY(IIM(2,3))
                TEMPRE(3) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
                TEMPIM(3) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
                DIFFRE = CARRAY(IRE(2,1))
     +          + COEFR3(1)*CARRAY(IRE(2,2))-COEFI3(1)*CARRAY(IIM(2,2))
     +          + COEFR3(2)*CARRAY(IRE(2,3))-COEFI3(2)*CARRAY(IIM(2,3))
                DIFFIM = CARRAY(IIM(2,1))
     +          + COEFI3(1)*CARRAY(IRE(2,2))+COEFR3(1)*CARRAY(IIM(2,2))
     +          + COEFI3(2)*CARRAY(IRE(2,3))+COEFR3(2)*CARRAY(IIM(2,3))
                CARRAY(IRE(2,2)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
                CARRAY(IIM(2,2)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
                CARRAY(IRE(2,1)) = TEMPRE(1)
                CARRAY(IIM(2,1)) = TEMPIM(1)

C               THIRD BUTTERFLY
                CARRAY(IRE(1,3)) = CARRAY(IRE(3,1))
     +                           + CARRAY(IRE(3,2))
     +                           + CARRAY(IRE(3,3))
                CARRAY(IIM(1,3)) = CARRAY(IIM(3,1))
     +                           + CARRAY(IIM(3,2))
     +                           + CARRAY(IIM(3,3))
                DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR3(1)*CARRAY(IRE(3,2))-COEFI3(1)*CARRAY(IIM(3,2))
     +          + COEFR3(2)*CARRAY(IRE(3,3))-COEFI3(2)*CARRAY(IIM(3,3))
                DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI3(1)*CARRAY(IRE(3,2))+COEFR3(1)*CARRAY(IIM(3,2))
     +          + COEFI3(2)*CARRAY(IRE(3,3))+COEFR3(2)*CARRAY(IIM(3,3))
                CARRAY(IRE(2,3)) = TWDLRE(1)*DIFFRE-TWDLIM(1)*DIFFIM
                CARRAY(IIM(2,3)) = TWDLIM(1)*DIFFRE+TWDLRE(1)*DIFFIM
                DIFFRE = CARRAY(IRE(3,1))
     +          + COEFR3(2)*CARRAY(IRE(3,2))-COEFI3(2)*CARRAY(IIM(3,2))
     +          + COEFR3(1)*CARRAY(IRE(3,3))-COEFI3(1)*CARRAY(IIM(3,3))
                DIFFIM = CARRAY(IIM(3,1))
     +          + COEFI3(2)*CARRAY(IRE(3,2))+COEFR3(2)*CARRAY(IIM(3,2))
     +          + COEFI3(1)*CARRAY(IRE(3,3))+COEFR3(1)*CARRAY(IIM(3,3))
                CARRAY(IRE(3,3)) = TWDLRE(2)*DIFFRE-TWDLIM(2)*DIFFIM
                CARRAY(IIM(3,3)) = TWDLIM(2)*DIFFRE+TWDLRE(2)*DIFFIM
                CARRAY(IRE(3,2)) = TEMPRE(3)
                CARRAY(IIM(3,2)) = TEMPIM(3)
                CARRAY(IRE(3,1)) = TEMPRE(2)
                CARRAY(IIM(3,1)) = TEMPIM(2)

                DO INPT = 1,3
                  DO IFLY = 1,3
                    IBASE(IFLY,INPT) = IBASE(IFLY,INPT)+NI3
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

        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = LABY3

      ENDDO
C     LOOP OVER STAGES


      RETURN
      END


