      SUBROUTINE FFTFRA(CARRAY,NX)

C     *************************************************************************
C
C     FFTFRA
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
C     ARBITRARY RADIX
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
C     FFACOM-------------------------------------------------------------------

      INTEGER NFTMAX
      PARAMETER(NFTMAX=1024)

      INTEGER NRDXMX,NRDXM1
      PARAMETER(NRDXMX=11,NRDXM1=NRDXMX-1)

      DOUBLE PRECISION COSTMA(NFTMAX),SINTMA(NFTMAX)
      DOUBLE PRECISION COEFRA(NRDXMX),COEFIA(NRDXMX)

      INTEGER NXFRAC,MPOWRA,MPHLFA,NIA,INCRMA,NRADXA,NRADM1

      COMMON/FFACOM/COSTMA,SINTMA,
     +              COEFRA,COEFIA,
     +              NXFRAC,MPOWRA,MPHLFA,NIA,INCRMA,NRADXA,NRADM1

C     FFACOM-------------------------------------------------------------------


C     ARGUMENTS
C     =========
      INTEGER NX
      DOUBLE PRECISION CARRAY(2*NX)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TWDLRE(NRDXM1),TWDLIM(NRDXM1)
      DOUBLE PRECISION TTMPRE(NRDXMX,NRDXMX),TTMPIM(NRDXMX,NRDXMX)
      DOUBLE PRECISION DIFFRE,DIFFIM
      INTEGER IBASE(NRDXMX,NRDXMX),JBASE(NRDXMX,NRDXMX)
      INTEGER IRE(NRDXMX,NRDXMX),IIM(NRDXMX,NRDXMX)
      INTEGER LA,NXONLA,LABYR
      INTEGER KK(NRDXM1)
      INTEGER IX,JX,KX,LX
      INTEGER IFLY,INPT,IC,JC,INDX,INCR
      INTEGER II,LAINC,LAINCA


C     BEGIN
C     =====

C     FIRST HALF
C     ----------
      LA = 1

C     LOOP OVER STAGES
      DO LX = 1,MPHLFA

        NXONLA = NX/LA
        DO IC = 1,NRADXA
          JBASE(1,IC) = (IC-1)*NXFRAC/LA
        ENDDO
        DO IC = 1,NRADM1
          KK(IC) = 1
        ENDDO

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRMA

          DO IC = 1,NRADM1
            TWDLRE(IC) = COSTMA(KK(IC))
            TWDLIM(IC) = SINTMA(KK(IC))
          ENDDO

C         LOOP OVER BLOCKS
          DO IX = KX+1,NX,NXONLA

C           LOCAL INDEXING
            DO INPT = 1,NRADXA
              IBASE(1,INPT) = JBASE(1,INPT)+IX
            ENDDO

C           LOOP OVER INCREMENTS
            DO II = 1, INCRMA

C             LOCAL INDEXING
              DO INPT = 1,NRADXA
                IIM(1,INPT) = 2*IBASE(1,INPT)
                IRE(1,INPT) = IIM(1,INPT)-1
              ENDDO

C             SINGLE BUTTERFLY
C             FIRST OUTPUT
              DIFFRE = CARRAY(IRE(1,1))
              DIFFIM = CARRAY(IIM(1,1))
              DO INPT = 2,NRADXA
                DIFFRE = DIFFRE + CARRAY(IRE(1,INPT))
                DIFFIM = DIFFIM + CARRAY(IIM(1,INPT))
              ENDDO
              TTMPRE(1,1) = DIFFRE
              TTMPIM(1,1) = DIFFIM
C             REST OF OUTPUTS
              DO IC = 2,NRADXA
                DIFFRE = CARRAY(IRE(1,1))
                DIFFIM = CARRAY(IIM(1,1))
                INDX = 0
                INCR = IC-1
                DO INPT = 2,NRADXA
                  INDX = MOD((INDX+INCR),NRADXA)
                  DIFFRE = DIFFRE + COEFRA(INDX)*CARRAY(IRE(1,INPT))
     +                            - COEFIA(INDX)*CARRAY(IIM(1,INPT))
                  DIFFIM = DIFFIM + COEFIA(INDX)*CARRAY(IRE(1,INPT))
     +                            + COEFRA(INDX)*CARRAY(IIM(1,INPT))
                ENDDO
                TTMPRE(1,IC) = TWDLRE(INCR)*DIFFRE-TWDLIM(INCR)*DIFFIM
                TTMPIM(1,IC) = TWDLIM(INCR)*DIFFRE+TWDLRE(INCR)*DIFFIM
              ENDDO
C             RESTORE FROM TEMPORARY STORAGE
              DO IC = 1,NRADXA
                CARRAY(IRE(1,IC)) = TTMPRE(1,IC)
                CARRAY(IIM(1,IC)) = TTMPIM(1,IC)
              ENDDO

              DO INPT = 1,NRADXA
                IBASE(1,INPT) = IBASE(1,INPT)+NIA
                IF(IBASE(1,INPT).GT.NX)
     +             IBASE(1,INPT) = IBASE(1,INPT)-NX
              ENDDO

            ENDDO
C           LOOP OVER INCREMENTS

          ENDDO
C         LOOP OVER BLOCKS

          KK(1) = KK(1)+LA
          DO IC = 2,NRADM1
            KK(IC) = IC*(KK(1)-1)+1
          ENDDO

        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = NRADXA*LA

      ENDDO
C     LOOP OVER STAGES


C     SECOND HALF
C     -----------
C     LOOP OVER STAGES
      DO LX = MPHLFA+1,MPOWRA

        NXONLA = NX/LA
        LABYR = LA*NRADXA
        LAINC = LA*INCRMA
        LAINCA = LABYR*INCRMA
        DO INPT = 1,NRADXA
          JBASE(1,INPT) = (INPT-1)*NXFRAC/LA
          DO IFLY = 2,NRADXA
            JBASE(IFLY,INPT) = JBASE(IFLY-1,INPT)+LAINC
          ENDDO
        ENDDO
        DO IC = 1,NRADM1
          KK(IC) = 1
        ENDDO

C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JBASE(1,2)-1, INCRMA

          DO IC = 1,NRADM1
            TWDLRE(IC) = COSTMA(KK(IC))
            TWDLIM(IC) = SINTMA(KK(IC))
          ENDDO

C         LOOP OVER SETS OF COUPLED BUTTERFLIES
          DO JX = KX,LAINC-1,NXONLA

C           LOOP OVER UNCOUPLED BLOCKS
            DO IX = JX+1,NX,LAINCA

C             LOCAL INDEXING
              DO INPT = 1,NRADXA
                DO IFLY = 1,NRADXA
                  IBASE(IFLY,INPT) = JBASE(IFLY,INPT)+IX
                ENDDO
              ENDDO

C             LOOP OVER INCREMENTS
              DO II = 1, INCRMA

C               LOCAL INDEXING
                DO INPT = 1,NRADXA
                  DO IFLY = 1,NRADXA
                    IIM(IFLY,INPT) = 2*IBASE(IFLY,INPT)
                    IRE(IFLY,INPT) = IIM(IFLY,INPT)-1
                  ENDDO
                ENDDO

C               COUPLED BUTTERFLIES
                DO JC = 1,NRADXA

C                 FIRST OUTPUT
                  DIFFRE = CARRAY(IRE(JC,1))
                  DIFFIM = CARRAY(IIM(JC,1))
                  DO INPT = 2,NRADXA
                    DIFFRE = DIFFRE + CARRAY(IRE(JC,INPT))
                    DIFFIM = DIFFIM + CARRAY(IIM(JC,INPT))
                  ENDDO
                  TTMPRE(JC,1) = DIFFRE
                  TTMPIM(JC,1) = DIFFIM
C                 REST OF OUTPUTS
                  DO IC = 2,NRADXA
                    DIFFRE = CARRAY(IRE(JC,1))
                    DIFFIM = CARRAY(IIM(JC,1))
                    INDX = 0
                    INCR = IC-1
                    DO INPT = 2,NRADXA
                      INDX = MOD((INDX+INCR),NRADXA)
                      DIFFRE = DIFFRE
     +                       + COEFRA(INDX)*CARRAY(IRE(JC,INPT))
     +                       - COEFIA(INDX)*CARRAY(IIM(JC,INPT))
                      DIFFIM = DIFFIM
     +                       + COEFIA(INDX)*CARRAY(IRE(JC,INPT))
     +                       + COEFRA(INDX)*CARRAY(IIM(JC,INPT))
                    ENDDO
                  TTMPRE(JC,IC)=TWDLRE(INCR)*DIFFRE-TWDLIM(INCR)*DIFFIM
                  TTMPIM(JC,IC)=TWDLIM(INCR)*DIFFRE+TWDLRE(INCR)*DIFFIM
                  ENDDO

                ENDDO

C               RESTORE FROM TEMPORARY STORAGE
C               TRANSPOSED
                DO JC = 1,NRADXA
                  DO IC = 1,NRADXA
                    CARRAY(IRE(IC,JC)) = TTMPRE(JC,IC)
                    CARRAY(IIM(IC,JC)) = TTMPIM(JC,IC)
                  ENDDO
                ENDDO

                DO INPT = 1,NRADXA
                  DO IFLY = 1,NRADXA
                    IBASE(IFLY,INPT) = IBASE(IFLY,INPT)+NIA
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

          KK(1) = KK(1)+LA
          DO IC = 2,NRADM1
            KK(IC) = IC*(KK(1)-1)+1
          ENDDO


        ENDDO
C       LOOP OVER TWIDDLE FACTORS

        LA = LABYR

      ENDDO
C     LOOP OVER STAGES


      RETURN
      END
