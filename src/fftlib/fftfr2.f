      SUBROUTINE FFTFR2(CARRAY,NX)

C     *************************************************************************
C
C     FFTFR2
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     20-MAY-1999:  CREATED
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT AN FFT IN 1D
C     USING TEMPERTON SELF-SORTING PRIME FACTOR ALGORITHM
C     RADIX 2
C     WITH ROTATIONS
C     SEPARATE INITIALISATION
C
C     REFERENCES
C     ----------
C     TEMPERTON C.: SIAM J SCI STAT COMP 12,4,808-823, 1991
C
C     *************************************************************************


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
      INTEGER NX
      DOUBLE PRECISION CARRAY(2*NX)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION WTEMPR,WTEMPI,ZTEMPR,ZTEMPI,DIFFRE,DIFFIM
      INTEGER LA,KK,NXONLA,LABY2
      INTEGER IX,JX,KX,LX
      INTEGER IA,IARE,IAIM,IB,IBRE,IBIM
      INTEGER IC,ICRE,ICIM,ID,IDRE,IDIM
      INTEGER JA,JB,JC,JD,II,LAINC,LAINC2


C     BEGIN
C     =====

C     FIRST HALF
C     ----------
      LA = 1
C     LOOP OVER STAGES
      DO LX = 1,MPHLF2
        NXONLA = NX/LA
        JA = 0
        JB = NXHALF/LA
        KK = 1
C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JB-1, INCRM2
          WTEMPR = COSTM2(KK)
          WTEMPI = SINTM2(KK)
C         LOOP OVER BUTTERFLIES HAVING THE SAME TWIDDLE FACTOR
          DO IX = KX+1,NX,NXONLA
            IA = JA+IX
            IB = JB+IX
C           LOOP OVER SHORT TRANSFORMS
            DO II = 1, INCRM2
              IAIM = 2*IA
              IARE = IAIM-1
              IBIM = 2*IB
              IBRE = IBIM-1
              DIFFRE = CARRAY(IARE)-CARRAY(IBRE)
              DIFFIM = CARRAY(IAIM)-CARRAY(IBIM)
              ZTEMPR = WTEMPR*DIFFRE-WTEMPI*DIFFIM
              ZTEMPI = WTEMPI*DIFFRE+WTEMPR*DIFFIM
              CARRAY(IARE) = CARRAY(IARE)+CARRAY(IBRE)
              CARRAY(IAIM) = CARRAY(IAIM)+CARRAY(IBIM)
              CARRAY(IBRE) = ZTEMPR
              CARRAY(IBIM) = ZTEMPI
              IA = IA+NI2
              IB = IB+NI2
              IF(IA.GT.NX)IA = IA-NX
              IF(IB.GT.NX)IB = IB-NX
            ENDDO
          ENDDO
          KK = KK+LA
        ENDDO
        LA = 2*LA
      ENDDO


C     SECOND HALF
C     -----------
C     LOOP OVER STAGES
      DO LX = MPHLF2+1,MPOWR2
        NXONLA = NX/LA
        LABY2 = LA*2
        LAINC = LA*INCRM2
        LAINC2 = LABY2*INCRM2
        JA = 0
        JB = NXHALF/LA
        JC = LAINC
        JD = JB+LAINC
        KK = 1
C       LOOP OVER TWIDDLE FACTORS
        DO KX = 0, JB-1, INCRM2
          WTEMPR = COSTM2(KK)
          WTEMPI = SINTM2(KK)
C         LOOP OVER BUTTERFLY PAIRS HAVING THE SAME TWIDDLE FACTOR
          DO JX = KX,LAINC-1,NXONLA
C           LOOP OVER BASE ADDRESSES OF BUTTERFLY PAIRS
            DO IX = JX+1,NX,LAINC2
              IA = JA+IX
              IB = JB+IX
              IC = JC+IX
              ID = JD+IX
C             LOOP OVER SHORT TRANSFORMS
              DO II = 1, INCRM2
                IAIM = 2*IA
                IARE = IAIM-1
                IBIM = 2*IB
                IBRE = IBIM-1
                ICIM = 2*IC
                ICRE = ICIM-1
                IDIM = 2*ID
                IDRE = IDIM-1
                DIFFRE = CARRAY(IARE)-CARRAY(IBRE)
                DIFFIM = CARRAY(IAIM)-CARRAY(IBIM)
                ZTEMPR = WTEMPR*DIFFRE-WTEMPI*DIFFIM
                ZTEMPI = WTEMPI*DIFFRE+WTEMPR*DIFFIM
                CARRAY(IARE) = CARRAY(IARE)+CARRAY(IBRE)
                CARRAY(IAIM) = CARRAY(IAIM)+CARRAY(IBIM)
                CARRAY(IBRE) = CARRAY(ICRE)+CARRAY(IDRE)
                CARRAY(IBIM) = CARRAY(ICIM)+CARRAY(IDIM)
                DIFFRE = CARRAY(ICRE)-CARRAY(IDRE)
                DIFFIM = CARRAY(ICIM)-CARRAY(IDIM)
                CARRAY(IDRE) = WTEMPR*DIFFRE-WTEMPI*DIFFIM
                CARRAY(IDIM) = WTEMPI*DIFFRE+WTEMPR*DIFFIM
                CARRAY(ICRE) = ZTEMPR
                CARRAY(ICIM) = ZTEMPI
                IA = IA+NI2
                IB = IB+NI2
                IC = IC+NI2
                ID = ID+NI2
                IF(IA.GT.NX)IA = IA-NX
                IF(IB.GT.NX)IB = IB-NX
                IF(IC.GT.NX)IC = IC-NX
                IF(ID.GT.NX)ID = ID-NX
              ENDDO
            ENDDO
          ENDDO
          KK = KK+LA
        ENDDO
        LA = LABY2
      ENDDO


      RETURN
      END

