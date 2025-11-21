      SUBROUTINE D2FDX2(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDX2
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     06-APR-2003:  RSC MODIFIED FOR SENGA2
 
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES SECOND X-DERIVATIVE OF SPECIFIED FUNCTION
C     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
C     EXPLICIT 8TH,6TH,4TH,4TH ORDER END CONDITIONS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     ARGUMENTS
C     =========
      DOUBLE PRECISION FUNCTN(NXBIGL:NXBIGR,NYBIGL:NYBIGR,NZBIGL:NZBIGR)
      DOUBLE PRECISION FDERIV(NXSIZE,NYSIZE,NZSIZE)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FDIFAP,FDIFBP,FDIFCP,FDIFDP,FDIFEP
      DOUBLE PRECISION FDIFAM,FDIFBM,FDIFCM,FDIFDM,FDIFEM
      INTEGER IC,JC,KC
      INTEGER ISTART,IFINIS
      INTEGER ICM5,ICM4,ICM3,ICM2,ICM1,ICCC,ICP1,ICP2,ICP3,ICP4,ICP5


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      ISTART = ISTAL
      IFINIS = ISTOL
      IF(NENDXL.EQ.NBOUND)ISTART = ISTAP5
      IF(NENDXR.EQ.NBOUND)IFINIS = ISTOM5

C     =========================================================================

C     INTERIOR SCHEME
C     ===============

C     TENTH ORDER EXPLICIT DIFFERENCES
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL

          ICM4 = ISTART-5
          ICM3 = ISTART-4
          ICM2 = ISTART-3
          ICM1 = ISTART-2
          ICCC = ISTART-1
          ICP1 = ISTART
          ICP2 = ISTART+1
          ICP3 = ISTART+2
          ICP4 = ISTART+3
          ICP5 = ISTART+4

          DO IC = ISTART,IFINIS

            ICM5 = ICM4
            ICM4 = ICM3
            ICM3 = ICM2
            ICM2 = ICM1
            ICM1 = ICCC
            ICCC = ICP1
            ICP1 = ICP2
            ICP2 = ICP3
            ICP3 = ICP4
            ICP4 = ICP5
            ICP5 = IC+5

            FDIFAP = FUNCTN(ICP1,JC,KC) - FUNCTN(ICCC,JC,KC) 
            FDIFAM = FUNCTN(ICCC,JC,KC) - FUNCTN(ICM1,JC,KC)
            FDIFBP = FUNCTN(ICP2,JC,KC) - FUNCTN(ICCC,JC,KC) 
            FDIFBM = FUNCTN(ICCC,JC,KC) - FUNCTN(ICM2,JC,KC) 
            FDIFCP = FUNCTN(ICP3,JC,KC) - FUNCTN(ICCC,JC,KC) 
            FDIFCM = FUNCTN(ICCC,JC,KC) - FUNCTN(ICM3,JC,KC) 
            FDIFDP = FUNCTN(ICP4,JC,KC) - FUNCTN(ICCC,JC,KC) 
            FDIFDM = FUNCTN(ICCC,JC,KC) - FUNCTN(ICM4,JC,KC) 
            FDIFEP = FUNCTN(ICP5,JC,KC) - FUNCTN(ICCC,JC,KC) 
            FDIFEM = FUNCTN(ICCC,JC,KC) - FUNCTN(ICM5,JC,KC) 

            FDERIV(IC,JC,KC) = ACOFSX*(FDIFAP-FDIFAM)
     +                       + BCOFSX*(FDIFBP-FDIFBM)
     +                       + CCOFSX*(FDIFCP-FDIFCM)
     +                       + DCOFSX*(FDIFDP-FDIFDM)
     +                       + ECOFSX*(FDIFEP-FDIFEM)

          ENDDO

        ENDDO
      ENDDO

C     =========================================================================

C     LH END
C     ======
      IF(NENDXL.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(ISTAP1,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFBP = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFCP = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFDP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFEP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAL,JC,KC) = ACFS1X*FDIFAP
     +                          + BCFS1X*FDIFBP
     +                          + CCFS1X*FDIFCP
     +                          + DCFS1X*FDIFDP
     +                          + ECFS1X*FDIFEP

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(ISTAL,JC,KC)  - FUNCTN(ISTAP1,JC,KC) 
            FDIFBP = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFCP = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFDP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFEP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDERIV(ISTAP1,JC,KC) = ACFS2X*FDIFAP
     +                           + BCFS2X*FDIFBP
     +                           + CCFS2X*FDIFCP
     +                           + DCFS2X*FDIFDP
     +                           + ECFS2X*FDIFEP

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFAM = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFBP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFBM = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP2,JC,KC) = ACFS3X*(FDIFAP-FDIFAM)
     +                           + BCFS3X*(FDIFBP-FDIFBM)

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFAM = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFBP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFBM = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFCP = FUNCTN(ISTAP6,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFCM = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP3,JC,KC) = ACFS4X*(FDIFAP-FDIFAM)
     +                           + BCFS4X*(FDIFBP-FDIFBM)
     +                           + CCFS4X*(FDIFCP-FDIFCM)

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFAM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFBP = FUNCTN(ISTAP6,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFBM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFCP = FUNCTN(ISTAP7,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFCM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFDP = FUNCTN(ISTAP8,JC,KC) - FUNCTN(ISTAP4,JC,KC) 
            FDIFDM = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP4,JC,KC) = ACFS5X*(FDIFAP-FDIFAM)
     +                           + BCFS5X*(FDIFBP-FDIFBM)
     +                           + CCFS5X*(FDIFCP-FDIFCM)
     +                           + DCFS5X*(FDIFDP-FDIFDM)

          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END
C     ======
      IF(NENDXR.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFAM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM5,JC,KC) 
            FDIFBP = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFBM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM6,JC,KC) 
            FDIFCP = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFCM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM7,JC,KC) 
            FDIFDP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM4,JC,KC) 
            FDIFDM = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM8,JC,KC) 
            FDERIV(ISTOM4,JC,KC) = ACFS5X*(FDIFAP-FDIFAM)
     +                           + BCFS5X*(FDIFBP-FDIFBM)
     +                           + CCFS5X*(FDIFCP-FDIFCM)
     +                           + DCFS5X*(FDIFDP-FDIFDM)
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFAM = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFBP = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFBM = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM5,JC,KC) 
            FDIFCP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM3,JC,KC) 
            FDIFCM = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM6,JC,KC) 
            FDERIV(ISTOM3,JC,KC) = ACFS4X*(FDIFAP-FDIFAM)
     +                           + BCFS4X*(FDIFBP-FDIFBM)
     +                           + CCFS4X*(FDIFCP-FDIFCM)
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM2,JC,KC) 
            FDIFAM = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFBP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM2,JC,KC) 
            FDIFBM = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDERIV(ISTOM2,JC,KC) = ACFS3X*(FDIFAP-FDIFAM)
     +                           + BCFS3X*(FDIFBP-FDIFBM)
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM1,JC,KC) 
            FDIFBP = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDIFCP = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDIFDP = FUNCTN(ISTOM4,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDIFEP = FUNCTN(ISTOM5,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDERIV(ISTOM1,JC,KC) = ACFS2X*FDIFAP
     +                           + BCFS2X*FDIFBP
     +                           + CCFS2X*FDIFCP
     +                           + DCFS2X*FDIFDP
     +                           + ECFS2X*FDIFEP
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(ISTOM1,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFBP = FUNCTN(ISTOM2,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFCP = FUNCTN(ISTOM3,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFDP = FUNCTN(ISTOM4,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDIFEP = FUNCTN(ISTOM5,JC,KC)  - FUNCTN(ISTOL,JC,KC) 
            FDERIV(ISTOL,JC,KC) = ACFS1X*FDIFAP
     +                          + BCFS1X*FDIFBP
     +                          + CCFS1X*FDIFCP
     +                          + DCFS1X*FDIFDP
     +                          + ECFS1X*FDIFEP

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     SCALING
C     =======
      DO KC = KSTAL, KSTOL
        DO JC = JSTAL, JSTOL
          DO IC = ISTAL, ISTOL

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDLX2

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
