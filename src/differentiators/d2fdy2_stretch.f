      SUBROUTINE D2FDY2(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDY2
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     11-APR-2003:  RSC MODIFIED FOR SENGA2
C     10-NOV-2013;  RSC MODIFIED FOR MESH STRETCHING
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES SECOND Y-DERIVATIVE OF SPECIFIED FUNCTION
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
      DOUBLE PRECISION FFIRST(NXSIZE,NYSIZE,NZSIZE)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FDIFAP,FDIFBP,FDIFCP,FDIFDP,FDIFEP
      DOUBLE PRECISION FDIFAM,FDIFBM,FDIFCM,FDIFDM,FDIFEM
      DOUBLE PRECISION FDIFFA,FDIFFB,FDIFFC,FDIFFD,FDIFFE
      INTEGER IC,JC,KC
      INTEGER JSTART,JFINIS
      INTEGER JCM5,JCM4,JCM3,JCM2,JCM1,JCCC,JCP1,JCP2,JCP3,JCP4,JCP5


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      JSTART = JSTAL
      JFINIS = JSTOL
      IF(NENDYL.EQ.NBOUND)JSTART = JSTAP5
      IF(NENDYR.EQ.NBOUND)JFINIS = JSTOM5

C     =========================================================================

C     INTERIOR SCHEME
C     ===============

C     TENTH ORDER EXPLICIT DIFFERENCES
      DO KC = KSTAL,KSTOL

        JCM4 = JSTART-5
        JCM3 = JSTART-4
        JCM2 = JSTART-3
        JCM1 = JSTART-2
        JCCC = JSTART-1
        JCP1 = JSTART
        JCP2 = JSTART+1
        JCP3 = JSTART+2
        JCP4 = JSTART+3
        JCP5 = JSTART+4

        DO JC = JSTART,JFINIS

          JCM5 = JCM4
          JCM4 = JCM3
          JCM3 = JCM2
          JCM2 = JCM1
          JCM1 = JCCC
          JCCC = JCP1
          JCP1 = JCP2
          JCP2 = JCP3
          JCP3 = JCP4
          JCP4 = JCP5
          JCP5 = JC+5

          DO IC = ISTAL,ISTOL

            FDIFAP = FUNCTN(IC,JCP1,KC) - FUNCTN(IC,JCCC,KC) 
            FDIFAM = FUNCTN(IC,JCCC,KC) - FUNCTN(IC,JCM1,KC) 
            FDIFBP = FUNCTN(IC,JCP2,KC) - FUNCTN(IC,JCCC,KC) 
            FDIFBM = FUNCTN(IC,JCCC,KC) - FUNCTN(IC,JCM2,KC) 
            FDIFCP = FUNCTN(IC,JCP3,KC) - FUNCTN(IC,JCCC,KC) 
            FDIFCM = FUNCTN(IC,JCCC,KC) - FUNCTN(IC,JCM3,KC) 
            FDIFDP = FUNCTN(IC,JCP4,KC) - FUNCTN(IC,JCCC,KC) 
            FDIFDM = FUNCTN(IC,JCCC,KC) - FUNCTN(IC,JCM4,KC) 
            FDIFEP = FUNCTN(IC,JCP5,KC) - FUNCTN(IC,JCCC,KC) 
            FDIFEM = FUNCTN(IC,JCCC,KC) - FUNCTN(IC,JCM5,KC) 

            FDERIV(IC,JC,KC) = ACOFSY*(FDIFAP-FDIFAM)
     +                       + BCOFSY*(FDIFBP-FDIFBM)
     +                       + CCOFSY*(FDIFCP-FDIFCM)
     +                       + DCOFSY*(FDIFDP-FDIFDM)
     +                       + ECOFSY*(FDIFEP-FDIFEM)

          ENDDO

        ENDDO

      ENDDO

C     =========================================================================

C     LH END
C     ======
      IF(NENDYL.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JSTAP1,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFBP = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFCP = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFDP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFEP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAL,KC) 
            FDERIV(IC,JSTAL,KC) = ACFS1Y*FDIFAP
     +                          + BCFS1Y*FDIFBP
     +                          + CCFS1Y*FDIFCP
     +                          + DCFS1Y*FDIFDP
     +                          + ECFS1Y*FDIFEP

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JSTAL,KC)  - FUNCTN(IC,JSTAP1,KC) 
            FDIFBP = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFCP = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFDP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFEP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDERIV(IC,JSTAP1,KC) = ACFS2Y*FDIFAP
     +                           + BCFS2Y*FDIFBP
     +                           + CCFS2Y*FDIFCP
     +                           + DCFS2Y*FDIFDP
     +                           + ECFS2Y*FDIFEP

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFAM = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAP1,KC)
            FDIFBP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFBM = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAL,KC)
            FDERIV(IC,JSTAP2,KC) = ACFS3Y*(FDIFAP-FDIFAM)
     +                           + BCFS3Y*(FDIFBP-FDIFBM)

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFAM = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFBP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFBM = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC)
            FDIFCP = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFCM = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAL,KC)
            FDERIV(IC,JSTAP3,KC) = ACFS4Y*(FDIFAP-FDIFAM)
     +                           + BCFS4Y*(FDIFBP-FDIFBM)
     +                           + CCFS4Y*(FDIFCP-FDIFCM)

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFAM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFBP = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFBM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFCP = FUNCTN(IC,JSTAP7,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFCM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP1,KC)
            FDIFDP = FUNCTN(IC,JSTAP8,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFDM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC)
            FDERIV(IC,JSTAP4,KC) = ACFS5Y*(FDIFAP-FDIFAM)
     +                           + BCFS5Y*(FDIFBP-FDIFBM)
     +                           + CCFS5Y*(FDIFCP-FDIFCM)
     +                           + DCFS5Y*(FDIFDP-FDIFDM)
 
          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END
C     ======
      IF(NENDYR.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDIFAM = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOM5,KC) 
            FDIFBP = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDIFBM = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOM6,KC) 
            FDIFCP = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDIFCM = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOM7,KC) 
            FDIFDP = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM4,KC) 
            FDIFDM = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOM8,KC) 
            FDERIV(IC,JSTOM4,KC) = ACFS5Y*(FDIFAP-FDIFAM)
     +                           + BCFS5Y*(FDIFBP-FDIFBM)
     +                           + CCFS5Y*(FDIFCP-FDIFCM)
     +                           + DCFS5Y*(FDIFDP-FDIFDM)
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFAM = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDIFBP = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFBM = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM5,KC) 
            FDIFCP = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM3,KC) 
            FDIFCM = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM6,KC) 
            FDERIV(IC,JSTOM3,KC) = ACFS4Y*(FDIFAP-FDIFAM)
     +                           + BCFS4Y*(FDIFBP-FDIFBM)
     +                           + CCFS4Y*(FDIFCP-FDIFCM)
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM2,KC) 
            FDIFAM = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFBP = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM2,KC) 
            FDIFBM = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDERIV(IC,JSTOM2,KC) = ACFS3Y*(FDIFAP-FDIFAM)
     +                           + BCFS3Y*(FDIFBP-FDIFBM)
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM1,KC) 
            FDIFBP = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFCP = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFDP = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFEP = FUNCTN(IC,JSTOM5,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDERIV(IC,JSTOM1,KC) = ACFS2Y*FDIFAP
     +                           + BCFS2Y*FDIFBP
     +                           + CCFS2Y*FDIFCP
     +                           + DCFS2Y*FDIFDP
     +                           + ECFS2Y*FDIFEP
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFBP = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFCP = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFDP = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFEP = FUNCTN(IC,JSTOM5,KC) - FUNCTN(IC,JSTOL,KC) 
            FDERIV(IC,JSTOL,KC) = ACFS1Y*FDIFAP
     +                          + BCFS1Y*FDIFBP
     +                          + CCFS1Y*FDIFCP
     +                          + DCFS1Y*FDIFDP
     +                          + ECFS1Y*FDIFEP
      
          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     RSC 10-NOV-2013 MESH STRETCHING

C     FIRST DERIVATIVE
C     ================

C     =========================================================================

C     END CONDITIONS
C     ==============

      JSTART = JSTAL
      JFINIS = JSTOL
      IF(NENDYL.EQ.NBOUND)JSTART = JSTAP5
      IF(NENDYR.EQ.NBOUND)JFINIS = JSTOM5

C     =========================================================================

C     INTERIOR SCHEME
C     ===============

C     TENTH ORDER EXPLICIT DIFFERENCES
      DO KC = KSTAL,KSTOL

        JCM4 = JSTART-5
        JCM3 = JSTART-4
        JCM2 = JSTART-3
        JCM1 = JSTART-2
        JCCC = JSTART-1
        JCP1 = JSTART
        JCP2 = JSTART+1
        JCP3 = JSTART+2
        JCP4 = JSTART+3
        JCP5 = JSTART+4

        DO JC = JSTART,JFINIS

          JCM5 = JCM4
          JCM4 = JCM3
          JCM3 = JCM2
          JCM2 = JCM1
          JCM1 = JCCC
          JCCC = JCP1
          JCP1 = JCP2
          JCP2 = JCP3
          JCP3 = JCP4
          JCP4 = JCP5
          JCP5 = JC+5

          DO IC = ISTAL,ISTOL

            FDIFFA = FUNCTN(IC,JCP1,KC) - FUNCTN(IC,JCM1,KC) 
            FDIFFB = FUNCTN(IC,JCP2,KC) - FUNCTN(IC,JCM2,KC) 
            FDIFFC = FUNCTN(IC,JCP3,KC) - FUNCTN(IC,JCM3,KC) 
            FDIFFD = FUNCTN(IC,JCP4,KC) - FUNCTN(IC,JCM4,KC) 
            FDIFFE = FUNCTN(IC,JCP5,KC) - FUNCTN(IC,JCM5,KC) 

            FFIRST(IC,JC,KC) = ACOFFY*FDIFFA
     +                       + BCOFFY*FDIFFB
     +                       + CCOFFY*FDIFFC
     +                       + DCOFFY*FDIFFD
     +                       + ECOFFY*FDIFFE

          ENDDO

        ENDDO

      ENDDO

C     =========================================================================

C     LH END
C     ======
      IF(NENDYL.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFFA = FUNCTN(IC,JSTAP1,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFFB = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFFC = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFFD = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC) 
            FFIRST(IC,JSTAL,KC) = ACOF1Y*FDIFFA
     +                          + BCOF1Y*FDIFFB
     +                          + CCOF1Y*FDIFFC
     +                          + DCOF1Y*FDIFFD

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTAL,KC)  - FUNCTN(IC,JSTAP1,KC) 
            FDIFFB = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFC = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFD = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP1,KC) 
            FFIRST(IC,JSTAP1,KC) = ACOF2Y*FDIFFA
     +                           + BCOF2Y*FDIFFB
     +                           + CCOF2Y*FDIFFC
     +                           + DCOF2Y*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFB = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC) 
            FFIRST(IC,JSTAP2,KC) = ACOF3Y*FDIFFA
     +                           + BCOF3Y*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFFB = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFC = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAL,KC) 
            FFIRST(IC,JSTAP3,KC) = ACOF4Y*FDIFFA
     +                           + BCOF4Y*FDIFFB
     +                           + CCOF4Y*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP3,KC) 
            FDIFFB = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAP2,KC) 
            FDIFFC = FUNCTN(IC,JSTAP7,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFD = FUNCTN(IC,JSTAP8,KC) - FUNCTN(IC,JSTAL,KC) 
            FFIRST(IC,JSTAP4,KC) = ACOF5Y*FDIFFA
     +                           + BCOF5Y*FDIFFB
     +                           + CCOF5Y*FDIFFC
     +                           + DCOF5Y*FDIFFD
      
          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END
C     ======
      IF(NENDYR.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM5,KC) 
            FDIFFB = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM6,KC) 
            FDIFFC = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM7,KC) 
            FDIFFD = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM8,KC) 
            FFIRST(IC,JSTOM4,KC) = ACOF5Y*FDIFFA
     +                           + BCOF5Y*FDIFFB
     +                           + CCOF5Y*FDIFFC
     +                           + DCOF5Y*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDIFFB = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM5,KC) 
            FDIFFC = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM6,KC) 
            FFIRST(IC,JSTOM3,KC) = ACOF4Y*FDIFFA
     +                           + BCOF4Y*FDIFFB
     +                           + CCOF4Y*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFFB = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM4,KC) 
            FFIRST(IC,JSTOM2,KC) = ACOF3Y*FDIFFA
     +                           + BCOF3Y*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFFB = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM2,KC) 
            FDIFFC = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFFD = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM4,KC) 
            FFIRST(IC,JSTOM1,KC) = ACOF2Y*FDIFFA
     +                           + BCOF2Y*FDIFFB
     +                           + CCOF2Y*FDIFFC
     +                           + DCOF2Y*FDIFFD
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFFA = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFFB = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM2,KC) 
            FDIFFC = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFFD = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM4,KC) 
            FFIRST(IC,JSTOL,KC) = ACOF1Y*FDIFFA
     +                          + BCOF1Y*FDIFFB
     +                          + CCOF1Y*FDIFFC
     +                          + DCOF1Y*FDIFFD
      
          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     SCALING
C     =======
C     RSC 10-NOV-2013 MESH STRETCHING

      DO KC = KSTAL, KSTOL
        DO JC = JSTAL, JSTOL
          DO IC = ISTAL, ISTOL

C            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDLY2

            FDERIV(IC,JC,KC) = (FDERIV(IC,JC,KC)*DGDHSQ(JC)
     +                       +  FFIRST(IC,JC,KC)*D2GDH2(JC))*OVDLY2

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
