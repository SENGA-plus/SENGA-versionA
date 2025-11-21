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
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES SECOND Y-DERIVATIVE OF SPECIFIED FUNCTION
C     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
C     EXPLICIT 8TH,6TH,4TH,4TH ORDER END CONDITIONS
C     WITH CRAY ARCHER PERFORMANCE IMPROVEMENTS
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
      INTEGER JSTART,JFINIS


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

      DO KC = KSTAL,KSTOL

C       =======================================================================

C       INTERIOR SCHEME
C       ===============

C       TENTH ORDER EXPLICIT DIFFERENCES
        DO JC = JSTART,JFINIS
          DO IC = ISTAL,ISTOL

            FDIFAP = FUNCTN(IC,JC+1,KC) - FUNCTN(IC,JC,KC) 
            FDIFAM = FUNCTN(IC,JC-1,KC) - FUNCTN(IC,JC,KC) 
            FDIFBP = FUNCTN(IC,JC+2,KC) - FUNCTN(IC,JC,KC) 
            FDIFBM = FUNCTN(IC,JC-2,KC) - FUNCTN(IC,JC,KC) 
            FDIFCP = FUNCTN(IC,JC+3,KC) - FUNCTN(IC,JC,KC) 
            FDIFCM = FUNCTN(IC,JC-3,KC) - FUNCTN(IC,JC,KC) 
            FDIFDP = FUNCTN(IC,JC+4,KC) - FUNCTN(IC,JC,KC) 
            FDIFDM = FUNCTN(IC,JC-4,KC) - FUNCTN(IC,JC,KC) 
            FDIFEP = FUNCTN(IC,JC+5,KC) - FUNCTN(IC,JC,KC) 
            FDIFEM = FUNCTN(IC,JC-5,KC) - FUNCTN(IC,JC,KC) 

            FDERIV(IC,JC,KC) = ACOFSY*(FDIFAP+FDIFAM)
     +                       + BCOFSY*(FDIFBP+FDIFBM)
     +                       + CCOFSY*(FDIFCP+FDIFCM)
     +                       + DCOFSY*(FDIFDP+FDIFDM)
     +                       + ECOFSY*(FDIFEP+FDIFEM)

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDLY2

          ENDDO
        ENDDO

C       =======================================================================

C       LH END
C       ======
        IF(NENDYL.EQ.NBOUND)THEN

C         EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JSTAP1,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFBP = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFCP = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFDP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFEP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAL,KC) 
            FDERIV(IC,JSTAL,KC) = (ACFS1Y*FDIFAP
     +                          +  BCFS1Y*FDIFBP
     +                          +  CCFS1Y*FDIFCP
     +                          +  DCFS1Y*FDIFDP
     +                          +  ECFS1Y*FDIFEP)*OVDLY2

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JSTAL,KC)  - FUNCTN(IC,JSTAP1,KC) 
            FDIFBP = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFCP = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFDP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFEP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDERIV(IC,JSTAP1,KC) = (ACFS2Y*FDIFAP
     +                           +  BCFS2Y*FDIFBP
     +                           +  CCFS2Y*FDIFCP
     +                           +  DCFS2Y*FDIFDP
     +                           +  ECFS2Y*FDIFEP)*OVDLY2

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFAM = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAP1,KC)
            FDIFBP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFBM = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAL,KC)
            FDERIV(IC,JSTAP2,KC) = (ACFS3Y*(FDIFAP-FDIFAM)
     +                           +  BCFS3Y*(FDIFBP-FDIFBM))*OVDLY2

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFAM = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFBP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFBM = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC)
            FDIFCP = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFCM = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAL,KC)
            FDERIV(IC,JSTAP3,KC) = (ACFS4Y*(FDIFAP-FDIFAM)
     +                           +  BCFS4Y*(FDIFBP-FDIFBM)
     +                           +  CCFS4Y*(FDIFCP-FDIFCM))*OVDLY2

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFAM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP3,KC)
            FDIFBP = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFBM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFCP = FUNCTN(IC,JSTAP7,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFCM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP1,KC)
            FDIFDP = FUNCTN(IC,JSTAP8,KC) - FUNCTN(IC,JSTAP4,KC)
            FDIFDM = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC)
            FDERIV(IC,JSTAP4,KC) = (ACFS5Y*(FDIFAP-FDIFAM)
     +                           +  BCFS5Y*(FDIFBP-FDIFBM)
     +                           +  CCFS5Y*(FDIFCP-FDIFCM)
     +                           +  DCFS5Y*(FDIFDP-FDIFDM))*OVDLY2
 
          ENDDO

        ENDIF 

C       =======================================================================

C       RH END
C       ======
        IF(NENDYR.EQ.NBOUND)THEN

C         EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
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
            FDERIV(IC,JSTOM4,KC) = (ACFS5Y*(FDIFAP-FDIFAM)
     +                           +  BCFS5Y*(FDIFBP-FDIFBM)
     +                           +  CCFS5Y*(FDIFCP-FDIFCM)
     +                           +  DCFS5Y*(FDIFDP-FDIFDM))*OVDLY2
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFAM = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDIFBP = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFBM = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM5,KC) 
            FDIFCP = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM3,KC) 
            FDIFCM = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM6,KC) 
            FDERIV(IC,JSTOM3,KC) = (ACFS4Y*(FDIFAP-FDIFAM)
     +                           +  BCFS4Y*(FDIFBP-FDIFBM)
     +                           +  CCFS4Y*(FDIFCP-FDIFCM))*OVDLY2
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM2,KC) 
            FDIFAM = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFBP = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM2,KC) 
            FDIFBM = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDERIV(IC,JSTOM2,KC) = (ACFS3Y*(FDIFAP-FDIFAM)
     +                           +  BCFS3Y*(FDIFBP-FDIFBM))*OVDLY2
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM1,KC) 
            FDIFBP = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFCP = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFDP = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFEP = FUNCTN(IC,JSTOM5,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDERIV(IC,JSTOM1,KC) = (ACFS2Y*FDIFAP
     +                           +  BCFS2Y*FDIFBP
     +                           +  CCFS2Y*FDIFCP
     +                           +  DCFS2Y*FDIFDP
     +                           +  ECFS2Y*FDIFEP)*OVDLY2
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFBP = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFCP = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFDP = FUNCTN(IC,JSTOM4,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFEP = FUNCTN(IC,JSTOM5,KC) - FUNCTN(IC,JSTOL,KC) 
            FDERIV(IC,JSTOL,KC) = (ACFS1Y*FDIFAP
     +                          +  BCFS1Y*FDIFBP
     +                          +  CCFS1Y*FDIFCP
     +                          +  DCFS1Y*FDIFDP
     +                          +  ECFS1Y*FDIFEP)*OVDLY2
      
          ENDDO

        ENDIF

C       =======================================================================

      ENDDO

C     =========================================================================


      RETURN
      END
