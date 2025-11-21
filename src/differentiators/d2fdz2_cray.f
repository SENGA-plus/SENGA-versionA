      SUBROUTINE D2FDZ2(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDZ2
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
C     EVALUATES SECOND Z-DERIVATIVE OF SPECIFIED FUNCTION
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
      INTEGER KSTART,KFINIS


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      KSTART = KSTAL
      KFINIS = KSTOL
      IF(NENDZL.EQ.NBOUND)KSTART = KSTAP5
      IF(NENDZR.EQ.NBOUND)KFINIS = KSTOM5

C     =========================================================================

      DO JC = JSTAL,JSTOL

C       =======================================================================

C       INTERIOR SCHEME
C       ===============

C       TENTH ORDER EXPLICIT DIFFERENCES
        DO KC = KSTART,KFINIS
          DO IC = ISTAL,ISTOL

            FDIFAP = FUNCTN(IC,JC,KC+1) - FUNCTN(IC,JC,KC) 
            FDIFAM = FUNCTN(IC,JC,KC-1) - FUNCTN(IC,JC,KC)
            FDIFBP = FUNCTN(IC,JC,KC+2) - FUNCTN(IC,JC,KC) 
            FDIFBM = FUNCTN(IC,JC,KC-2) - FUNCTN(IC,JC,KC) 
            FDIFCP = FUNCTN(IC,JC,KC+3) - FUNCTN(IC,JC,KC) 
            FDIFCM = FUNCTN(IC,JC,KC-3) - FUNCTN(IC,JC,KC) 
            FDIFDP = FUNCTN(IC,JC,KC+4) - FUNCTN(IC,JC,KC) 
            FDIFDM = FUNCTN(IC,JC,KC-4) - FUNCTN(IC,JC,KC) 
            FDIFEP = FUNCTN(IC,JC,KC+5) - FUNCTN(IC,JC,KC) 
            FDIFEM = FUNCTN(IC,JC,KC-5) - FUNCTN(IC,JC,KC) 

            FDERIV(IC,JC,KC) = ACOFSZ*(FDIFAP+FDIFAM)
     +                       + BCOFSZ*(FDIFBP+FDIFBM)
     +                       + CCOFSZ*(FDIFCP+FDIFCM)
     +                       + DCOFSZ*(FDIFDP+FDIFDM)
     +                       + ECOFSZ*(FDIFEP+FDIFEM)

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDLZ2

          ENDDO
        ENDDO

C       =======================================================================

C       LH END
C       ======
        IF(NENDZL.EQ.NBOUND)THEN

C         EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JC,KSTAP1) - FUNCTN(IC,JC,KSTAL) 
            FDIFBP = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAL) 
            FDIFCP = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAL) 
            FDIFDP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAL) 
            FDIFEP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAL) = (ACFS1Z*FDIFAP
     +                          +  BCFS1Z*FDIFBP
     +                          +  CCFS1Z*FDIFCP
     +                          +  DCFS1Z*FDIFDP
     +                          +  ECFS1Z*FDIFEP)*OVDLZ2

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JC,KSTAL)  - FUNCTN(IC,JC,KSTAP1) 
            FDIFBP = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAP1) 
            FDIFCP = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP1) 
            FDIFDP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP1) 
            FDIFEP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP1) 
            FDERIV(IC,JC,KSTAP1) = (ACFS2Z*FDIFAP
     +                           +  BCFS2Z*FDIFBP
     +                           +  CCFS2Z*FDIFCP
     +                           +  DCFS2Z*FDIFDP
     +                           +  ECFS2Z*FDIFEP)*OVDLZ2

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP2) 
            FDIFAM = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAP1) 
            FDIFBP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP2) 
            FDIFBM = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP2) = (ACFS3Z*(FDIFAP-FDIFAM)
     +                           +  BCFS3Z*(FDIFBP-FDIFBM))*OVDLZ2

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP3) 
            FDIFAM = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP2) 
            FDIFBP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP3) 
            FDIFBM = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP1) 
            FDIFCP = FUNCTN(IC,JC,KSTAP6) - FUNCTN(IC,JC,KSTAP3) 
            FDIFCM = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP3) = (ACFS4Z*(FDIFAP-FDIFAM)
     +                           +  BCFS4Z*(FDIFBP-FDIFBM)
     +                           +  CCFS4Z*(FDIFCP-FDIFCM))*OVDLZ2

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP4) 
            FDIFAM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP3) 
            FDIFBP = FUNCTN(IC,JC,KSTAP6) - FUNCTN(IC,JC,KSTAP4) 
            FDIFBM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP2) 
            FDIFCP = FUNCTN(IC,JC,KSTAP7) - FUNCTN(IC,JC,KSTAP4) 
            FDIFCM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP1) 
            FDIFDP = FUNCTN(IC,JC,KSTAP8) - FUNCTN(IC,JC,KSTAP4) 
            FDIFDM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP4) = (ACFS5Z*(FDIFAP-FDIFAM)
     +                           +  BCFS5Z*(FDIFBP-FDIFBM)
     +                           +  CCFS5Z*(FDIFCP-FDIFCM)
     +                           +  DCFS5Z*(FDIFDP-FDIFDM))*OVDLZ2

          ENDDO

        ENDIF 

C       =======================================================================

C       RH END
C       ======
        IF(NENDZR.EQ.NBOUND)THEN

C         EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM4) 
            FDIFAM = FUNCTN(IC,JC,KSTOM4) - FUNCTN(IC,JC,KSTOM5) 
            FDIFBP = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM4) 
            FDIFBM = FUNCTN(IC,JC,KSTOM4) - FUNCTN(IC,JC,KSTOM6) 
            FDIFCP = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM4) 
            FDIFCM = FUNCTN(IC,JC,KSTOM4) - FUNCTN(IC,JC,KSTOM7) 
            FDIFDP = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM4) 
            FDIFDM = FUNCTN(IC,JC,KSTOM4) - FUNCTN(IC,JC,KSTOM8) 
            FDERIV(IC,JC,KSTOM4) = (ACFS5Z*(FDIFAP-FDIFAM)
     +                           +  BCFS5Z*(FDIFBP-FDIFBM)
     +                           +  CCFS5Z*(FDIFCP-FDIFCM)
     +                           +  DCFS5Z*(FDIFDP-FDIFDM))*OVDLZ2
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM3) 
            FDIFAM = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM4) 
            FDIFBP = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM3) 
            FDIFBM = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM5) 
            FDIFCP = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM3) 
            FDIFCM = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM6) 
            FDERIV(IC,JC,KSTOM3) = (ACFS4Z*(FDIFAP-FDIFAM)
     +                           +  BCFS4Z*(FDIFBP-FDIFBM)
     +                           +  CCFS4Z*(FDIFCP-FDIFCM))*OVDLZ2
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM2) 
            FDIFAM = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM3) 
            FDIFBP = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM2) 
            FDIFBM = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM4) 
            FDERIV(IC,JC,KSTOM2) = (ACFS3Z*(FDIFAP-FDIFAM)
     +                           +  BCFS3Z*(FDIFBP-FDIFBM))*OVDLZ2
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM1) 
            FDIFBP = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM1) 
            FDIFCP = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM1) 
            FDIFDP = FUNCTN(IC,JC,KSTOM4) - FUNCTN(IC,JC,KSTOM1) 
            FDIFEP = FUNCTN(IC,JC,KSTOM5) - FUNCTN(IC,JC,KSTOM1) 
            FDERIV(IC,JC,KSTOM1) = (ACFS2Z*FDIFAP
     +                           +  BCFS2Z*FDIFBP
     +                           +  CCFS2Z*FDIFCP
     +                           +  DCFS2Z*FDIFDP
     +                           +  ECFS2Z*FDIFEP)*OVDLZ2
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JC,KSTOM1)  - FUNCTN(IC,JC,KSTOL) 
            FDIFBP = FUNCTN(IC,JC,KSTOM2)  - FUNCTN(IC,JC,KSTOL) 
            FDIFCP = FUNCTN(IC,JC,KSTOM3)  - FUNCTN(IC,JC,KSTOL) 
            FDIFDP = FUNCTN(IC,JC,KSTOM4)  - FUNCTN(IC,JC,KSTOL) 
            FDIFEP = FUNCTN(IC,JC,KSTOM5)  - FUNCTN(IC,JC,KSTOL) 
            FDERIV(IC,JC,KSTOL) = (ACFS1Z*FDIFAP
     +                          +  BCFS1Z*FDIFBP
     +                          +  CCFS1Z*FDIFCP
     +                          +  DCFS1Z*FDIFDP
     +                          +  ECFS1Z*FDIFEP)*OVDLZ2

          ENDDO

        ENDIF

C       =======================================================================

      ENDDO

C     =========================================================================


      RETURN
      END
