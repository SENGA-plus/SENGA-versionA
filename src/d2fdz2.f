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
      INTEGER KCM5,KCM4,KCM3,KCM2,KCM1,KCCC,KCP1,KCP2,KCP3,KCP4,KCP5


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

C     INTERIOR SCHEME
C     ===============

C     TENTH ORDER EXPLICIT DIFFERENCES
      KCM4 = KSTART-5
      KCM3 = KSTART-4
      KCM2 = KSTART-3
      KCM1 = KSTART-2
      KCCC = KSTART-1
      KCP1 = KSTART
      KCP2 = KSTART+1
      KCP3 = KSTART+2
      KCP4 = KSTART+3
      KCP5 = KSTART+4

      DO KC = KSTART,KFINIS

        KCM5 = KCM4
        KCM4 = KCM3
        KCM3 = KCM2
        KCM2 = KCM1
        KCM1 = KCCC
        KCCC = KCP1
        KCP1 = KCP2
        KCP2 = KCP3
        KCP3 = KCP4
        KCP4 = KCP5
        KCP5 = KC+5

        DO JC = JSTAL,JSTOL

          DO IC = ISTAL,ISTOL

            FDIFAP = FUNCTN(IC,JC,KCP1) - FUNCTN(IC,JC,KCCC) 
            FDIFAM = FUNCTN(IC,JC,KCCC) - FUNCTN(IC,JC,KCM1)
            FDIFBP = FUNCTN(IC,JC,KCP2) - FUNCTN(IC,JC,KCCC) 
            FDIFBM = FUNCTN(IC,JC,KCCC) - FUNCTN(IC,JC,KCM2) 
            FDIFCP = FUNCTN(IC,JC,KCP3) - FUNCTN(IC,JC,KCCC) 
            FDIFCM = FUNCTN(IC,JC,KCCC) - FUNCTN(IC,JC,KCM3) 
            FDIFDP = FUNCTN(IC,JC,KCP4) - FUNCTN(IC,JC,KCCC) 
            FDIFDM = FUNCTN(IC,JC,KCCC) - FUNCTN(IC,JC,KCM4) 
            FDIFEP = FUNCTN(IC,JC,KCP5) - FUNCTN(IC,JC,KCCC) 
            FDIFEM = FUNCTN(IC,JC,KCCC) - FUNCTN(IC,JC,KCM5) 

            FDERIV(IC,JC,KC) = ACOFSZ*(FDIFAP-FDIFAM)
     +                       + BCOFSZ*(FDIFBP-FDIFBM)
     +                       + CCOFSZ*(FDIFCP-FDIFCM)
     +                       + DCOFSZ*(FDIFDP-FDIFDM)
     +                       + ECOFSZ*(FDIFEP-FDIFEM)

          ENDDO

        ENDDO
      ENDDO

C     =========================================================================

C     LH END
C     ======
      IF(NENDZL.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JC,KSTAP1) - FUNCTN(IC,JC,KSTAL) 
            FDIFBP = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAL) 
            FDIFCP = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAL) 
            FDIFDP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAL) 
            FDIFEP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAL) = ACFS1Z*FDIFAP
     +                          + BCFS1Z*FDIFBP
     +                          + CCFS1Z*FDIFCP
     +                          + DCFS1Z*FDIFDP
     +                          + ECFS1Z*FDIFEP

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JC,KSTAL)  - FUNCTN(IC,JC,KSTAP1) 
            FDIFBP = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAP1) 
            FDIFCP = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP1) 
            FDIFDP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP1) 
            FDIFEP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP1) 
            FDERIV(IC,JC,KSTAP1) = ACFS2Z*FDIFAP
     +                           + BCFS2Z*FDIFBP
     +                           + CCFS2Z*FDIFCP
     +                           + DCFS2Z*FDIFDP
     +                           + ECFS2Z*FDIFEP

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP2) 
            FDIFAM = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAP1) 
            FDIFBP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP2) 
            FDIFBM = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP2) = ACFS3Z*(FDIFAP-FDIFAM)
     +                           + BCFS3Z*(FDIFBP-FDIFBM)

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP3) 
            FDIFAM = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP2) 
            FDIFBP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP3) 
            FDIFBM = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP1) 
            FDIFCP = FUNCTN(IC,JC,KSTAP6) - FUNCTN(IC,JC,KSTAP3) 
            FDIFCM = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP3) = ACFS4Z*(FDIFAP-FDIFAM)
     +                           + BCFS4Z*(FDIFBP-FDIFBM)
     +                           + CCFS4Z*(FDIFCP-FDIFCM)

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP4) 
            FDIFAM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP3) 
            FDIFBP = FUNCTN(IC,JC,KSTAP6) - FUNCTN(IC,JC,KSTAP4) 
            FDIFBM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP2) 
            FDIFCP = FUNCTN(IC,JC,KSTAP7) - FUNCTN(IC,JC,KSTAP4) 
            FDIFCM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP1) 
            FDIFDP = FUNCTN(IC,JC,KSTAP8) - FUNCTN(IC,JC,KSTAP4) 
            FDIFDM = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP4) = ACFS5Z*(FDIFAP-FDIFAM)
     +                           + BCFS5Z*(FDIFBP-FDIFBM)
     +                           + CCFS5Z*(FDIFCP-FDIFCM)
     +                           + DCFS5Z*(FDIFDP-FDIFDM)

          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END
C     ======
      IF(NENDZR.EQ.NBOUND)THEN

C       EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
        DO JC = JSTAL,JSTOL
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
            FDERIV(IC,JC,KSTOM4) = ACFS5Z*(FDIFAP-FDIFAM)
     +                           + BCFS5Z*(FDIFBP-FDIFBM)
     +                           + CCFS5Z*(FDIFCP-FDIFCM)
     +                           + DCFS5Z*(FDIFDP-FDIFDM)
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM3) 
            FDIFAM = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM4) 
            FDIFBP = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM3) 
            FDIFBM = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM5) 
            FDIFCP = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM3) 
            FDIFCM = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM6) 
            FDERIV(IC,JC,KSTOM3) = ACFS4Z*(FDIFAP-FDIFAM)
     +                           + BCFS4Z*(FDIFBP-FDIFBM)
     +                           + CCFS4Z*(FDIFCP-FDIFCM)
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFAP = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM2) 
            FDIFAM = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM3) 
            FDIFBP = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM2) 
            FDIFBM = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM4) 
            FDERIV(IC,JC,KSTOM2) = ACFS3Z*(FDIFAP-FDIFAM)
     +                           + BCFS3Z*(FDIFBP-FDIFBM)
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFAP = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM1) 
            FDIFBP = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM1) 
            FDIFCP = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM1) 
            FDIFDP = FUNCTN(IC,JC,KSTOM4) - FUNCTN(IC,JC,KSTOM1) 
            FDIFEP = FUNCTN(IC,JC,KSTOM5) - FUNCTN(IC,JC,KSTOM1) 
            FDERIV(IC,JC,KSTOM1) = ACFS2Z*FDIFAP
     +                           + BCFS2Z*FDIFBP
     +                           + CCFS2Z*FDIFCP
     +                           + DCFS2Z*FDIFDP
     +                           + ECFS2Z*FDIFEP
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFAP = FUNCTN(IC,JC,KSTOM1)  - FUNCTN(IC,JC,KSTOL) 
            FDIFBP = FUNCTN(IC,JC,KSTOM2)  - FUNCTN(IC,JC,KSTOL) 
            FDIFCP = FUNCTN(IC,JC,KSTOM3)  - FUNCTN(IC,JC,KSTOL) 
            FDIFDP = FUNCTN(IC,JC,KSTOM4)  - FUNCTN(IC,JC,KSTOL) 
            FDIFEP = FUNCTN(IC,JC,KSTOM5)  - FUNCTN(IC,JC,KSTOL) 
            FDERIV(IC,JC,KSTOL) = ACFS1Z*FDIFAP
     +                          + BCFS1Z*FDIFBP
     +                          + CCFS1Z*FDIFCP
     +                          + DCFS1Z*FDIFDP
     +                          + ECFS1Z*FDIFEP

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     SCALING
C     =======
      DO KC = KSTAL, KSTOL
        DO JC = JSTAL, JSTOL
          DO IC = ISTAL, ISTOL

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDLZ2

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
