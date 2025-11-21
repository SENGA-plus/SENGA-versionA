      SUBROUTINE D2FDXZ(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDXZ
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     15-APR-2003:  RSC MODIFIED FOR SENGA2
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES SECOND XZ-DERIVATIVE OF SPECIFIED FUNCTION
C     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
C     EXPLICIT 8TH,6TH,4TH,4TH,4TH COMPATIBLE ORDER END CONDITIONS
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
      DOUBLE PRECISION FSTORA(3,3),FSTORB(3,3),FSTORC(2,2)
      DOUBLE PRECISION FDIFFA,FDIFFB,FDIFFC,FDIFFD,FDIFFE
      INTEGER IC,JC,KC
      INTEGER IS,KS,ISM1,KSM1
      INTEGER ISTART,IFINIS,KSTART,KFINIS
      INTEGER ICM4,ICM3,ICM2,ICM1,ICCC,ICP1,ICP2,ICP3,ICP4
      INTEGER KCM4,KCM3,KCM2,KCM1,KCP1,KCP2,KCP3,KCP4


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      ISTART = ISTAL
      IFINIS = ISTOL
      KSTART = KSTAL
      KFINIS = KSTOL
      IF(NENDXL.EQ.NBOUND)ISTART = ISTAP5
      IF(NENDXR.EQ.NBOUND)IFINIS = ISTOM5
      IF(NENDZL.EQ.NBOUND)KSTART = KSTAP5
      IF(NENDZR.EQ.NBOUND)KFINIS = KSTOM5

C     =========================================================================

      DO KC = KSTART,KFINIS

C       =======================================================================

C       INTERIOR SCHEME
C       ===============

C       TENTH ORDER EXPLICIT DIFFERENCES
        DO JC = JSTAL,JSTOL
          DO IC = ISTART,IFINIS
             
            FDIFFA = FUNCTN(IC+1,JC,KC+1) - FUNCTN(IC+1,JC,KC-1) 
     +             - FUNCTN(IC-1,JC,KC+1) + FUNCTN(IC-1,JC,KC-1) 
            FDIFFB = FUNCTN(IC+2,JC,KC+2) - FUNCTN(IC+2,JC,KC-2) 
     +             - FUNCTN(IC-2,JC,KC+2) + FUNCTN(IC-2,JC,KC-2) 
            FDIFFC = FUNCTN(IC+3,JC,KC+3) - FUNCTN(IC+3,JC,KC-3) 
     +             - FUNCTN(IC-3,JC,KC+3) + FUNCTN(IC-3,JC,KC-3) 
            FDIFFD = FUNCTN(IC+4,JC,KC+4) - FUNCTN(IC+4,JC,KC-4) 
     +             - FUNCTN(IC-4,JC,KC+4) + FUNCTN(IC-4,JC,KC-4) 
            FDIFFE = FUNCTN(IC+5,JC,KC+5) - FUNCTN(IC+5,JC,KC-5) 
     +             - FUNCTN(IC-5,JC,KC+5) + FUNCTN(IC-5,JC,KC-5) 

            FDERIV(IC,JC,KC) = ACOFXZ*FDIFFA
     +                       + BCOFXZ*FDIFFB
     +                       + CCOFXZ*FDIFFC
     +                       + DCOFXZ*FDIFFD
     +                       + ECOFXZ*FDIFFE

          ENDDO
        ENDDO

C       =======================================================================

C       LH END X-DIRECTION
C       ==================
        IF(NENDXL.EQ.NBOUND)THEN

C         TAKE SECOND XZ-DERIVATIVE IN X-LEFT INNER HALO
C         EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
          DO JC = JSTAL,JSTOL

C           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTAP1,JC,KC+1) - FUNCTN(ISTAP1,JC,KC-1)
     +                - FUNCTN(ISTAL,JC,KC+1)  + FUNCTN(ISTAL,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAP1,JC,KC+2) - FUNCTN(ISTAP1,JC,KC-2)
     +                - FUNCTN(ISTAL,JC,KC+2)  + FUNCTN(ISTAL,JC,KC-2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTAP2,JC,KC+1) - FUNCTN(ISTAP2,JC,KC-1)
     +                - FUNCTN(ISTAL,JC,KC+1)  + FUNCTN(ISTAL,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAP2,JC,KC+2) - FUNCTN(ISTAP2,JC,KC-2)
     +                - FUNCTN(ISTAL,JC,KC+2)  + FUNCTN(ISTAL,JC,KC-2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTAP3,JC,KC+1) - FUNCTN(ISTAP3,JC,KC-1)
     +                - FUNCTN(ISTAL,JC,KC+1)  + FUNCTN(ISTAL,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAP3,JC,KC+2) - FUNCTN(ISTAP3,JC,KC-2)
     +                - FUNCTN(ISTAL,JC,KC+2)  + FUNCTN(ISTAL,JC,KC-2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTAP4,JC,KC+1) - FUNCTN(ISTAP4,JC,KC-1)
     +                - FUNCTN(ISTAL,JC,KC+1)  + FUNCTN(ISTAL,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAP4,JC,KC+2) - FUNCTN(ISTAP4,JC,KC-2)
     +                - FUNCTN(ISTAL,JC,KC+2)  + FUNCTN(ISTAL,JC,KC-2))
            FDERIV(ISTAL,JC,KC) = ACF1XZ*FDIFFA
     +                          + BCF1XZ*FDIFFB
     +                          + CCF1XZ*FDIFFC
     +                          + DCF1XZ*FDIFFD
     
C           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTAL,JC,KC+1)  - FUNCTN(ISTAL,JC,KC-1)
     +                - FUNCTN(ISTAP1,JC,KC+1) + FUNCTN(ISTAP1,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAL,JC,KC+2)  - FUNCTN(ISTAL,JC,KC-2)
     +                - FUNCTN(ISTAP1,JC,KC+2) + FUNCTN(ISTAP1,JC,KC-2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTAP2,JC,KC+1) - FUNCTN(ISTAP2,JC,KC-1)
     +                - FUNCTN(ISTAP1,JC,KC+1) + FUNCTN(ISTAP1,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAP2,JC,KC+2) - FUNCTN(ISTAP2,JC,KC-2)
     +                - FUNCTN(ISTAP1,JC,KC+2) + FUNCTN(ISTAP1,JC,KC-2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTAP3,JC,KC+1) - FUNCTN(ISTAP3,JC,KC-1)
     +                - FUNCTN(ISTAP1,JC,KC+1) + FUNCTN(ISTAP1,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAP3,JC,KC+2) - FUNCTN(ISTAP3,JC,KC-2)
     +                - FUNCTN(ISTAP1,JC,KC+2) + FUNCTN(ISTAP1,JC,KC-2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTAP4,JC,KC+1) - FUNCTN(ISTAP4,JC,KC-1)
     +                - FUNCTN(ISTAP1,JC,KC+1) + FUNCTN(ISTAP1,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTAP4,JC,KC+2) - FUNCTN(ISTAP4,JC,KC-2)
     +                - FUNCTN(ISTAP1,JC,KC+2) + FUNCTN(ISTAP1,JC,KC-2))
            FDERIV(ISTAP1,JC,KC) = ACF2XZ*FDIFFA
     +                           + BCF2XZ*FDIFFB
     +                           + CCF2XZ*FDIFFC
     +                           + DCF2XZ*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP3,JC,KC+1) - FUNCTN(ISTAP3,JC,KC-1) 
     +             - FUNCTN(ISTAP1,JC,KC+1) + FUNCTN(ISTAP1,JC,KC-1) 
            FDIFFB = FUNCTN(ISTAP4,JC,KC+2) - FUNCTN(ISTAP4,JC,KC-2) 
     +             - FUNCTN(ISTAL,JC,KC+2)  + FUNCTN(ISTAL,JC,KC-2) 
            FDERIV(ISTAP2,JC,KC) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP4,JC,KC+1) - FUNCTN(ISTAP4,JC,KC-1) 
     +             - FUNCTN(ISTAP2,JC,KC+1) + FUNCTN(ISTAP2,JC,KC-1) 
            FDIFFB = FUNCTN(ISTAP5,JC,KC+2) - FUNCTN(ISTAP5,JC,KC-2) 
     +             - FUNCTN(ISTAP1,JC,KC+2) + FUNCTN(ISTAP1,JC,KC-2) 
            FDIFFC = FUNCTN(ISTAP6,JC,KC+3) - FUNCTN(ISTAP6,JC,KC-3) 
     +             - FUNCTN(ISTAL,JC,KC+3)  + FUNCTN(ISTAL,JC,KC-3) 
            FDERIV(ISTAP3,JC,KC) = ACF4XZ*FDIFFA
     +                           + BCF4XZ*FDIFFB
     +                           + CCF4XZ*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP5,JC,KC+1) - FUNCTN(ISTAP5,JC,KC-1) 
     +             - FUNCTN(ISTAP3,JC,KC+1) + FUNCTN(ISTAP3,JC,KC-1) 
            FDIFFB = FUNCTN(ISTAP6,JC,KC+2) - FUNCTN(ISTAP6,JC,KC-2) 
     +             - FUNCTN(ISTAP2,JC,KC+2) + FUNCTN(ISTAP2,JC,KC-2) 
            FDIFFC = FUNCTN(ISTAP7,JC,KC+3) - FUNCTN(ISTAP7,JC,KC-3) 
     +             - FUNCTN(ISTAP1,JC,KC+3) + FUNCTN(ISTAP1,JC,KC-3) 
            FDIFFD = FUNCTN(ISTAP8,JC,KC+4) - FUNCTN(ISTAP8,JC,KC-4) 
     +             - FUNCTN(ISTAL,JC,KC+4)  + FUNCTN(ISTAL,JC,KC-4) 
            FDERIV(ISTAP4,JC,KC) = ACF5XZ*FDIFFA
     +                           + BCF5XZ*FDIFFB
     +                           + CCF5XZ*FDIFFC
     +                           + DCF5XZ*FDIFFD
      
          ENDDO

        ENDIF 

C       =======================================================================

C       RH END X-DIRECTION
C       ==================
        IF(NENDXR.EQ.NBOUND)THEN

C         TAKE SECOND XZ-DERIVATIVE IN X-RIGHT INNER HALO
C         EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
          DO JC = JSTAL,JSTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM3,JC,KC+1) - FUNCTN(ISTOM3,JC,KC-1) 
     +             - FUNCTN(ISTOM5,JC,KC+1) + FUNCTN(ISTOM5,JC,KC-1) 
            FDIFFB = FUNCTN(ISTOM2,JC,KC+2) - FUNCTN(ISTOM2,JC,KC-2) 
     +             - FUNCTN(ISTOM6,JC,KC+2) + FUNCTN(ISTOM6,JC,KC-2) 
            FDIFFC = FUNCTN(ISTOM1,JC,KC+3) - FUNCTN(ISTOM1,JC,KC-3) 
     +             - FUNCTN(ISTOM7,JC,KC+3) + FUNCTN(ISTOM7,JC,KC-3) 
            FDIFFD = FUNCTN(ISTOL,JC,KC+4)  - FUNCTN(ISTOL,JC,KC-4) 
     +             - FUNCTN(ISTOM8,JC,KC+4) + FUNCTN(ISTOM8,JC,KC-4) 
            FDERIV(ISTOM4,JC,KC) = ACF5XZ*FDIFFA
     +                           + BCF5XZ*FDIFFB
     +                           + CCF5XZ*FDIFFC
     +                           + DCF5XZ*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM2,JC,KC+1) - FUNCTN(ISTOM2,JC,KC-1) 
     +             - FUNCTN(ISTOM4,JC,KC+1) + FUNCTN(ISTOM4,JC,KC-1) 
            FDIFFB = FUNCTN(ISTOM1,JC,KC+2) - FUNCTN(ISTOM1,JC,KC-2) 
     +             - FUNCTN(ISTOM5,JC,KC+2) + FUNCTN(ISTOM5,JC,KC-2) 
            FDIFFC = FUNCTN(ISTOL,JC,KC+3)  - FUNCTN(ISTOL,JC,KC-3) 
     +             - FUNCTN(ISTOM6,JC,KC+3) + FUNCTN(ISTOM6,JC,KC-3) 
            FDERIV(ISTOM3,JC,KC) = ACF4XZ*FDIFFA
     +                           + BCF4XZ*FDIFFB
     +                           + CCF4XZ*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM1,JC,KC+1) - FUNCTN(ISTOM1,JC,KC-1) 
     +             - FUNCTN(ISTOM3,JC,KC+1) + FUNCTN(ISTOM3,JC,KC-1) 
            FDIFFB = FUNCTN(ISTOL,JC,KC+2)  - FUNCTN(ISTOL,JC,KC-2) 
     +             - FUNCTN(ISTOM4,JC,KC+2) + FUNCTN(ISTOM4,JC,KC-2) 
            FDERIV(ISTOM2,JC,KC) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KC+1) - FUNCTN(ISTOM1,JC,KC-1)
     +                - FUNCTN(ISTOL,JC,KC+1)  + FUNCTN(ISTOL,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KC+2) - FUNCTN(ISTOM1,JC,KC-2)
     +                - FUNCTN(ISTOL,JC,KC+2)  + FUNCTN(ISTOL,JC,KC-2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KC+1) - FUNCTN(ISTOM1,JC,KC-1)
     +                - FUNCTN(ISTOM2,JC,KC+1) + FUNCTN(ISTOM2,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KC+2) - FUNCTN(ISTOM1,JC,KC-2)
     +                - FUNCTN(ISTOM2,JC,KC+2) + FUNCTN(ISTOM2,JC,KC-2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KC+1) - FUNCTN(ISTOM1,JC,KC-1)
     +                - FUNCTN(ISTOM3,JC,KC+1) + FUNCTN(ISTOM3,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KC+2) - FUNCTN(ISTOM1,JC,KC-2)
     +                - FUNCTN(ISTOM3,JC,KC+2) + FUNCTN(ISTOM3,JC,KC-2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KC+1) - FUNCTN(ISTOM1,JC,KC-1)
     +                - FUNCTN(ISTOM4,JC,KC+1) + FUNCTN(ISTOM4,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KC+2) - FUNCTN(ISTOM1,JC,KC-2)
     +                - FUNCTN(ISTOM4,JC,KC+2) + FUNCTN(ISTOM4,JC,KC-2))
            FDERIV(ISTOM1,JC,KC) = ACF2XZ*FDIFFA
     +                           + BCF2XZ*FDIFFB
     +                           + CCF2XZ*FDIFFC
     +                           + DCF2XZ*FDIFFD

C           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KC+1)  - FUNCTN(ISTOL,JC,KC-1)
     +                - FUNCTN(ISTOM1,JC,KC+1) + FUNCTN(ISTOM1,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KC+2)  - FUNCTN(ISTOL,JC,KC-2)
     +                - FUNCTN(ISTOM1,JC,KC+2) + FUNCTN(ISTOM1,JC,KC-2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KC+1)  - FUNCTN(ISTOL,JC,KC-1)
     +                - FUNCTN(ISTOM2,JC,KC+1) + FUNCTN(ISTOM2,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KC+2)  - FUNCTN(ISTOL,JC,KC-2)
     +                - FUNCTN(ISTOM2,JC,KC+2) + FUNCTN(ISTOM2,JC,KC-2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KC+1)  - FUNCTN(ISTOL,JC,KC-1)
     +                - FUNCTN(ISTOM3,JC,KC+1) + FUNCTN(ISTOM3,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KC+2)  - FUNCTN(ISTOL,JC,KC-2)
     +                - FUNCTN(ISTOM3,JC,KC+2) + FUNCTN(ISTOM3,JC,KC-2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KC+1)  - FUNCTN(ISTOL,JC,KC-1)
     +                - FUNCTN(ISTOM4,JC,KC+1) + FUNCTN(ISTOM4,JC,KC-1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KC+2)  - FUNCTN(ISTOL,JC,KC-2)
     +                - FUNCTN(ISTOM4,JC,KC+2) + FUNCTN(ISTOM4,JC,KC-2))
            FDERIV(ISTOL,JC,KC) = ACF1XZ*FDIFFA
     +                          + BCF1XZ*FDIFFB
     +                          + CCF1XZ*FDIFFC
     +                          + DCF1XZ*FDIFFD

          ENDDO

        ENDIF

C       =======================================================================

      ENDDO

C     =========================================================================

C     LH END Z-DIRECTION
C     ==================
      IF(NENDZL.EQ.NBOUND)THEN

C       TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
C       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        DO JC = JSTAL,JSTOL
          DO IC = ISTART,IFINIS

            ICM4 = IC-4
            ICM3 = IC-3
            ICM2 = IC-2
            ICM1 = IC-1
            ICCC = IC
            ICP1 = IC+1
            ICP2 = IC+2
            ICP3 = IC+3
            ICP4 = IC+4

C           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP1) - FUNCTN(ICM1,JC,KSTAP1)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP1) - FUNCTN(ICM2,JC,KSTAP1)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP2) - FUNCTN(ICM1,JC,KSTAP2)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP2) - FUNCTN(ICM2,JC,KSTAP2)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP3) - FUNCTN(ICM1,JC,KSTAP3)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP3) - FUNCTN(ICM2,JC,KSTAP3)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP4) - FUNCTN(ICM1,JC,KSTAP4)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP4) - FUNCTN(ICM2,JC,KSTAP4)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
            FDERIV(IC,JC,KSTAL) = ACF1XZ*FDIFFA
     +                          + BCF1XZ*FDIFFB
     +                          + CCF1XZ*FDIFFC
     +                          + DCF1XZ*FDIFFD
     
C           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAL)  - FUNCTN(ICM1,JC,KSTAL)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAL)  - FUNCTN(ICM2,JC,KSTAL)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP2) - FUNCTN(ICM1,JC,KSTAP2)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP2) - FUNCTN(ICM2,JC,KSTAP2)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP3) - FUNCTN(ICM1,JC,KSTAP3)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP3) - FUNCTN(ICM2,JC,KSTAP3)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP4) - FUNCTN(ICM1,JC,KSTAP4)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP4) - FUNCTN(ICM2,JC,KSTAP4)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
            FDERIV(IC,JC,KSTAP1) = ACF2XZ*FDIFFA
     +                           + BCF2XZ*FDIFFB
     +                           + CCF2XZ*FDIFFC
     +                           + DCF2XZ*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JC,KSTAP3) - FUNCTN(ICM1,JC,KSTAP3) 
     +             - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1) 
            FDIFFB = FUNCTN(ICP2,JC,KSTAP4) - FUNCTN(ICM2,JC,KSTAP4) 
     +             - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP2) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JC,KSTAP4) - FUNCTN(ICM1,JC,KSTAP4) 
     +             - FUNCTN(ICP1,JC,KSTAP2) + FUNCTN(ICM1,JC,KSTAP2) 
            FDIFFB = FUNCTN(ICP2,JC,KSTAP5) - FUNCTN(ICM2,JC,KSTAP5) 
     +             - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1) 
            FDIFFC = FUNCTN(ICP3,JC,KSTAP6) - FUNCTN(ICM3,JC,KSTAP6) 
     +             - FUNCTN(ICP3,JC,KSTAL)  + FUNCTN(ICM3,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP3) = ACF4XZ*FDIFFA
     +                           + BCF4XZ*FDIFFB
     +                           + CCF4XZ*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JC,KSTAP5) - FUNCTN(ICM1,JC,KSTAP5) 
     +             - FUNCTN(ICP1,JC,KSTAP3) + FUNCTN(ICM1,JC,KSTAP3) 
            FDIFFB = FUNCTN(ICP2,JC,KSTAP6) - FUNCTN(ICM2,JC,KSTAP6) 
     +             - FUNCTN(ICP2,JC,KSTAP2) + FUNCTN(ICM2,JC,KSTAP2) 
            FDIFFC = FUNCTN(ICP3,JC,KSTAP7) - FUNCTN(ICM3,JC,KSTAP7) 
     +             - FUNCTN(ICP3,JC,KSTAP1) + FUNCTN(ICM3,JC,KSTAP1) 
            FDIFFD = FUNCTN(ICP4,JC,KSTAP8) - FUNCTN(ICM4,JC,KSTAP8) 
     +             - FUNCTN(ICP4,JC,KSTAL)  + FUNCTN(ICM4,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP4) = ACF5XZ*FDIFFA
     +                           + BCF5XZ*FDIFFB
     +                           + CCF5XZ*FDIFFC
     +                           + DCF5XZ*FDIFFD
      
          ENDDO
        ENDDO

C       LH IN X LH IN Z CORNER
C       ======================
        IF(NENDXL.EQ.NBOUND)THEN

          DO JC = JSTAL,JSTOL

C           LH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTAP1,JC,KSTAP1) - FUNCTN(ISTAP1,JC,KSTAL)
     +             - FUNCTN(ISTAL,JC,KSTAP1)  + FUNCTN(ISTAL,JC,KSTAL)
            FDIFFB = FUNCTN(ISTAP2,JC,KSTAP2) - FUNCTN(ISTAP2,JC,KSTAL)
     +             - FUNCTN(ISTAL,JC,KSTAP2)  + FUNCTN(ISTAL,JC,KSTAL)
            FDIFFC = FUNCTN(ISTAP3,JC,KSTAP3) - FUNCTN(ISTAP3,JC,KSTAL)
     +             - FUNCTN(ISTAL,JC,KSTAP3)  + FUNCTN(ISTAL,JC,KSTAL)
            FDIFFD = FUNCTN(ISTAP4,JC,KSTAP4) - FUNCTN(ISTAP4,JC,KSTAL)
     +             - FUNCTN(ISTAL,JC,KSTAP4)  + FUNCTN(ISTAL,JC,KSTAL)
            FDERIV(ISTAL,JC,KSTAL) = ACC1XZ*FDIFFA
     +                             + BCC1XZ*FDIFFB
     +                             + CCC1XZ*FDIFFC
     +                             + DCC1XZ*FDIFFD

C           LH+1 LH+1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTAL,JC,KSTAL)   - FUNCTN(ISTAL,JC,KSTAP1)
     +             - FUNCTN(ISTAP1,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAP1)
            FDIFFB = FUNCTN(ISTAP2,JC,KSTAP2) - FUNCTN(ISTAP2,JC,KSTAP1)
     +             - FUNCTN(ISTAP1,JC,KSTAP2) + FUNCTN(ISTAP1,JC,KSTAP1)
            FDIFFC = FUNCTN(ISTAP3,JC,KSTAP3) - FUNCTN(ISTAP3,JC,KSTAP1)
     +             - FUNCTN(ISTAP1,JC,KSTAP3) + FUNCTN(ISTAP1,JC,KSTAP1)
            FDIFFD = FUNCTN(ISTAP4,JC,KSTAP4) - FUNCTN(ISTAP4,JC,KSTAP1)
     +             - FUNCTN(ISTAP1,JC,KSTAP4) + FUNCTN(ISTAP1,JC,KSTAP1)
            FDERIV(ISTAP1,JC,KSTAP1) = ACC2XZ*FDIFFA
     +                               + BCC2XZ*FDIFFB
     +                               + CCC2XZ*FDIFFC
     +                               + DCC2XZ*FDIFFD

C           LH LH+1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTAP1,JC,KSTAL)  - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAL,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTAP1,JC,KSTAP2) - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP2)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTAP1,JC,KSTAP3) - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP3)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTAP1,JC,KSTAP4) - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP4)  + FUNCTN(ISTAL,JC,KSTAP1))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTAP2,JC,KSTAL)  - FUNCTN(ISTAP2,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAL,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTAP2) - FUNCTN(ISTAP2,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP2)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTAP2,JC,KSTAP3) - FUNCTN(ISTAP2,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP3)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTAP2,JC,KSTAP4) - FUNCTN(ISTAP2,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP4)  + FUNCTN(ISTAL,JC,KSTAP1))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTAP3,JC,KSTAL)  - FUNCTN(ISTAP3,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAL,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTAP3,JC,KSTAP2) - FUNCTN(ISTAP3,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP2)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTAP3) - FUNCTN(ISTAP3,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP3)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTAP3,JC,KSTAP4) - FUNCTN(ISTAP3,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP4)  + FUNCTN(ISTAL,JC,KSTAP1))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTAP4,JC,KSTAL)  - FUNCTN(ISTAP4,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAL,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTAP4,JC,KSTAP2) - FUNCTN(ISTAP4,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP2)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTAP4,JC,KSTAP3) - FUNCTN(ISTAP4,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP3)  + FUNCTN(ISTAL,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTAP4) - FUNCTN(ISTAP4,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAP4)  + FUNCTN(ISTAL,JC,KSTAP1))
            FDERIV(ISTAL,JC,KSTAP1) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           LH+1 LH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTAP1)  - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAP1,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTAP1) - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAP2,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTAP1) - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAP3,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTAP1) - FUNCTN(ISTAP1,JC,KSTAP1)
     +            - FUNCTN(ISTAP4,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTAP2)  - FUNCTN(ISTAP1,JC,KSTAP2)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAP1,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTAP2) - FUNCTN(ISTAP1,JC,KSTAP2)
     +            - FUNCTN(ISTAP2,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTAP2) - FUNCTN(ISTAP1,JC,KSTAP2)
     +            - FUNCTN(ISTAP3,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTAP2) - FUNCTN(ISTAP1,JC,KSTAP2)
     +            - FUNCTN(ISTAP4,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTAP3)  - FUNCTN(ISTAP1,JC,KSTAP3)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAP1,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTAP3) - FUNCTN(ISTAP1,JC,KSTAP3)
     +            - FUNCTN(ISTAP2,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTAP3) - FUNCTN(ISTAP1,JC,KSTAP3)
     +            - FUNCTN(ISTAP3,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTAP3) - FUNCTN(ISTAP1,JC,KSTAP3)
     +            - FUNCTN(ISTAP4,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTAP4)  - FUNCTN(ISTAP1,JC,KSTAP4)
     +            - FUNCTN(ISTAL,JC,KSTAL)   + FUNCTN(ISTAP1,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTAP4) - FUNCTN(ISTAP1,JC,KSTAP4)
     +            - FUNCTN(ISTAP2,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTAP4) - FUNCTN(ISTAP1,JC,KSTAP4)
     +            - FUNCTN(ISTAP3,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTAP4) - FUNCTN(ISTAP1,JC,KSTAP4)
     +            - FUNCTN(ISTAP4,JC,KSTAL)  + FUNCTN(ISTAP1,JC,KSTAL))
            FDERIV(ISTAP1,JC,KSTAL) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           LH EDGE IN Z
            DO IC = ISTAP2,ISTAP4

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP1) - FUNCTN(ICM1,JC,KSTAP1)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP1) - FUNCTN(ICM2,JC,KSTAP1)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP2) - FUNCTN(ICM1,JC,KSTAP2)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP2) - FUNCTN(ICM2,JC,KSTAP2)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP3) - FUNCTN(ICM1,JC,KSTAP3)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP3) - FUNCTN(ICM2,JC,KSTAP3)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP4) - FUNCTN(ICM1,JC,KSTAP4)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP4) - FUNCTN(ICM2,JC,KSTAP4)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDERIV(IC,JC,KSTAL) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAL)  - FUNCTN(ICM1,JC,KSTAL)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAL)  - FUNCTN(ICM2,JC,KSTAL)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP2) - FUNCTN(ICM1,JC,KSTAP2)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP2) - FUNCTN(ICM2,JC,KSTAP2)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP3) - FUNCTN(ICM1,JC,KSTAP3)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP3) - FUNCTN(ICM2,JC,KSTAP3)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP4) - FUNCTN(ICM1,JC,KSTAP4)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP4) - FUNCTN(ICM2,JC,KSTAP4)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDERIV(IC,JC,KSTAP1) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

            ENDDO

C           LH EDGE IN X
            DO KC = KSTAP2,KSTAP4

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTAP1,JC,KCP1) - FUNCTN(ISTAP1,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP1,JC,KCP2) - FUNCTN(ISTAP1,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTAP2,JC,KCP1) - FUNCTN(ISTAP2,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP2,JC,KCP2) - FUNCTN(ISTAP2,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTAP3,JC,KCP1) - FUNCTN(ISTAP3,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP3,JC,KCP2) - FUNCTN(ISTAP3,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTAP4,JC,KCP1) - FUNCTN(ISTAP4,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP4,JC,KCP2) - FUNCTN(ISTAP4,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDERIV(ISTAL,JC,KC) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTAL,JC,KCP1)  - FUNCTN(ISTAL,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAL,JC,KCP2)  - FUNCTN(ISTAL,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTAP2,JC,KCP1) - FUNCTN(ISTAP2,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP2,JC,KCP2) - FUNCTN(ISTAP2,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTAP3,JC,KCP1) - FUNCTN(ISTAP3,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP3,JC,KCP2) - FUNCTN(ISTAP3,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTAP4,JC,KCP1) - FUNCTN(ISTAP4,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP4,JC,KCP2) - FUNCTN(ISTAP4,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDERIV(ISTAP1,JC,KC) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTAP2,KSTAP4

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              IS = 0
              DO IC = ISTAP2,ISTAP4

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JC,KCP1) - FUNCTN(ICP1,JC,KCM1)
     +                 - FUNCTN(ICM1,JC,KCP1) + FUNCTN(ICM1,JC,KCM1)
                FDIFFB = FUNCTN(ICP2,JC,KCP2) - FUNCTN(ICP2,JC,KCM2)
     +                 - FUNCTN(ICM2,JC,KCP2) + FUNCTN(ICM2,JC,KCM2)
                FDERIV(IC,JC,KC) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB
                FSTORA(IS,KS) = FDIFFA
                FSTORB(IS,KS) = FDIFFB
 
              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 1
            DO KC = KSTAP3,KSTAP4

              KSM1 = KS
              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              IS = 1
              DO IC = ISTAP3,ISTAP4

                ISM1 = IS
                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JC,KCP3) - FUNCTN(ICP3,JC,KCM3)
     +                 - FUNCTN(ICM3,JC,KCP3) + FUNCTN(ICM3,JC,KCM3)
                FDERIV(IC,JC,KC) = ACF4XZ*FSTORA(IS,KS)
     +                           + BCF4XZ*FSTORB(IS,KS)
     +                           + CCF4XZ*FDIFFC
                FSTORC(ISM1,KSM1) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 3
            IS = 3
            KSM1 = 2
            ISM1 = 2
            KC = KSTAP4
            IC = ISTAP4
            KCM4 = KC-4
            KCP4 = KC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JC,KCP4) - FUNCTN(ICP4,JC,KCM4)
     +             - FUNCTN(ICM4,JC,KCP4) + FUNCTN(ICM4,JC,KCM4)
            FDERIV(IC,JC,KC) = ACF5XZ*FSTORA(IS,KS)
     +                       + BCF5XZ*FSTORB(IS,KS)
     +                       + CCF5XZ*FSTORC(ISM1,KSM1)
     +                       + DCF5XZ*FDIFFD

          ENDDO

        ENDIF

C       RH IN X LH IN Z CORNER
C       ======================
        IF(NENDXR.EQ.NBOUND)THEN

          DO JC = JSTAL,JSTOL

C           RH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTOL,JC,KSTAP1)  - FUNCTN(ISTOL,JC,KSTAL)
     +             - FUNCTN(ISTOM1,JC,KSTAP1) + FUNCTN(ISTOM1,JC,KSTAL)
            FDIFFB = FUNCTN(ISTOL,JC,KSTAP2)  - FUNCTN(ISTOL,JC,KSTAL)
     +             - FUNCTN(ISTOM2,JC,KSTAP2) + FUNCTN(ISTOM2,JC,KSTAL)
            FDIFFC = FUNCTN(ISTOL,JC,KSTAP3)  - FUNCTN(ISTOL,JC,KSTAL)
     +             - FUNCTN(ISTOM3,JC,KSTAP3) + FUNCTN(ISTOM3,JC,KSTAL)
            FDIFFD = FUNCTN(ISTOL,JC,KSTAP4)  - FUNCTN(ISTOL,JC,KSTAL)
     +             - FUNCTN(ISTOM4,JC,KSTAP4) + FUNCTN(ISTOM4,JC,KSTAL)
            FDERIV(ISTOL,JC,KSTAL) = ACC1XZ*FDIFFA
     +                             + BCC1XZ*FDIFFB
     +                             + CCC1XZ*FDIFFC
     +                             + DCC1XZ*FDIFFD

C           RH-1 LH+1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTOM1,JC,KSTAL)  - FUNCTN(ISTOM1,JC,KSTAP1)
     +             - FUNCTN(ISTOL,JC,KSTAL)   + FUNCTN(ISTOL,JC,KSTAP1)
            FDIFFB = FUNCTN(ISTOM1,JC,KSTAP2) - FUNCTN(ISTOM1,JC,KSTAP1)
     +             - FUNCTN(ISTOM2,JC,KSTAP2) + FUNCTN(ISTOM2,JC,KSTAP1)
            FDIFFC = FUNCTN(ISTOM1,JC,KSTAP3) - FUNCTN(ISTOM1,JC,KSTAP1)
     +             - FUNCTN(ISTOM3,JC,KSTAP3) + FUNCTN(ISTOM3,JC,KSTAP1)
            FDIFFD = FUNCTN(ISTOM1,JC,KSTAP4) - FUNCTN(ISTOM1,JC,KSTAP1)
     +             - FUNCTN(ISTOM4,JC,KSTAP4) + FUNCTN(ISTOM4,JC,KSTAP1)
            FDERIV(ISTOM1,JC,KSTAP1) = ACC2XZ*FDIFFA
     +                               + BCC2XZ*FDIFFB
     +                               + CCC2XZ*FDIFFC
     +                               + DCC2XZ*FDIFFD

C           RH LH+1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTAL)   - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM1,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTOL,JC,KSTAP2)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAP2) + FUNCTN(ISTOM1,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTOL,JC,KSTAP3)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAP3) + FUNCTN(ISTOM1,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTOL,JC,KSTAP4)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAP4) + FUNCTN(ISTOM1,JC,KSTAP1))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTAL)   - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM2,JC,KSTAL)  + FUNCTN(ISTOM2,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTOL,JC,KSTAP2)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM2,JC,KSTAP2) + FUNCTN(ISTOM2,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTOL,JC,KSTAP3)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM2,JC,KSTAP3) + FUNCTN(ISTOM2,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTOL,JC,KSTAP4)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM2,JC,KSTAP4) + FUNCTN(ISTOM2,JC,KSTAP1))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTAL)   - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM3,JC,KSTAL)  + FUNCTN(ISTOM3,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTOL,JC,KSTAP2)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM3,JC,KSTAP2) + FUNCTN(ISTOM3,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTOL,JC,KSTAP3)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM3,JC,KSTAP3) + FUNCTN(ISTOM3,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTOL,JC,KSTAP4)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM3,JC,KSTAP4) + FUNCTN(ISTOM3,JC,KSTAP1))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTAL)   - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM4,JC,KSTAL)  + FUNCTN(ISTOM4,JC,KSTAP1))
     +    + BCF2XZ*(FUNCTN(ISTOL,JC,KSTAP2)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM4,JC,KSTAP2) + FUNCTN(ISTOM4,JC,KSTAP1))
     +    + CCF2XZ*(FUNCTN(ISTOL,JC,KSTAP3)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM4,JC,KSTAP3) + FUNCTN(ISTOM4,JC,KSTAP1))
     +    + DCF2XZ*(FUNCTN(ISTOL,JC,KSTAP4)  - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM4,JC,KSTAP4) + FUNCTN(ISTOM4,JC,KSTAP1))
            FDERIV(ISTOL,JC,KSTAP1) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           RH-1 LH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTOM1,JC,KSTAP1) - FUNCTN(ISTOL,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOL,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP1) - FUNCTN(ISTOM2,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM2,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP1) - FUNCTN(ISTOM3,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM3,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP1) - FUNCTN(ISTOM4,JC,KSTAP1)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM4,JC,KSTAL))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTOM1,JC,KSTAP2) - FUNCTN(ISTOL,JC,KSTAP2)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOL,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP2) - FUNCTN(ISTOM2,JC,KSTAP2)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM2,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP2) - FUNCTN(ISTOM3,JC,KSTAP2)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM3,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP2) - FUNCTN(ISTOM4,JC,KSTAP2)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM4,JC,KSTAL))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTOM1,JC,KSTAP3) - FUNCTN(ISTOL,JC,KSTAP3)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOL,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP3) - FUNCTN(ISTOM2,JC,KSTAP3)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM2,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP3) - FUNCTN(ISTOM3,JC,KSTAP3)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM3,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP3) - FUNCTN(ISTOM4,JC,KSTAP3)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM4,JC,KSTAL))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTOM1,JC,KSTAP4) - FUNCTN(ISTOL,JC,KSTAP4)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOL,JC,KSTAL))
     +    + BCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP4) - FUNCTN(ISTOM2,JC,KSTAP4)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM2,JC,KSTAL))
     +    + CCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP4) - FUNCTN(ISTOM3,JC,KSTAP4)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM3,JC,KSTAL))
     +    + DCF2XZ*(FUNCTN(ISTOM1,JC,KSTAP4) - FUNCTN(ISTOM4,JC,KSTAP4)
     +            - FUNCTN(ISTOM1,JC,KSTAL)  + FUNCTN(ISTOM4,JC,KSTAL))
            FDERIV(ISTOM1,JC,KSTAL) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           LH EDGE IN Z
            DO IC = ISTOM4,ISTOM2

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP1) - FUNCTN(ICM1,JC,KSTAP1)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP1) - FUNCTN(ICM2,JC,KSTAP1)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP2) - FUNCTN(ICM1,JC,KSTAP2)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP2) - FUNCTN(ICM2,JC,KSTAP2)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP3) - FUNCTN(ICM1,JC,KSTAP3)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP3) - FUNCTN(ICM2,JC,KSTAP3)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP4) - FUNCTN(ICM1,JC,KSTAP4)
     +                - FUNCTN(ICP1,JC,KSTAL)  + FUNCTN(ICM1,JC,KSTAL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP4) - FUNCTN(ICM2,JC,KSTAP4)
     +                - FUNCTN(ICP2,JC,KSTAL)  + FUNCTN(ICM2,JC,KSTAL))
              FDERIV(IC,JC,KSTAL) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAL)  - FUNCTN(ICM1,JC,KSTAL)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAL)  - FUNCTN(ICM2,JC,KSTAL)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP2) - FUNCTN(ICM1,JC,KSTAP2)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP2) - FUNCTN(ICM2,JC,KSTAP2)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP3) - FUNCTN(ICM1,JC,KSTAP3)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP3) - FUNCTN(ICM2,JC,KSTAP3)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTAP4) - FUNCTN(ICM1,JC,KSTAP4)
     +                - FUNCTN(ICP1,JC,KSTAP1) + FUNCTN(ICM1,JC,KSTAP1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTAP4) - FUNCTN(ICM2,JC,KSTAP4)
     +                - FUNCTN(ICP2,JC,KSTAP1) + FUNCTN(ICM2,JC,KSTAP1))
              FDERIV(IC,JC,KSTAP1) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

            ENDDO

C           RH EDGE IN X
            DO KC = KSTAP2,KSTAP4

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM1,JC,KCP1) + FUNCTN(ISTOM1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM1,JC,KCP2) + FUNCTN(ISTOM1,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM2,JC,KCP1) + FUNCTN(ISTOM2,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM2,JC,KCP2) + FUNCTN(ISTOM2,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM3,JC,KCP1) + FUNCTN(ISTOM3,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM3,JC,KCP2) + FUNCTN(ISTOM3,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM4,JC,KCP1) + FUNCTN(ISTOM4,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM4,JC,KCP2) + FUNCTN(ISTOM4,JC,KCM2))
              FDERIV(ISTOL,JC,KC) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD
    
C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOL,JC,KCP1)  + FUNCTN(ISTOL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOL,JC,KCP2)  + FUNCTN(ISTOL,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOM2,JC,KCP1) + FUNCTN(ISTOM2,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOM2,JC,KCP2) + FUNCTN(ISTOM2,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOM3,JC,KCP1) + FUNCTN(ISTOM3,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOM3,JC,KCP2) + FUNCTN(ISTOM3,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOM4,JC,KCP1) + FUNCTN(ISTOM4,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOM4,JC,KCP2) + FUNCTN(ISTOM4,JC,KCM2))
              FDERIV(ISTOM1,JC,KC) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTAP2,KSTAP4

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              IS = 0
              DO IC = ISTOM4,ISTOM2

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JC,KCP1) - FUNCTN(ICP1,JC,KCM1)
     +                 - FUNCTN(ICM1,JC,KCP1) + FUNCTN(ICM1,JC,KCM1)
                FDIFFB = FUNCTN(ICP2,JC,KCP2) - FUNCTN(ICP2,JC,KCM2)
     +                 - FUNCTN(ICM2,JC,KCP2) + FUNCTN(ICM2,JC,KCM2)
                FDERIV(IC,JC,KC) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB
                FSTORA(IS,KS) = FDIFFA
                FSTORB(IS,KS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 1
            DO KC = KSTAP3,KSTAP4

              KSM1 = KS
              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              IS = 0
              DO IC = ISTOM4,ISTOM3

                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JC,KCP3) - FUNCTN(ICP3,JC,KCM3)
     +                 - FUNCTN(ICM3,JC,KCP3) + FUNCTN(ICM3,JC,KCM3)
                FDERIV(IC,JC,KC) = ACF4XZ*FSTORA(IS,KS)
     +                           + BCF4XZ*FSTORB(IS,KS)
     +                           + CCF4XZ*FDIFFC
                FSTORC(IS,KSM1) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 3
            IS = 1
            KSM1 = 2
            KC = KSTAP4
            IC = ISTOM4
            KCM4 = KC-4
            KCP4 = KC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JC,KCP4) - FUNCTN(ICP4,JC,KCM4)
     +             - FUNCTN(ICM4,JC,KCP4) + FUNCTN(ICM4,JC,KCM4)
            FDERIV(IC,JC,KC) = ACF5XZ*FSTORA(IS,KS)
     +                       + BCF5XZ*FSTORB(IS,KS)
     +                       + CCF5XZ*FSTORC(IS,KSM1)
     +                       + DCF5XZ*FDIFFD

          ENDDO

        ENDIF

      ENDIF 

C     =========================================================================

C     RH END Z-DIRECTION
C     ==================
      IF(NENDZR.EQ.NBOUND)THEN

C       TAKE SECOND XZ-DERIVATIVE IN Z-RIGHT INNER HALO
C       EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        DO JC = JSTAL,JSTOL
          DO IC = ISTART,IFINIS

            ICM4 = IC-4
            ICM3 = IC-3
            ICM2 = IC-2
            ICM1 = IC-1
            ICCC = IC
            ICP1 = IC+1
            ICP2 = IC+2
            ICP3 = IC+2
            ICP4 = IC+4

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JC,KSTOM3) - FUNCTN(ICM1,JC,KSTOM3) 
     +             - FUNCTN(ICP1,JC,KSTOM5) + FUNCTN(ICM1,JC,KSTOM5) 
            FDIFFB = FUNCTN(ICP2,JC,KSTOM2) - FUNCTN(ICM2,JC,KSTOM2) 
     +             - FUNCTN(ICP2,JC,KSTOM6) + FUNCTN(ICM2,JC,KSTOM6) 
            FDIFFC = FUNCTN(ICP3,JC,KSTOM1) - FUNCTN(ICM3,JC,KSTOM1) 
     +             - FUNCTN(ICP3,JC,KSTOM7) + FUNCTN(ICM3,JC,KSTOM7) 
            FDIFFD = FUNCTN(ICP4,JC,KSTOL)  - FUNCTN(ICM4,JC,KSTOL) 
     +             - FUNCTN(ICP4,JC,KSTOM8) + FUNCTN(ICM4,JC,KSTOM8) 
            FDERIV(IC,JC,KSTOM4) = ACF5XZ*FDIFFA
     +                           + BCF5XZ*FDIFFB
     +                           + CCF5XZ*FDIFFC
     +                           + DCF5XZ*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JC,KSTOM2) - FUNCTN(ICM1,JC,KSTOM2) 
     +             - FUNCTN(ICP1,JC,KSTOM4) + FUNCTN(ICM1,JC,KSTOM4) 
            FDIFFB = FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1) 
     +             - FUNCTN(ICP2,JC,KSTOM5) + FUNCTN(ICM2,JC,KSTOM5) 
            FDIFFC = FUNCTN(ICP3,JC,KSTOL)  - FUNCTN(ICM3,JC,KSTOL) 
     +             - FUNCTN(ICP3,JC,KSTOM6) + FUNCTN(ICM3,JC,KSTOM6) 
            FDERIV(IC,JC,KSTOM3) = ACF4XZ*FDIFFA
     +                           + BCF4XZ*FDIFFB
     +                           + CCF4XZ*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1) 
     +             - FUNCTN(ICP1,JC,KSTOM3) + FUNCTN(ICM1,JC,KSTOM3) 
            FDIFFB = FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL) 
     +             - FUNCTN(ICP2,JC,KSTOM4) + FUNCTN(ICM2,JC,KSTOM4) 
            FDERIV(IC,JC,KSTOM2) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOL)  + FUNCTN(ICM1,JC,KSTOL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOL)  + FUNCTN(ICM2,JC,KSTOL))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM2) + FUNCTN(ICM1,JC,KSTOM2))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM2) + FUNCTN(ICM2,JC,KSTOM2))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM3) + FUNCTN(ICM1,JC,KSTOM3))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM3) + FUNCTN(ICM2,JC,KSTOM3))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM4) + FUNCTN(ICM1,JC,KSTOM4))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM4) + FUNCTN(ICM2,JC,KSTOM4))
            FDERIV(IC,JC,KSTOM1) = ACF2XZ*FDIFFA
     +                           + BCF2XZ*FDIFFB
     +                           + CCF2XZ*FDIFFC
     +                           + DCF2XZ*FDIFFD

C           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM1) + FUNCTN(ICM1,JC,KSTOM1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM1) + FUNCTN(ICM2,JC,KSTOM1))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM2) + FUNCTN(ICM1,JC,KSTOM2))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM2) + FUNCTN(ICM2,JC,KSTOM2))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM3) + FUNCTN(ICM1,JC,KSTOM3))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM3) + FUNCTN(ICM2,JC,KSTOM3))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM4) + FUNCTN(ICM1,JC,KSTOM4))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM4) + FUNCTN(ICM2,JC,KSTOM4))
            FDERIV(IC,JC,KSTOL) = ACF1XZ*FDIFFA
     +                          + BCF1XZ*FDIFFB
     +                          + CCF1XZ*FDIFFC
     +                          + DCF1XZ*FDIFFD

          ENDDO
        ENDDO

C       LH IN X RH IN Z CORNER
C       ======================
        IF(NENDXL.EQ.NBOUND)THEN

          DO JC = JSTAL,JSTOL

C           LH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTAP1,JC,KSTOL) - FUNCTN(ISTAP1,JC,KSTOM1)
     +             - FUNCTN(ISTAL,JC,KSTOL)  + FUNCTN(ISTAL,JC,KSTOM1) 
            FDIFFB = FUNCTN(ISTAP2,JC,KSTOL) - FUNCTN(ISTAP2,JC,KSTOM2)
     +             - FUNCTN(ISTAL,JC,KSTOL)  + FUNCTN(ISTAL,JC,KSTOM2) 
            FDIFFC = FUNCTN(ISTAP3,JC,KSTOL) - FUNCTN(ISTAP3,JC,KSTOM3)
     +             - FUNCTN(ISTAL,JC,KSTOL)  + FUNCTN(ISTAL,JC,KSTOM3) 
            FDIFFD = FUNCTN(ISTAP4,JC,KSTOL) - FUNCTN(ISTAP4,JC,KSTOM4)
     +             - FUNCTN(ISTAL,JC,KSTOL)  + FUNCTN(ISTAL,JC,KSTOM4) 
            FDERIV(ISTAL,JC,KSTOL) = ACC1XZ*FDIFFA
     +                             + BCC1XZ*FDIFFB
     +                             + CCC1XZ*FDIFFC
     +                             + DCC1XZ*FDIFFD

C           LH+1 RH-1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTAL,JC,KSTOM1)  - FUNCTN(ISTAL,JC,KSTOL)
     +             - FUNCTN(ISTAP1,JC,KSTOM1) + FUNCTN(ISTAP1,JC,KSTOL)
            FDIFFB = FUNCTN(ISTAP2,JC,KSTOM1) - FUNCTN(ISTAP2,JC,KSTOM2)
     +             - FUNCTN(ISTAP1,JC,KSTOM1) + FUNCTN(ISTAP1,JC,KSTOM2)
            FDIFFC = FUNCTN(ISTAP3,JC,KSTOM1) - FUNCTN(ISTAP3,JC,KSTOM3)
     +             - FUNCTN(ISTAP1,JC,KSTOM1) + FUNCTN(ISTAP1,JC,KSTOM3)
            FDIFFD = FUNCTN(ISTAP4,JC,KSTOM1) - FUNCTN(ISTAP4,JC,KSTOM4)
     +             - FUNCTN(ISTAP1,JC,KSTOM1) + FUNCTN(ISTAP1,JC,KSTOM4)
            FDERIV(ISTAP1,JC,KSTOM1) = ACC2XZ*FDIFFA
     +                               + BCC2XZ*FDIFFB
     +                               + CCC2XZ*FDIFFC
     +                               + DCC2XZ*FDIFFD

C           LH RH-1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTAP1,JC,KSTOM1) - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTAP1,JC,KSTOM1) - FUNCTN(ISTAP1,JC,KSTOM2)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM2))
     +    + CCF2XZ*(FUNCTN(ISTAP1,JC,KSTOM1) - FUNCTN(ISTAP1,JC,KSTOM3)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM3))
     +    + DCF2XZ*(FUNCTN(ISTAP1,JC,KSTOM1) - FUNCTN(ISTAP1,JC,KSTOM4)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM4))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTAP2,JC,KSTOM1) - FUNCTN(ISTAP2,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTOM1) - FUNCTN(ISTAP2,JC,KSTOM2)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM2))
     +    + CCF2XZ*(FUNCTN(ISTAP2,JC,KSTOM1) - FUNCTN(ISTAP2,JC,KSTOM3)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM3))
     +    + DCF2XZ*(FUNCTN(ISTAP2,JC,KSTOM1) - FUNCTN(ISTAP2,JC,KSTOM4)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM4))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTAP3,JC,KSTOM1) - FUNCTN(ISTAP3,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTAP3,JC,KSTOM1) - FUNCTN(ISTAP3,JC,KSTOM2)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM2))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTOM1) - FUNCTN(ISTAP3,JC,KSTOM3)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM3))
     +    + DCF2XZ*(FUNCTN(ISTAP3,JC,KSTOM1) - FUNCTN(ISTAP3,JC,KSTOM4)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM4))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTAP4,JC,KSTOM1) - FUNCTN(ISTAP4,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTAP4,JC,KSTOM1) - FUNCTN(ISTAP4,JC,KSTOM2)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM2))
     +    + CCF2XZ*(FUNCTN(ISTAP4,JC,KSTOM1) - FUNCTN(ISTAP4,JC,KSTOM3)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM3))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTOM1) - FUNCTN(ISTAP4,JC,KSTOM4)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAL,JC,KSTOM4))
            FDERIV(ISTAL,JC,KSTOM1) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           LH+1 RH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTOL)   - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM1)  + FUNCTN(ISTAP1,JC,KSTOM1))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP2,JC,KSTOM1) + FUNCTN(ISTAP1,JC,KSTOM1))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP3,JC,KSTOM1) + FUNCTN(ISTAP1,JC,KSTOM1))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP4,JC,KSTOM1) + FUNCTN(ISTAP1,JC,KSTOM1))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTOL)   - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM2)  + FUNCTN(ISTAP1,JC,KSTOM2))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP2,JC,KSTOM2) + FUNCTN(ISTAP1,JC,KSTOM2))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP3,JC,KSTOM2) + FUNCTN(ISTAP1,JC,KSTOM2))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP4,JC,KSTOM2) + FUNCTN(ISTAP1,JC,KSTOM2))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTOL)   - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM3)  + FUNCTN(ISTAP1,JC,KSTOM3))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP2,JC,KSTOM3) + FUNCTN(ISTAP1,JC,KSTOM3))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP3,JC,KSTOM3) + FUNCTN(ISTAP1,JC,KSTOM3))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP4,JC,KSTOM3) + FUNCTN(ISTAP1,JC,KSTOM3))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTAL,JC,KSTOL)   - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAL,JC,KSTOM4)  + FUNCTN(ISTAP1,JC,KSTOM4))
     +    + BCF2XZ*(FUNCTN(ISTAP2,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP2,JC,KSTOM4) + FUNCTN(ISTAP1,JC,KSTOM4))
     +    + CCF2XZ*(FUNCTN(ISTAP3,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP3,JC,KSTOM4) + FUNCTN(ISTAP1,JC,KSTOM4))
     +    + DCF2XZ*(FUNCTN(ISTAP4,JC,KSTOL)  - FUNCTN(ISTAP1,JC,KSTOL)
     +            - FUNCTN(ISTAP4,JC,KSTOM4) + FUNCTN(ISTAP1,JC,KSTOM4))
            FDERIV(ISTAP1,JC,KSTOL) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           RH EDGE IN Z
            DO IC = ISTAP2,ISTAP4

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOL)  + FUNCTN(ICM1,JC,KSTOL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOL)  + FUNCTN(ICM2,JC,KSTOL))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM2) + FUNCTN(ICM1,JC,KSTOM2))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM2) + FUNCTN(ICM2,JC,KSTOM2))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM3) + FUNCTN(ICM1,JC,KSTOM3))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM3) + FUNCTN(ICM2,JC,KSTOM3))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM4) + FUNCTN(ICM1,JC,KSTOM4))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM4) + FUNCTN(ICM2,JC,KSTOM4))
              FDERIV(IC,JC,KSTOM1) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM1) + FUNCTN(ICM1,JC,KSTOM1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM1) + FUNCTN(ICM2,JC,KSTOM1))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM2) + FUNCTN(ICM1,JC,KSTOM2))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM2) + FUNCTN(ICM2,JC,KSTOM2))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM3) + FUNCTN(ICM1,JC,KSTOM3))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM3) + FUNCTN(ICM2,JC,KSTOM3))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM4) + FUNCTN(ICM1,JC,KSTOM4))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM4) + FUNCTN(ICM2,JC,KSTOM4))
              FDERIV(IC,JC,KSTOL) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD

            ENDDO

C           LH EDGE IN X
            DO KC = KSTOM4,KSTOM2

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTAP1,JC,KCP1) - FUNCTN(ISTAP1,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP1,JC,KCP2) - FUNCTN(ISTAP1,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTAP2,JC,KCP1) - FUNCTN(ISTAP2,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP2,JC,KCP2) - FUNCTN(ISTAP2,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTAP3,JC,KCP1) - FUNCTN(ISTAP3,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP3,JC,KCP2) - FUNCTN(ISTAP3,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTAP4,JC,KCP1) - FUNCTN(ISTAP4,JC,KCM1)
     +                - FUNCTN(ISTAL,JC,KCP1)  + FUNCTN(ISTAL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP4,JC,KCP2) - FUNCTN(ISTAP4,JC,KCM2)
     +                - FUNCTN(ISTAL,JC,KCP2)  + FUNCTN(ISTAL,JC,KCM2))
              FDERIV(ISTAL,JC,KC) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTAL,JC,KCP1)  - FUNCTN(ISTAL,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAL,JC,KCP2)  - FUNCTN(ISTAL,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTAP2,JC,KCP1) - FUNCTN(ISTAP2,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP2,JC,KCP2) - FUNCTN(ISTAP2,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTAP3,JC,KCP1) - FUNCTN(ISTAP3,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP3,JC,KCP2) - FUNCTN(ISTAP3,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTAP4,JC,KCP1) - FUNCTN(ISTAP4,JC,KCM1)
     +                - FUNCTN(ISTAP1,JC,KCP1) + FUNCTN(ISTAP1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTAP4,JC,KCP2) - FUNCTN(ISTAP4,JC,KCM2)
     +                - FUNCTN(ISTAP1,JC,KCP2) + FUNCTN(ISTAP1,JC,KCM2))
              FDERIV(ISTAP1,JC,KC) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM2

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              IS = 0
              DO IC = ISTAP2,ISTAP4

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JC,KCP1) - FUNCTN(ICP1,JC,KCM1)
     +                 - FUNCTN(ICM1,JC,KCP1) + FUNCTN(ICM1,JC,KCM1)
                FDIFFB = FUNCTN(ICP2,JC,KCP2) - FUNCTN(ICP2,JC,KCM2)
     +                 - FUNCTN(ICM2,JC,KCP2) + FUNCTN(ICM2,JC,KCM2)
                FDERIV(IC,JC,KC) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB
                FSTORA(IS,KS) = FDIFFA
                FSTORB(IS,KS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM3

              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              IS = 1
              DO IC = ISTAP3,ISTAP4

                ISM1 = IS
                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JC,KCP3) - FUNCTN(ICP3,JC,KCM3)
     +                 - FUNCTN(ICM3,JC,KCP3) + FUNCTN(ICM3,JC,KCM3)
                FDERIV(IC,JC,KC) = ACF4XZ*FSTORA(IS,KS)
     +                           + BCF4XZ*FSTORB(IS,KS)
     +                           + CCF4XZ*FDIFFC
                FSTORC(ISM1,KS) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 1
            IS = 3
            ISM1 = 2
            KC = KSTOM4
            IC = ISTAP4
            KCM4 = KC-4
            KCP4 = KC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JC,KCP4) - FUNCTN(ICP4,JC,KCM4)
     +             - FUNCTN(ICM4,JC,KCP4) + FUNCTN(ICM4,JC,KCM4)
            FDERIV(IC,JC,KC) = ACF5XZ*FSTORA(IS,KS)
     +                       + BCF5XZ*FSTORB(IS,KS)
     +                       + CCF5XZ*FSTORC(ISM1,KS)
     +                       + DCF5XZ*FDIFFD

          ENDDO

        ENDIF

C       RH IN X RH IN Z CORNER
C       ======================
        IF(NENDXR.EQ.NBOUND)THEN

          DO JC = JSTAL,JSTOL

C           RH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTOM1,JC,KSTOM1) - FUNCTN(ISTOM1,JC,KSTOL)
     +             - FUNCTN(ISTOL,JC,KSTOM1)  + FUNCTN(ISTOL,JC,KSTOL)
            FDIFFB = FUNCTN(ISTOM2,JC,KSTOM2) - FUNCTN(ISTOM2,JC,KSTOL)
     +             - FUNCTN(ISTOL,JC,KSTOM2)  + FUNCTN(ISTOL,JC,KSTOL)
            FDIFFC = FUNCTN(ISTOM3,JC,KSTOM3) - FUNCTN(ISTOM3,JC,KSTOL)
     +             - FUNCTN(ISTOL,JC,KSTOM3)  + FUNCTN(ISTOL,JC,KSTOL)
            FDIFFD = FUNCTN(ISTOM4,JC,KSTOM4) - FUNCTN(ISTOM4,JC,KSTOL)
     +             - FUNCTN(ISTOL,JC,KSTOM4)  + FUNCTN(ISTOL,JC,KSTOL)
            FDERIV(ISTOL,JC,KSTOL) = ACC1XZ*FDIFFA
     +                             + BCC1XZ*FDIFFB
     +                             + CCC1XZ*FDIFFC
     +                             + DCC1XZ*FDIFFD

C           RH-1 RH-1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTOL,JC,KSTOL)   - FUNCTN(ISTOL,JC,KSTOM1)
     +             - FUNCTN(ISTOM1,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOM1)
            FDIFFB = FUNCTN(ISTOM2,JC,KSTOM2) - FUNCTN(ISTOM2,JC,KSTOM1)
     +             - FUNCTN(ISTOM1,JC,KSTOM2) + FUNCTN(ISTOM1,JC,KSTOM1)
            FDIFFC = FUNCTN(ISTOM3,JC,KSTOM3) - FUNCTN(ISTOM3,JC,KSTOM1)
     +             - FUNCTN(ISTOM1,JC,KSTOM3) + FUNCTN(ISTOM1,JC,KSTOM1)
            FDIFFD = FUNCTN(ISTOM4,JC,KSTOM4) - FUNCTN(ISTOM4,JC,KSTOM1)
     +             - FUNCTN(ISTOM1,JC,KSTOM4) + FUNCTN(ISTOM1,JC,KSTOM1)
            FDERIV(ISTOM1,JC,KSTOM1) = ACC2XZ*FDIFFA
     +                               + BCC2XZ*FDIFFB
     +                               + CCC2XZ*FDIFFC
     +                               + DCC2XZ*FDIFFD

C           RH RH-1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTOM1,JC,KSTOL)  - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOL,JC,KSTOM1))
     +    + BCF2XZ*(FUNCTN(ISTOM1,JC,KSTOM2) - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM2)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + CCF2XZ*(FUNCTN(ISTOM1,JC,KSTOM3) - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM3)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + DCF2XZ*(FUNCTN(ISTOM1,JC,KSTOM4) - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM4)  + FUNCTN(ISTOL,JC,KSTOM1))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTOM2,JC,KSTOL)  - FUNCTN(ISTOM2,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOL,JC,KSTOM1))
     +    + BCF2XZ*(FUNCTN(ISTOM2,JC,KSTOM2) - FUNCTN(ISTOM2,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM2)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + CCF2XZ*(FUNCTN(ISTOM2,JC,KSTOM3) - FUNCTN(ISTOM2,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM3)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + DCF2XZ*(FUNCTN(ISTOM2,JC,KSTOM4) - FUNCTN(ISTOM2,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM4)  + FUNCTN(ISTOL,JC,KSTOM1))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTOM3,JC,KSTOL)  - FUNCTN(ISTOM3,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOL,JC,KSTOM1))
     +    + BCF2XZ*(FUNCTN(ISTOM3,JC,KSTOM2) - FUNCTN(ISTOM3,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM2)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + CCF2XZ*(FUNCTN(ISTOM3,JC,KSTOM3) - FUNCTN(ISTOM3,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM3)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + DCF2XZ*(FUNCTN(ISTOM3,JC,KSTOM4) - FUNCTN(ISTOM3,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM4)  + FUNCTN(ISTOL,JC,KSTOM1))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTOM4,JC,KSTOL)  - FUNCTN(ISTOM4,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOL,JC,KSTOM1))
     +    + BCF2XZ*(FUNCTN(ISTOM4,JC,KSTOM2) - FUNCTN(ISTOM4,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM2)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + CCF2XZ*(FUNCTN(ISTOM4,JC,KSTOM3) - FUNCTN(ISTOM4,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM3)  + FUNCTN(ISTOL,JC,KSTOM1))
     +    + DCF2XZ*(FUNCTN(ISTOM4,JC,KSTOM4) - FUNCTN(ISTOM4,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOM4)  + FUNCTN(ISTOL,JC,KSTOM1))
            FDERIV(ISTOL,JC,KSTOM1) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           RH+1 RH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTOM1)  - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOM1,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTOM2,JC,KSTOM1) - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOM2,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + CCF2XZ*(FUNCTN(ISTOM3,JC,KSTOM1) - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOM3,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + DCF2XZ*(FUNCTN(ISTOM4,JC,KSTOM1) - FUNCTN(ISTOM1,JC,KSTOM1)
     +            - FUNCTN(ISTOM4,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
            FDIFFB =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTOM2)  - FUNCTN(ISTOM1,JC,KSTOM2)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOM1,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTOM2,JC,KSTOM2) - FUNCTN(ISTOM1,JC,KSTOM2)
     +            - FUNCTN(ISTOM2,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + CCF2XZ*(FUNCTN(ISTOM3,JC,KSTOM2) - FUNCTN(ISTOM1,JC,KSTOM2)
     +            - FUNCTN(ISTOM3,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + DCF2XZ*(FUNCTN(ISTOM4,JC,KSTOM2) - FUNCTN(ISTOM1,JC,KSTOM2)
     +            - FUNCTN(ISTOM4,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
            FDIFFC =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTOM3)  - FUNCTN(ISTOM1,JC,KSTOM3)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOM1,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTOM2,JC,KSTOM3) - FUNCTN(ISTOM1,JC,KSTOM3)
     +            - FUNCTN(ISTOM2,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + CCF2XZ*(FUNCTN(ISTOM3,JC,KSTOM3) - FUNCTN(ISTOM1,JC,KSTOM3)
     +            - FUNCTN(ISTOM3,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + DCF2XZ*(FUNCTN(ISTOM4,JC,KSTOM3) - FUNCTN(ISTOM1,JC,KSTOM3)
     +            - FUNCTN(ISTOM4,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
            FDIFFD =
     +      ACF2XZ*(FUNCTN(ISTOL,JC,KSTOM4)  - FUNCTN(ISTOM1,JC,KSTOM4)
     +            - FUNCTN(ISTOL,JC,KSTOL)   + FUNCTN(ISTOM1,JC,KSTOL))
     +    + BCF2XZ*(FUNCTN(ISTOM2,JC,KSTOM4) - FUNCTN(ISTOM1,JC,KSTOM4)
     +            - FUNCTN(ISTOM2,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + CCF2XZ*(FUNCTN(ISTOM3,JC,KSTOM4) - FUNCTN(ISTOM1,JC,KSTOM4)
     +            - FUNCTN(ISTOM3,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
     +    + DCF2XZ*(FUNCTN(ISTOM4,JC,KSTOM4) - FUNCTN(ISTOM1,JC,KSTOM4)
     +            - FUNCTN(ISTOM4,JC,KSTOL)  + FUNCTN(ISTOM1,JC,KSTOL))
            FDERIV(ISTOM1,JC,KSTOL) = ACF1XZ*FDIFFA
     +                              + BCF1XZ*FDIFFB
     +                              + CCF1XZ*FDIFFC
     +                              + DCF1XZ*FDIFFD

C           RH EDGE IN Z
            DO IC = ISTOM4,ISTOM2

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOL)  + FUNCTN(ICM1,JC,KSTOL))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOL)  + FUNCTN(ICM2,JC,KSTOL))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM2) + FUNCTN(ICM1,JC,KSTOM2))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM2) + FUNCTN(ICM2,JC,KSTOM2))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM3) + FUNCTN(ICM1,JC,KSTOM3))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM3) + FUNCTN(ICM2,JC,KSTOM3))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOM1) - FUNCTN(ICM1,JC,KSTOM1)
     +                - FUNCTN(ICP1,JC,KSTOM4) + FUNCTN(ICM1,JC,KSTOM4))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOM1) - FUNCTN(ICM2,JC,KSTOM1)
     +                - FUNCTN(ICP2,JC,KSTOM4) + FUNCTN(ICM2,JC,KSTOM4))
              FDERIV(IC,JC,KSTOM1) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM1) + FUNCTN(ICM1,JC,KSTOM1))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM1) + FUNCTN(ICM2,JC,KSTOM1))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM2) + FUNCTN(ICM1,JC,KSTOM2))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM2) + FUNCTN(ICM2,JC,KSTOM2))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM3) + FUNCTN(ICM1,JC,KSTOM3))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM3) + FUNCTN(ICM2,JC,KSTOM3))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JC,KSTOL)  - FUNCTN(ICM1,JC,KSTOL)
     +                - FUNCTN(ICP1,JC,KSTOM4) + FUNCTN(ICM1,JC,KSTOM4))
     +        + BCOFX1*(FUNCTN(ICP2,JC,KSTOL)  - FUNCTN(ICM2,JC,KSTOL)
     +                - FUNCTN(ICP2,JC,KSTOM4) + FUNCTN(ICM2,JC,KSTOM4))
              FDERIV(IC,JC,KSTOL) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD

            ENDDO

C           RH EDGE IN X
            DO KC = KSTOM4,KSTOM2

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM1,JC,KCP1) + FUNCTN(ISTOM1,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM1,JC,KCP2) + FUNCTN(ISTOM1,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM2,JC,KCP1) + FUNCTN(ISTOM2,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM2,JC,KCP2) + FUNCTN(ISTOM2,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM3,JC,KCP1) + FUNCTN(ISTOM3,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM3,JC,KCP2) + FUNCTN(ISTOM3,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTOL,JC,KCP1)  - FUNCTN(ISTOL,JC,KCM1)
     +                - FUNCTN(ISTOM4,JC,KCP1) + FUNCTN(ISTOM4,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOL,JC,KCP2)  - FUNCTN(ISTOL,JC,KCM2)
     +                - FUNCTN(ISTOM4,JC,KCP2) + FUNCTN(ISTOM4,JC,KCM2))
              FDERIV(ISTOL,JC,KC) = ACF1XZ*FDIFFA
     +                            + BCF1XZ*FDIFFB
     +                            + CCF1XZ*FDIFFC
     +                            + DCF1XZ*FDIFFD
     
C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOL,JC,KCP1)  + FUNCTN(ISTOL,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOL,JC,KCP2)  + FUNCTN(ISTOL,JC,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOM2,JC,KCP1) + FUNCTN(ISTOM2,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOM2,JC,KCP2) + FUNCTN(ISTOM2,JC,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOM3,JC,KCP1) + FUNCTN(ISTOM3,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOM3,JC,KCP2) + FUNCTN(ISTOM3,JC,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(ISTOM1,JC,KCP1) - FUNCTN(ISTOM1,JC,KCM1)
     +                - FUNCTN(ISTOM4,JC,KCP1) + FUNCTN(ISTOM4,JC,KCM1))
     +        + BCOFZ1*(FUNCTN(ISTOM1,JC,KCP2) - FUNCTN(ISTOM1,JC,KCM2)
     +                - FUNCTN(ISTOM4,JC,KCP2) + FUNCTN(ISTOM4,JC,KCM2))
              FDERIV(ISTOM1,JC,KC) = ACF2XZ*FDIFFA
     +                             + BCF2XZ*FDIFFB
     +                             + CCF2XZ*FDIFFC
     +                             + DCF2XZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM2

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              IS = 0
              DO IC = ISTOM4,ISTOM2

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JC,KCP1) - FUNCTN(ICP1,JC,KCM1)
     +                 - FUNCTN(ICM1,JC,KCP1) + FUNCTN(ICM1,JC,KCM1)
                FDIFFB = FUNCTN(ICP2,JC,KCP2) - FUNCTN(ICP2,JC,KCM2)
     +                 - FUNCTN(ICM2,JC,KCP2) + FUNCTN(ICM2,JC,KCM2)
                FDERIV(IC,JC,KC) = ACF3XZ*FDIFFA
     +                           + BCF3XZ*FDIFFB
                FSTORA(IS,KS) = FDIFFA
                FSTORB(IS,KS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM3

              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              IS = 0
              DO IC = ISTOM4,ISTOM3

                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JC,KCP3) - FUNCTN(ICP3,JC,KCM3)
     +                 - FUNCTN(ICM3,JC,KCP3) + FUNCTN(ICM3,JC,KCM3)
                FDERIV(IC,JC,KC) = ACF4XZ*FSTORA(IS,KS)
     +                           + BCF4XZ*FSTORB(IS,KS)
     +                           + CCF4XZ*FDIFFC
                FSTORC(IS,KS) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 1
            IS = 1
            KC = KSTOM4
            IC = ISTOM4
            KCM4 = KC-4
            KCP4 = KC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JC,KCP4) - FUNCTN(ICP4,JC,KCM4)
     +             - FUNCTN(ICM4,JC,KCP4) + FUNCTN(ICM4,JC,KCM4)
            FDERIV(IC,JC,KC) = ACF5XZ*FSTORA(IS,KS)
     +                       + BCF5XZ*FSTORB(IS,KS)
     +                       + CCF5XZ*FSTORC(IS,KS)
     +                       + DCF5XZ*FDIFFD

          ENDDO

        ENDIF

      ENDIF

C     =========================================================================

C     SCALING
C     =======
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDELX*OVDELZ

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
