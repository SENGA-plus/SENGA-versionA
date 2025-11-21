      SUBROUTINE D2FDYZ(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDYZ
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
C     EVALUATES SECOND YZ-DERIVATIVE OF SPECIFIED FUNCTION
C     EXPLICIT 10TH ORDER FINITE DIFFERENCE METHOD
C     EXPLICIT 8TH,6TH,4TH,4TH,4TH COMPATIBLE ORDER END CONDITIONS
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
      INTEGER JS,KS,JSM1,KSM1
      INTEGER JSTART,JFINIS,KSTART,KFINIS
      INTEGER JCM5,JCM4,JCM3,JCM2,JCM1,JCCC,JCP1,JCP2,JCP3,JCP4,JCP5
      INTEGER KCM5,KCM4,KCM3,KCM2,KCM1,KCCC,KCP1,KCP2,KCP3,KCP4,KCP5


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      JSTART = JSTAL
      JFINIS = JSTOL
      KSTART = KSTAL
      KFINIS = KSTOL
      IF(NENDYL.EQ.NBOUND)JSTART = JSTAP5
      IF(NENDYR.EQ.NBOUND)JFINIS = JSTOM5
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
             
            FDIFFA = FUNCTN(IC,JCP1,KCP1) - FUNCTN(IC,JCP1,KCM1) 
     +             - FUNCTN(IC,JCM1,KCP1) + FUNCTN(IC,JCM1,KCM1) 
            FDIFFB = FUNCTN(IC,JCP2,KCP2) - FUNCTN(IC,JCP2,KCM2) 
     +             - FUNCTN(IC,JCM2,KCP2) + FUNCTN(IC,JCM2,KCM2) 
            FDIFFC = FUNCTN(IC,JCP3,KCP3) - FUNCTN(IC,JCP3,KCM3) 
     +             - FUNCTN(IC,JCM3,KCP3) + FUNCTN(IC,JCM3,KCM3) 
            FDIFFD = FUNCTN(IC,JCP4,KCP4) - FUNCTN(IC,JCP4,KCM4) 
     +             - FUNCTN(IC,JCM4,KCP4) + FUNCTN(IC,JCM4,KCM4) 
            FDIFFE = FUNCTN(IC,JCP5,KCP5) - FUNCTN(IC,JCP5,KCM5) 
     +             - FUNCTN(IC,JCM5,KCP5) + FUNCTN(IC,JCM5,KCM5) 

            FDERIV(IC,JC,KC) = ACOFYZ*FDIFFA
     +                       + BCOFYZ*FDIFFB
     +                       + CCOFYZ*FDIFFC
     +                       + DCOFYZ*FDIFFD
     +                       + ECOFYZ*FDIFFE

          ENDDO

        ENDDO

      ENDDO

C     =========================================================================

C     LH END Y-DIRECTION
C     ==================
      IF(NENDYL.EQ.NBOUND)THEN

C       TAKE SECOND YZ-DERIVATIVE IN Y-LEFT INNER HALO
C       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        KCM3 = KSTART-4
        KCM2 = KSTART-3
        KCM1 = KSTART-2
        KCCC = KSTART-1
        KCP1 = KSTART
        KCP2 = KSTART+1
        KCP3 = KSTART+2
        KCP4 = KSTART+3

        DO KC = KSTART,KFINIS

          KCM4 = KCM3
          KCM3 = KCM2
          KCM2 = KCM1
          KCM1 = KCCC
          KCCC = KCP1
          KCP1 = KCP2
          KCP2 = KCP3
          KCP3 = KCP4
          KCP4 = KC+4

          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTAP1,KCP1) - FUNCTN(IC,JSTAP1,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP1,KCP2) - FUNCTN(IC,JSTAP1,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTAP2,KCP1) - FUNCTN(IC,JSTAP2,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP2,KCP2) - FUNCTN(IC,JSTAP2,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTAP3,KCP1) - FUNCTN(IC,JSTAP3,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP3,KCP2) - FUNCTN(IC,JSTAP3,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTAP4,KCP1) - FUNCTN(IC,JSTAP4,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP4,KCP2) - FUNCTN(IC,JSTAP4,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
            FDERIV(IC,JSTAL,KC) = ACF1YZ*FDIFFA
     +                          + BCF1YZ*FDIFFB
     +                          + CCF1YZ*FDIFFC
     +                          + DCF1YZ*FDIFFD
     
C           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTAL,KCP1)  - FUNCTN(IC,JSTAL,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAL,KCP2)  - FUNCTN(IC,JSTAL,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTAP2,KCP1) - FUNCTN(IC,JSTAP2,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP2,KCP2) - FUNCTN(IC,JSTAP2,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTAP3,KCP1) - FUNCTN(IC,JSTAP3,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP3,KCP2) - FUNCTN(IC,JSTAP3,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTAP4,KCP1) - FUNCTN(IC,JSTAP4,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP4,KCP2) - FUNCTN(IC,JSTAP4,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
            FDERIV(IC,JSTAP1,KC) = ACF2YZ*FDIFFA
     +                           + BCF2YZ*FDIFFB
     +                           + CCF2YZ*FDIFFC
     +                           + DCF2YZ*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP3,KCP1) - FUNCTN(IC,JSTAP3,KCM1) 
     +             - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1) 
            FDIFFB = FUNCTN(IC,JSTAP4,KCP2) - FUNCTN(IC,JSTAP4,KCM2) 
     +             - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2) 
            FDERIV(IC,JSTAP2,KC) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP4,KCP1) - FUNCTN(IC,JSTAP4,KCM1) 
     +             - FUNCTN(IC,JSTAP2,KCP1) + FUNCTN(IC,JSTAP2,KCM1) 
            FDIFFB = FUNCTN(IC,JSTAP5,KCP2) - FUNCTN(IC,JSTAP5,KCM2) 
     +             - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2) 
            FDIFFC = FUNCTN(IC,JSTAP6,KCP3) - FUNCTN(IC,JSTAP6,KCM3) 
     +             - FUNCTN(IC,JSTAL,KCP3)  + FUNCTN(IC,JSTAL,KCM3) 
            FDERIV(IC,JSTAP3,KC) = ACF4YZ*FDIFFA
     +                           + BCF4YZ*FDIFFB
     +                           + CCF4YZ*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP5,KCP1) - FUNCTN(IC,JSTAP5,KCM1) 
     +             - FUNCTN(IC,JSTAP3,KCP1) + FUNCTN(IC,JSTAP3,KCM1) 
            FDIFFB = FUNCTN(IC,JSTAP6,KCP2) - FUNCTN(IC,JSTAP6,KCM2) 
     +             - FUNCTN(IC,JSTAP2,KCP2) + FUNCTN(IC,JSTAP2,KCM2) 
            FDIFFC = FUNCTN(IC,JSTAP7,KCP3) - FUNCTN(IC,JSTAP7,KCM3) 
     +             - FUNCTN(IC,JSTAP1,KCP3) + FUNCTN(IC,JSTAP1,KCM3) 
            FDIFFD = FUNCTN(IC,JSTAP8,KCP4) - FUNCTN(IC,JSTAP8,KCM4) 
     +             - FUNCTN(IC,JSTAL,KCP4)  + FUNCTN(IC,JSTAL,KCM4) 
            FDERIV(IC,JSTAP4,KC) = ACF5YZ*FDIFFA
     +                           + BCF5YZ*FDIFFB
     +                           + CCF5YZ*FDIFFC
     +                           + DCF5YZ*FDIFFD
      
          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END Y-DIRECTION
C     ==================
      IF(NENDYR.EQ.NBOUND)THEN

C       TAKE SECOND YZ-DERIVATIVE IN Y-RIGHT INNER HALO
C       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        KCM3 = KSTART-4
        KCM2 = KSTART-3
        KCM1 = KSTART-2
        KCCC = KSTART-1
        KCP1 = KSTART
        KCP2 = KSTART+1
        KCP3 = KSTART+2
        KCP4 = KSTART+3

        DO KC = KSTART,KFINIS

          KCM4 = KCM3
          KCM3 = KCM2
          KCM2 = KCM1
          KCM1 = KCCC
          KCCC = KCP1
          KCP1 = KCP2
          KCP2 = KCP3
          KCP3 = KCP4
          KCP4 = KC+4

          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM3,KCP1) - FUNCTN(IC,JSTOM3,KCM1) 
     +             - FUNCTN(IC,JSTOM5,KCP1) + FUNCTN(IC,JSTOM5,KCM1) 
            FDIFFB = FUNCTN(IC,JSTOM2,KCP2) - FUNCTN(IC,JSTOM2,KCM2) 
     +             - FUNCTN(IC,JSTOM6,KCP2) + FUNCTN(IC,JSTOM6,KCM2) 
            FDIFFC = FUNCTN(IC,JSTOM1,KCP3) - FUNCTN(IC,JSTOM1,KCM3) 
     +             - FUNCTN(IC,JSTOM7,KCP3) + FUNCTN(IC,JSTOM7,KCM3) 
            FDIFFD = FUNCTN(IC,JSTOL,KCP4)  - FUNCTN(IC,JSTOL,KCM4) 
     +             - FUNCTN(IC,JSTOM8,KCP4) + FUNCTN(IC,JSTOM8,KCM4) 
            FDERIV(IC,JSTOM4,KC) = ACF5YZ*FDIFFA
     +                           + BCF5YZ*FDIFFB
     +                           + CCF5YZ*FDIFFC
     +                           + DCF5YZ*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM2,KCP1) - FUNCTN(IC,JSTOM2,KCM1) 
     +             - FUNCTN(IC,JSTOM4,KCP1) + FUNCTN(IC,JSTOM4,KCM1) 
            FDIFFB = FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2) 
     +             - FUNCTN(IC,JSTOM5,KCP2) + FUNCTN(IC,JSTOM5,KCM2) 
            FDIFFC = FUNCTN(IC,JSTOL,KCP3)  - FUNCTN(IC,JSTOL,KCM3) 
     +             - FUNCTN(IC,JSTOM6,KCP3) + FUNCTN(IC,JSTOM6,KCM3) 
            FDERIV(IC,JSTOM3,KC) = ACF4YZ*FDIFFA
     +                           + BCF4YZ*FDIFFB
     +                           + CCF4YZ*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1) 
     +             - FUNCTN(IC,JSTOM3,KCP1) + FUNCTN(IC,JSTOM3,KCM1) 
            FDIFFB = FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2) 
     +             - FUNCTN(IC,JSTOM4,KCP2) + FUNCTN(IC,JSTOM4,KCM2) 
            FDERIV(IC,JSTOM2,KC) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOL,KCP1)  + FUNCTN(IC,JSTOL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOL,KCP2)  + FUNCTN(IC,JSTOL,KCM2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM2,KCP1) + FUNCTN(IC,JSTOM2,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM2,KCP2) + FUNCTN(IC,JSTOM2,KCM2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM3,KCP1) + FUNCTN(IC,JSTOM3,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM3,KCP2) + FUNCTN(IC,JSTOM3,KCM2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM4,KCP1) + FUNCTN(IC,JSTOM4,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM4,KCP2) + FUNCTN(IC,JSTOM4,KCM2))
            FDERIV(IC,JSTOM1,KC) = ACF2YZ*FDIFFA
     +                           + BCF2YZ*FDIFFB
     +                           + CCF2YZ*FDIFFC
     +                           + DCF2YZ*FDIFFD

C           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM1,KCP1) + FUNCTN(IC,JSTOM1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM1,KCP2) + FUNCTN(IC,JSTOM1,KCM2))
            FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM2,KCP1) + FUNCTN(IC,JSTOM2,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM2,KCP2) + FUNCTN(IC,JSTOM2,KCM2))
            FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM3,KCP1) + FUNCTN(IC,JSTOM3,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM3,KCP2) + FUNCTN(IC,JSTOM3,KCM2))
            FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM4,KCP1) + FUNCTN(IC,JSTOM4,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM4,KCP2) + FUNCTN(IC,JSTOM4,KCM2))
            FDERIV(IC,JSTOL,KC) = ACF1YZ*FDIFFA
     +                          + BCF1YZ*FDIFFB
     +                          + CCF1YZ*FDIFFC
     +                          + DCF1YZ*FDIFFD

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     LH END Z-DIRECTION
C     ==================
      IF(NENDZL.EQ.NBOUND)THEN

C       TAKE SECOND XZ-DERIVATIVE IN Z-LEFT INNER HALO
C       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        JCM3 = JSTART-4
        JCM2 = JSTART-3
        JCM1 = JSTART-2
        JCCC = JSTART-1
        JCP1 = JSTART
        JCP2 = JSTART+1
        JCP3 = JSTART+2
        JCP4 = JSTART+3

        DO JC = JSTART,JFINIS

          JCM4 = JCM3
          JCM3 = JCM2
          JCM2 = JCM1
          JCM1 = JCCC
          JCCC = JCP1
          JCP1 = JCP2
          JCP2 = JCP3
          JCP3 = JCP4
          JCP4 = JC+4

          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP1) - FUNCTN(IC,JCM1,KSTAP1)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP1) - FUNCTN(IC,JCM2,KSTAP1)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
            FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP2) - FUNCTN(IC,JCM1,KSTAP2)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP2) - FUNCTN(IC,JCM2,KSTAP2)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
            FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP3) - FUNCTN(IC,JCM1,KSTAP3)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP3) - FUNCTN(IC,JCM2,KSTAP3)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
            FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP4) - FUNCTN(IC,JCM1,KSTAP4)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP4) - FUNCTN(IC,JCM2,KSTAP4)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
            FDERIV(IC,JC,KSTAL) = ACF1YZ*FDIFFA
     +                          + BCF1YZ*FDIFFB
     +                          + CCF1YZ*FDIFFC
     +                          + DCF1YZ*FDIFFD
     
C           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAL)  - FUNCTN(IC,JCM1,KSTAL)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAL)  - FUNCTN(IC,JCM2,KSTAL)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
            FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP2) - FUNCTN(IC,JCM1,KSTAP2)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP2) - FUNCTN(IC,JCM2,KSTAP2)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
            FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP3) - FUNCTN(IC,JCM1,KSTAP3)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP3) - FUNCTN(IC,JCM2,KSTAP3)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
            FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP4) - FUNCTN(IC,JCM1,KSTAP4)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP4) - FUNCTN(IC,JCM2,KSTAP4)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
            FDERIV(IC,JC,KSTAP1) = ACF2YZ*FDIFFA
     +                           + BCF2YZ*FDIFFB
     +                           + CCF2YZ*FDIFFC
     +                           + DCF2YZ*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JCP1,KSTAP3) - FUNCTN(IC,JCM1,KSTAP3) 
     +             - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1) 
            FDIFFB = FUNCTN(IC,JCP2,KSTAP4) - FUNCTN(IC,JCM2,KSTAP4) 
     +             - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL) 
            FDERIV(IC,JC,KSTAP2) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JCP1,KSTAP4) - FUNCTN(IC,JCM1,KSTAP4) 
     +             - FUNCTN(IC,JCP1,KSTAP2) + FUNCTN(IC,JCM1,KSTAP2) 
            FDIFFB = FUNCTN(IC,JCP2,KSTAP5) - FUNCTN(IC,JCM2,KSTAP5) 
     +             - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1) 
            FDIFFC = FUNCTN(IC,JCP3,KSTAP6) - FUNCTN(IC,JCM3,KSTAP6) 
     +             - FUNCTN(IC,JCP3,KSTAL)  + FUNCTN(IC,JCM3,KSTAL) 
            FDERIV(IC,JC,KSTAP3) = ACF4YZ*FDIFFA
     +                           + BCF4YZ*FDIFFB
     +                           + CCF4YZ*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JCP1,KSTAP5) - FUNCTN(IC,JCM1,KSTAP5) 
     +             - FUNCTN(IC,JCP1,KSTAP3) + FUNCTN(IC,JCM1,KSTAP3) 
            FDIFFB = FUNCTN(IC,JCP2,KSTAP6) - FUNCTN(IC,JCM2,KSTAP6) 
     +             - FUNCTN(IC,JCP2,KSTAP2) + FUNCTN(IC,JCM2,KSTAP2) 
            FDIFFC = FUNCTN(IC,JCP3,KSTAP7) - FUNCTN(IC,JCM3,KSTAP7) 
     +             - FUNCTN(IC,JCP3,KSTAP1) + FUNCTN(IC,JCM3,KSTAP1) 
            FDIFFD = FUNCTN(IC,JCP4,KSTAP8) - FUNCTN(IC,JCM4,KSTAP8) 
     +             - FUNCTN(IC,JCP4,KSTAL)  + FUNCTN(IC,JCM4,KSTAL) 
            FDERIV(IC,JC,KSTAP4) = ACF5YZ*FDIFFA
     +                           + BCF5YZ*FDIFFB
     +                           + CCF5YZ*FDIFFC
     +                           + DCF5YZ*FDIFFD
      
          ENDDO
        ENDDO

C       LH IN Y LH IN Z CORNER
C       ======================
        IF(NENDYL.EQ.NBOUND)THEN

          DO IC = ISTAL,ISTOL

C           LH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(IC,JSTAP1,KSTAP1) - FUNCTN(IC,JSTAP1,KSTAL)
     +             - FUNCTN(IC,JSTAL,KSTAP1)  + FUNCTN(IC,JSTAL,KSTAL)
            FDIFFB = FUNCTN(IC,JSTAP2,KSTAP2) - FUNCTN(IC,JSTAP2,KSTAL)
     +             - FUNCTN(IC,JSTAL,KSTAP2)  + FUNCTN(IC,JSTAL,KSTAL)
            FDIFFC = FUNCTN(IC,JSTAP3,KSTAP3) - FUNCTN(IC,JSTAP3,KSTAL)
     +             - FUNCTN(IC,JSTAL,KSTAP3)  + FUNCTN(IC,JSTAL,KSTAL)
            FDIFFD = FUNCTN(IC,JSTAP4,KSTAP4) - FUNCTN(IC,JSTAP4,KSTAL)
     +             - FUNCTN(IC,JSTAL,KSTAP4)  + FUNCTN(IC,JSTAL,KSTAL)
            FDERIV(IC,JSTAL,KSTAL) = ACC1YZ*FDIFFA
     +                             + BCC1YZ*FDIFFB
     +                             + CCC1YZ*FDIFFC
     +                             + DCC1YZ*FDIFFD

C           LH+1 LH+1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTAL,KSTAL)   - FUNCTN(IC,JSTAL,KSTAP1)
     +             - FUNCTN(IC,JSTAP1,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAP1)
            FDIFFB = FUNCTN(IC,JSTAP2,KSTAP2) - FUNCTN(IC,JSTAP2,KSTAP1)
     +             - FUNCTN(IC,JSTAP1,KSTAP2) + FUNCTN(IC,JSTAP1,KSTAP1)
            FDIFFC = FUNCTN(IC,JSTAP3,KSTAP3) - FUNCTN(IC,JSTAP3,KSTAP1)
     +             - FUNCTN(IC,JSTAP1,KSTAP3) + FUNCTN(IC,JSTAP1,KSTAP1)
            FDIFFD = FUNCTN(IC,JSTAP4,KSTAP4) - FUNCTN(IC,JSTAP4,KSTAP1)
     +             - FUNCTN(IC,JSTAP1,KSTAP4) + FUNCTN(IC,JSTAP1,KSTAP1)
            FDERIV(IC,JSTAP1,KSTAP1) = ACC2YZ*FDIFFA
     +                               + BCC2YZ*FDIFFB
     +                               + CCC2YZ*FDIFFC
     +                               + DCC2YZ*FDIFFD

C           LH LH+1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTAP1,KSTAL)  - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAL,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP1,KSTAP2) - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP2)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP1,KSTAP3) - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP3)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP1,KSTAP4) - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP4)  + FUNCTN(IC,JSTAL,KSTAP1))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTAP2,KSTAL)  - FUNCTN(IC,JSTAP2,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAL,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTAP2) - FUNCTN(IC,JSTAP2,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP2)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP2,KSTAP3) - FUNCTN(IC,JSTAP2,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP3)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP2,KSTAP4) - FUNCTN(IC,JSTAP2,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP4)  + FUNCTN(IC,JSTAL,KSTAP1))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTAP3,KSTAL)  - FUNCTN(IC,JSTAP3,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAL,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP3,KSTAP2) - FUNCTN(IC,JSTAP3,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP2)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTAP3) - FUNCTN(IC,JSTAP3,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP3)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP3,KSTAP4) - FUNCTN(IC,JSTAP3,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP4)  + FUNCTN(IC,JSTAL,KSTAP1))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTAP4,KSTAL)  - FUNCTN(IC,JSTAP4,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAL,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP4,KSTAP2) - FUNCTN(IC,JSTAP4,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP2)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP4,KSTAP3) - FUNCTN(IC,JSTAP4,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP3)  + FUNCTN(IC,JSTAL,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTAP4) - FUNCTN(IC,JSTAP4,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAP4)  + FUNCTN(IC,JSTAL,KSTAP1))
            FDERIV(IC,JSTAL,KSTAP1) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           LH+1 LH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTAP1)  - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAP1,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTAP1) - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAP2,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTAP1) - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAP3,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTAP1) - FUNCTN(IC,JSTAP1,KSTAP1)
     +            - FUNCTN(IC,JSTAP4,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTAP2)  - FUNCTN(IC,JSTAP1,KSTAP2)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAP1,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTAP2) - FUNCTN(IC,JSTAP1,KSTAP2)
     +            - FUNCTN(IC,JSTAP2,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTAP2) - FUNCTN(IC,JSTAP1,KSTAP2)
     +            - FUNCTN(IC,JSTAP3,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTAP2) - FUNCTN(IC,JSTAP1,KSTAP2)
     +            - FUNCTN(IC,JSTAP4,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTAP3)  - FUNCTN(IC,JSTAP1,KSTAP3)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAP1,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTAP3) - FUNCTN(IC,JSTAP1,KSTAP3)
     +            - FUNCTN(IC,JSTAP2,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTAP3) - FUNCTN(IC,JSTAP1,KSTAP3)
     +            - FUNCTN(IC,JSTAP3,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTAP3) - FUNCTN(IC,JSTAP1,KSTAP3)
     +            - FUNCTN(IC,JSTAP4,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTAP4)  - FUNCTN(IC,JSTAP1,KSTAP4)
     +            - FUNCTN(IC,JSTAL,KSTAL)   + FUNCTN(IC,JSTAP1,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTAP4) - FUNCTN(IC,JSTAP1,KSTAP4)
     +            - FUNCTN(IC,JSTAP2,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTAP4) - FUNCTN(IC,JSTAP1,KSTAP4)
     +            - FUNCTN(IC,JSTAP3,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTAP4) - FUNCTN(IC,JSTAP1,KSTAP4)
     +            - FUNCTN(IC,JSTAP4,KSTAL)  + FUNCTN(IC,JSTAP1,KSTAL))
            FDERIV(IC,JSTAP1,KSTAL) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           LH EDGE IN Z
            DO JC = JSTAP2,JSTAP4

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP1) - FUNCTN(IC,JCM1,KSTAP1)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP1) - FUNCTN(IC,JCM2,KSTAP1)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP2) - FUNCTN(IC,JCM1,KSTAP2)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP2) - FUNCTN(IC,JCM2,KSTAP2)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP3) - FUNCTN(IC,JCM1,KSTAP3)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP3) - FUNCTN(IC,JCM2,KSTAP3)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP4) - FUNCTN(IC,JCM1,KSTAP4)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP4) - FUNCTN(IC,JCM2,KSTAP4)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDERIV(IC,JC,KSTAL) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAL)  - FUNCTN(IC,JCM1,KSTAL)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAL)  - FUNCTN(IC,JCM2,KSTAL)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP2) - FUNCTN(IC,JCM1,KSTAP2)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP2) - FUNCTN(IC,JCM2,KSTAP2)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP3) - FUNCTN(IC,JCM1,KSTAP3)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP3) - FUNCTN(IC,JCM2,KSTAP3)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP4) - FUNCTN(IC,JCM1,KSTAP4)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP4) - FUNCTN(IC,JCM2,KSTAP4)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDERIV(IC,JC,KSTAP1) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

            ENDDO

C           LH EDGE IN Y
            DO KC = KSTAP2,KSTAP4

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTAP1,KCP1) - FUNCTN(IC,JSTAP1,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP1,KCP2) - FUNCTN(IC,JSTAP1,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTAP2,KCP1) - FUNCTN(IC,JSTAP2,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP2,KCP2) - FUNCTN(IC,JSTAP2,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTAP3,KCP1) - FUNCTN(IC,JSTAP3,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP3,KCP2) - FUNCTN(IC,JSTAP3,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTAP4,KCP1) - FUNCTN(IC,JSTAP4,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP4,KCP2) - FUNCTN(IC,JSTAP4,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDERIV(IC,JSTAL,KC) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTAL,KCP1)  - FUNCTN(IC,JSTAL,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAL,KCP2)  - FUNCTN(IC,JSTAL,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTAP2,KCP1) - FUNCTN(IC,JSTAP2,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP2,KCP2) - FUNCTN(IC,JSTAP2,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTAP3,KCP1) - FUNCTN(IC,JSTAP3,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP3,KCP2) - FUNCTN(IC,JSTAP3,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTAP4,KCP1) - FUNCTN(IC,JSTAP4,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP4,KCP2) - FUNCTN(IC,JSTAP4,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDERIV(IC,JSTAP1,KC) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTAP2,KSTAP4

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              JS = 0
              DO JC = JSTAP2,JSTAP4

                JS = JS+1
                JCM2 = JC-2
                JCM1 = JC-1
                JCP1 = JC+1
                JCP2 = JC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(IC,JCP1,KCP1) - FUNCTN(IC,JCP1,KCM1)
     +                 - FUNCTN(IC,JCM1,KCP1) + FUNCTN(IC,JCM1,KCM1)
                FDIFFB = FUNCTN(IC,JCP2,KCP2) - FUNCTN(IC,JCP2,KCM2)
     +                 - FUNCTN(IC,JCM2,KCP2) + FUNCTN(IC,JCM2,KCM2)
                FDERIV(IC,JC,KC) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB
                FSTORA(JS,KS) = FDIFFA
                FSTORB(JS,KS) = FDIFFB
 
              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 1
            DO KC = KSTAP3,KSTAP4

              KSM1 = KS
              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              JS = 1
              DO JC = JSTAP3,JSTAP4

                JSM1 = JS
                JS = JS+1
                JCM3 = JC-3
                JCP3 = JC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(IC,JCP3,KCP3) - FUNCTN(IC,JCP3,KCM3)
     +                 - FUNCTN(IC,JCM3,KCP3) + FUNCTN(IC,JCM3,KCM3)
                FDERIV(IC,JC,KC) = ACF4YZ*FSTORA(JS,KS)
     +                           + BCF4YZ*FSTORB(JS,KS)
     +                           + CCF4YZ*FDIFFC
                FSTORC(JSM1,KSM1) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 3
            JS = 3
            KSM1 = 2
            JSM1 = 2
            KC = KSTAP4
            JC = JSTAP4
            KCM4 = KC-4
            KCP4 = KC+4
            JCM4 = JC-4
            JCP4 = JC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(IC,JCP4,KCP4) - FUNCTN(IC,JCP4,KCM4)
     +             - FUNCTN(IC,JCM4,KCP4) + FUNCTN(IC,JCM4,KCM4)
            FDERIV(IC,JC,KC) = ACF5YZ*FSTORA(JS,KS)
     +                       + BCF5YZ*FSTORB(JS,KS)
     +                       + CCF5YZ*FSTORC(JSM1,KSM1)
     +                       + DCF5YZ*FDIFFD

          ENDDO

        ENDIF

C       RH IN Y LH IN Z CORNER
C       ======================
        IF(NENDYR.EQ.NBOUND)THEN

          DO IC = ISTAL,ISTOL

C           RH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(IC,JSTOL,KSTAP1)  - FUNCTN(IC,JSTOL,KSTAL)
     +             - FUNCTN(IC,JSTOM1,KSTAP1) + FUNCTN(IC,JSTOM1,KSTAL)
            FDIFFB = FUNCTN(IC,JSTOL,KSTAP2)  - FUNCTN(IC,JSTOL,KSTAL)
     +             - FUNCTN(IC,JSTOM2,KSTAP2) + FUNCTN(IC,JSTOM2,KSTAL)
            FDIFFC = FUNCTN(IC,JSTOL,KSTAP3)  - FUNCTN(IC,JSTOL,KSTAL)
     +             - FUNCTN(IC,JSTOM3,KSTAP3) + FUNCTN(IC,JSTOM3,KSTAL)
            FDIFFD = FUNCTN(IC,JSTOL,KSTAP4)  - FUNCTN(IC,JSTOL,KSTAL)
     +             - FUNCTN(IC,JSTOM4,KSTAP4) + FUNCTN(IC,JSTOM4,KSTAL)
            FDERIV(IC,JSTOL,KSTAL) = ACC1YZ*FDIFFA
     +                             + BCC1YZ*FDIFFB
     +                             + CCC1YZ*FDIFFC
     +                             + DCC1YZ*FDIFFD

C           RH-1 LH+1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTOM1,KSTAL)  - FUNCTN(IC,JSTOM1,KSTAP1)
     +             - FUNCTN(IC,JSTOL,KSTAL)   + FUNCTN(IC,JSTOL,KSTAP1)
            FDIFFB = FUNCTN(IC,JSTOM1,KSTAP2) - FUNCTN(IC,JSTOM1,KSTAP1)
     +             - FUNCTN(IC,JSTOM2,KSTAP2) + FUNCTN(IC,JSTOM2,KSTAP1)
            FDIFFC = FUNCTN(IC,JSTOM1,KSTAP3) - FUNCTN(IC,JSTOM1,KSTAP1)
     +             - FUNCTN(IC,JSTOM3,KSTAP3) + FUNCTN(IC,JSTOM3,KSTAP1)
            FDIFFD = FUNCTN(IC,JSTOM1,KSTAP4) - FUNCTN(IC,JSTOM1,KSTAP1)
     +             - FUNCTN(IC,JSTOM4,KSTAP4) + FUNCTN(IC,JSTOM4,KSTAP1)
            FDERIV(IC,JSTOM1,KSTAP1) = ACC2YZ*FDIFFA
     +                               + BCC2YZ*FDIFFB
     +                               + CCC2YZ*FDIFFC
     +                               + DCC2YZ*FDIFFD

C           RH LH+1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTAL)   - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM1,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOL,KSTAP2)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAP2) + FUNCTN(IC,JSTOM1,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOL,KSTAP3)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAP3) + FUNCTN(IC,JSTOM1,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOL,KSTAP4)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAP4) + FUNCTN(IC,JSTOM1,KSTAP1))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTAL)   - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM2,KSTAL)  + FUNCTN(IC,JSTOM2,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOL,KSTAP2)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM2,KSTAP2) + FUNCTN(IC,JSTOM2,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOL,KSTAP3)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM2,KSTAP3) + FUNCTN(IC,JSTOM2,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOL,KSTAP4)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM2,KSTAP4) + FUNCTN(IC,JSTOM2,KSTAP1))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTAL)   - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM3,KSTAL)  + FUNCTN(IC,JSTOM3,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOL,KSTAP2)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM3,KSTAP2) + FUNCTN(IC,JSTOM3,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOL,KSTAP3)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM3,KSTAP3) + FUNCTN(IC,JSTOM3,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOL,KSTAP4)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM3,KSTAP4) + FUNCTN(IC,JSTOM3,KSTAP1))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTAL)   - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM4,KSTAL)  + FUNCTN(IC,JSTOM4,KSTAP1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOL,KSTAP2)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM4,KSTAP2) + FUNCTN(IC,JSTOM4,KSTAP1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOL,KSTAP3)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM4,KSTAP3) + FUNCTN(IC,JSTOM4,KSTAP1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOL,KSTAP4)  - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM4,KSTAP4) + FUNCTN(IC,JSTOM4,KSTAP1))
            FDERIV(IC,JSTOL,KSTAP1) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           RH-1 LH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTOM1,KSTAP1) - FUNCTN(IC,JSTOL,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOL,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP1) - FUNCTN(IC,JSTOM2,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM2,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP1) - FUNCTN(IC,JSTOM3,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM3,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP1) - FUNCTN(IC,JSTOM4,KSTAP1)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM4,KSTAL))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTOM1,KSTAP2) - FUNCTN(IC,JSTOL,KSTAP2)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOL,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP2) - FUNCTN(IC,JSTOM2,KSTAP2)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM2,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP2) - FUNCTN(IC,JSTOM3,KSTAP2)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM3,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP2) - FUNCTN(IC,JSTOM4,KSTAP2)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM4,KSTAL))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTOM1,KSTAP3) - FUNCTN(IC,JSTOL,KSTAP3)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOL,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP3) - FUNCTN(IC,JSTOM2,KSTAP3)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM2,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP3) - FUNCTN(IC,JSTOM3,KSTAP3)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM3,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP3) - FUNCTN(IC,JSTOM4,KSTAP3)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM4,KSTAL))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTOM1,KSTAP4) - FUNCTN(IC,JSTOL,KSTAP4)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOL,KSTAL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP4) - FUNCTN(IC,JSTOM2,KSTAP4)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM2,KSTAL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP4) - FUNCTN(IC,JSTOM3,KSTAP4)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM3,KSTAL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM1,KSTAP4) - FUNCTN(IC,JSTOM4,KSTAP4)
     +            - FUNCTN(IC,JSTOM1,KSTAL)  + FUNCTN(IC,JSTOM4,KSTAL))
            FDERIV(IC,JSTOM1,KSTAL) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           LH EDGE IN Z
            DO JC = JSTOM4,JSTOM2

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP1) - FUNCTN(IC,JCM1,KSTAP1)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP1) - FUNCTN(IC,JCM2,KSTAP1)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP2) - FUNCTN(IC,JCM1,KSTAP2)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP2) - FUNCTN(IC,JCM2,KSTAP2)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP3) - FUNCTN(IC,JCM1,KSTAP3)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP3) - FUNCTN(IC,JCM2,KSTAP3)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP4) - FUNCTN(IC,JCM1,KSTAP4)
     +                - FUNCTN(IC,JCP1,KSTAL)  + FUNCTN(IC,JCM1,KSTAL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP4) - FUNCTN(IC,JCM2,KSTAP4)
     +                - FUNCTN(IC,JCP2,KSTAL)  + FUNCTN(IC,JCM2,KSTAL))
              FDERIV(IC,JC,KSTAL) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAL)  - FUNCTN(IC,JCM1,KSTAL)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAL)  - FUNCTN(IC,JCM2,KSTAL)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP2) - FUNCTN(IC,JCM1,KSTAP2)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP2) - FUNCTN(IC,JCM2,KSTAP2)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP3) - FUNCTN(IC,JCM1,KSTAP3)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP3) - FUNCTN(IC,JCM2,KSTAP3)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTAP4) - FUNCTN(IC,JCM1,KSTAP4)
     +                - FUNCTN(IC,JCP1,KSTAP1) + FUNCTN(IC,JCM1,KSTAP1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTAP4) - FUNCTN(IC,JCM2,KSTAP4)
     +                - FUNCTN(IC,JCP2,KSTAP1) + FUNCTN(IC,JCM2,KSTAP1))
              FDERIV(IC,JC,KSTAP1) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

            ENDDO

C           RH EDGE IN Y
            DO KC = KSTAP2,KSTAP4

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM1,KCP1) + FUNCTN(IC,JSTOM1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM1,KCP2) + FUNCTN(IC,JSTOM1,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM2,KCP1) + FUNCTN(IC,JSTOM2,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM2,KCP2) + FUNCTN(IC,JSTOM2,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM3,KCP1) + FUNCTN(IC,JSTOM3,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM3,KCP2) + FUNCTN(IC,JSTOM3,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM4,KCP1) + FUNCTN(IC,JSTOM4,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM4,KCP2) + FUNCTN(IC,JSTOM4,KCM2))
              FDERIV(IC,JSTOL,KC) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD
    
C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOL,KCP1)  + FUNCTN(IC,JSTOL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOL,KCP2)  + FUNCTN(IC,JSTOL,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM2,KCP1) + FUNCTN(IC,JSTOM2,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM2,KCP2) + FUNCTN(IC,JSTOM2,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM3,KCP1) + FUNCTN(IC,JSTOM3,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM3,KCP2) + FUNCTN(IC,JSTOM3,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM4,KCP1) + FUNCTN(IC,JSTOM4,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM4,KCP2) + FUNCTN(IC,JSTOM4,KCM2))
              FDERIV(IC,JSTOM1,KC) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTAP2,KSTAP4

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              JS = 0
              DO JC = JSTOM4,JSTOM2

                JS = JS+1
                JCM2 = JC-2
                JCM1 = JC-1
                JCP1 = JC+1
                JCP2 = JC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(IC,JCP1,KCP1) - FUNCTN(IC,JCP1,KCM1)
     +                 - FUNCTN(IC,JCM1,KCP1) + FUNCTN(IC,JCM1,KCM1)
                FDIFFB = FUNCTN(IC,JCP2,KCP2) - FUNCTN(IC,JCP2,KCM2)
     +                 - FUNCTN(IC,JCM2,KCP2) + FUNCTN(IC,JCM2,KCM2)
                FDERIV(IC,JC,KC) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB
                FSTORA(JS,KS) = FDIFFA
                FSTORB(JS,KS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 1
            DO KC = KSTAP3,KSTAP4

              KSM1 = KS
              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              JS = 0
              DO JC = JSTOM4,JSTOM3

                JS = JS+1
                JCM3 = JC-3
                JCP3 = JC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(IC,JCP3,KCP3) - FUNCTN(IC,JCP3,KCM3)
     +                 - FUNCTN(IC,JCM3,KCP3) + FUNCTN(IC,JCM3,KCM3)
                FDERIV(IC,JC,KC) = ACF4YZ*FSTORA(JS,KS)
     +                           + BCF4YZ*FSTORB(JS,KS)
     +                           + CCF4YZ*FDIFFC
                FSTORC(JS,KSM1) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 3
            JS = 1
            KSM1 = 2
            KC = KSTAP4
            JC = JSTOM4
            KCM4 = KC-4
            KCP4 = KC+4
            JCM4 = JC-4
            JCP4 = JC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(IC,JCP4,KCP4) - FUNCTN(IC,JCP4,KCM4)
     +             - FUNCTN(IC,JCM4,KCP4) + FUNCTN(IC,JCM4,KCM4)
            FDERIV(IC,JC,KC) = ACF5YZ*FSTORA(JS,KS)
     +                       + BCF5YZ*FSTORB(JS,KS)
     +                       + CCF5YZ*FSTORC(JS,KSM1)
     +                       + DCF5YZ*FDIFFD

          ENDDO

        ENDIF

      ENDIF 

C     =========================================================================

C     RH END Z-DIRECTION
C     ==================
      IF(NENDZR.EQ.NBOUND)THEN

C       TAKE SECOND YZ-DERIVATIVE IN Z-RIGHT INNER HALO
C       EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT

        JCM3 = JSTART-4
        JCM2 = JSTART-3
        JCM1 = JSTART-2
        JCCC = JSTART-1
        JCP1 = JSTART
        JCP2 = JSTART+1
        JCP3 = JSTART+2
        JCP4 = JSTART+3

        DO JC = JSTART,JFINIS

          JCM4 = JCM3
          JCM3 = JCM2
          JCM2 = JCM1
          JCM1 = JCCC
          JCCC = JCP1
          JCP1 = JCP2
          JCP2 = JCP3
          JCP3 = JCP4
          JCP4 = JC+4

          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JCP1,KSTOM3) - FUNCTN(IC,JCM1,KSTOM3) 
     +             - FUNCTN(IC,JCP1,KSTOM5) + FUNCTN(IC,JCM1,KSTOM5) 
            FDIFFB = FUNCTN(IC,JCP2,KSTOM2) - FUNCTN(IC,JCM2,KSTOM2) 
     +             - FUNCTN(IC,JCP2,KSTOM6) + FUNCTN(IC,JCM2,KSTOM6) 
            FDIFFC = FUNCTN(IC,JCP3,KSTOM1) - FUNCTN(IC,JCM3,KSTOM1) 
     +             - FUNCTN(IC,JCP3,KSTOM7) + FUNCTN(IC,JCM3,KSTOM7) 
            FDIFFD = FUNCTN(IC,JCP4,KSTOL)  - FUNCTN(IC,JCM4,KSTOL) 
     +             - FUNCTN(IC,JCP4,KSTOM8) + FUNCTN(IC,JCM4,KSTOM8) 
            FDERIV(IC,JC,KSTOM4) = ACF5YZ*FDIFFA
     +                           + BCF5YZ*FDIFFB
     +                           + CCF5YZ*FDIFFC
     +                           + DCF5YZ*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JCP1,KSTOM2) - FUNCTN(IC,JCM1,KSTOM2) 
     +             - FUNCTN(IC,JCP1,KSTOM4) + FUNCTN(IC,JCM1,KSTOM4) 
            FDIFFB = FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1) 
     +             - FUNCTN(IC,JCP2,KSTOM5) + FUNCTN(IC,JCM2,KSTOM5) 
            FDIFFC = FUNCTN(IC,JCP3,KSTOL)  - FUNCTN(IC,JCM3,KSTOL) 
     +             - FUNCTN(IC,JCP3,KSTOM6) + FUNCTN(IC,JCM3,KSTOM6) 
            FDERIV(IC,JC,KSTOM3) = ACF4YZ*FDIFFA
     +                           + BCF4YZ*FDIFFB
     +                           + CCF4YZ*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1) 
     +             - FUNCTN(IC,JCP1,KSTOM3) + FUNCTN(IC,JCM1,KSTOM3) 
            FDIFFB = FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL) 
     +             - FUNCTN(IC,JCP2,KSTOM4) + FUNCTN(IC,JCM2,KSTOM4) 
            FDERIV(IC,JC,KSTOM2) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOL)  + FUNCTN(IC,JCM1,KSTOL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOL)  + FUNCTN(IC,JCM2,KSTOL))
            FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM2) + FUNCTN(IC,JCM1,KSTOM2))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM2) + FUNCTN(IC,JCM2,KSTOM2))
            FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM3) + FUNCTN(IC,JCM1,KSTOM3))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM3) + FUNCTN(IC,JCM2,KSTOM3))
            FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM4) + FUNCTN(IC,JCM1,KSTOM4))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM4) + FUNCTN(IC,JCM2,KSTOM4))
            FDERIV(IC,JC,KSTOM1) = ACF2YZ*FDIFFA
     +                           + BCF2YZ*FDIFFB
     +                           + CCF2YZ*FDIFFC
     +                           + DCF2YZ*FDIFFD

C           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM1) + FUNCTN(IC,JCM1,KSTOM1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM1) + FUNCTN(IC,JCM2,KSTOM1))
            FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM2) + FUNCTN(IC,JCM1,KSTOM2))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM2) + FUNCTN(IC,JCM2,KSTOM2))
            FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM3) + FUNCTN(IC,JCM1,KSTOM3))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM3) + FUNCTN(IC,JCM2,KSTOM3))
            FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM4) + FUNCTN(IC,JCM1,KSTOM4))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM4) + FUNCTN(IC,JCM2,KSTOM4))
            FDERIV(IC,JC,KSTOL) = ACF1YZ*FDIFFA
     +                          + BCF1YZ*FDIFFB
     +                          + CCF1YZ*FDIFFC
     +                          + DCF1YZ*FDIFFD

          ENDDO
        ENDDO

C       LH IN Y RH IN Z CORNER
C       ======================
        IF(NENDYL.EQ.NBOUND)THEN

          DO IC = ISTAL,ISTOL

C           LH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(IC,JSTAP1,KSTOL) - FUNCTN(IC,JSTAP1,KSTOM1)
     +             - FUNCTN(IC,JSTAL,KSTOL)  + FUNCTN(IC,JSTAL,KSTOM1) 
            FDIFFB = FUNCTN(IC,JSTAP2,KSTOL) - FUNCTN(IC,JSTAP2,KSTOM2)
     +             - FUNCTN(IC,JSTAL,KSTOL)  + FUNCTN(IC,JSTAL,KSTOM2) 
            FDIFFC = FUNCTN(IC,JSTAP3,KSTOL) - FUNCTN(IC,JSTAP3,KSTOM3)
     +             - FUNCTN(IC,JSTAL,KSTOL)  + FUNCTN(IC,JSTAL,KSTOM3) 
            FDIFFD = FUNCTN(IC,JSTAP4,KSTOL) - FUNCTN(IC,JSTAP4,KSTOM4)
     +             - FUNCTN(IC,JSTAL,KSTOL)  + FUNCTN(IC,JSTAL,KSTOM4) 
            FDERIV(IC,JSTAL,KSTOL) = ACC1YZ*FDIFFA
     +                             + BCC1YZ*FDIFFB
     +                             + CCC1YZ*FDIFFC
     +                             + DCC1YZ*FDIFFD

C           LH+1 RH-1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTAL,KSTOM1)  - FUNCTN(IC,JSTAL,KSTOL)
     +             - FUNCTN(IC,JSTAP1,KSTOM1) + FUNCTN(IC,JSTAP1,KSTOL)
            FDIFFB = FUNCTN(IC,JSTAP2,KSTOM1) - FUNCTN(IC,JSTAP2,KSTOM2)
     +             - FUNCTN(IC,JSTAP1,KSTOM1) + FUNCTN(IC,JSTAP1,KSTOM2)
            FDIFFC = FUNCTN(IC,JSTAP3,KSTOM1) - FUNCTN(IC,JSTAP3,KSTOM3)
     +             - FUNCTN(IC,JSTAP1,KSTOM1) + FUNCTN(IC,JSTAP1,KSTOM3)
            FDIFFD = FUNCTN(IC,JSTAP4,KSTOM1) - FUNCTN(IC,JSTAP4,KSTOM4)
     +             - FUNCTN(IC,JSTAP1,KSTOM1) + FUNCTN(IC,JSTAP1,KSTOM4)
            FDERIV(IC,JSTAP1,KSTOM1) = ACC2YZ*FDIFFA
     +                               + BCC2YZ*FDIFFB
     +                               + CCC2YZ*FDIFFC
     +                               + DCC2YZ*FDIFFD

C           LH RH-1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTAP1,KSTOM1) - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP1,KSTOM1) - FUNCTN(IC,JSTAP1,KSTOM2)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM2))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP1,KSTOM1) - FUNCTN(IC,JSTAP1,KSTOM3)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM3))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP1,KSTOM1) - FUNCTN(IC,JSTAP1,KSTOM4)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM4))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTAP2,KSTOM1) - FUNCTN(IC,JSTAP2,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTOM1) - FUNCTN(IC,JSTAP2,KSTOM2)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM2))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP2,KSTOM1) - FUNCTN(IC,JSTAP2,KSTOM3)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM3))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP2,KSTOM1) - FUNCTN(IC,JSTAP2,KSTOM4)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM4))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTAP3,KSTOM1) - FUNCTN(IC,JSTAP3,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP3,KSTOM1) - FUNCTN(IC,JSTAP3,KSTOM2)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM2))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTOM1) - FUNCTN(IC,JSTAP3,KSTOM3)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM3))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP3,KSTOM1) - FUNCTN(IC,JSTAP3,KSTOM4)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM4))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTAP4,KSTOM1) - FUNCTN(IC,JSTAP4,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP4,KSTOM1) - FUNCTN(IC,JSTAP4,KSTOM2)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM2))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP4,KSTOM1) - FUNCTN(IC,JSTAP4,KSTOM3)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM3))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTOM1) - FUNCTN(IC,JSTAP4,KSTOM4)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAL,KSTOM4))
            FDERIV(IC,JSTAL,KSTOM1) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           LH+1 RH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTOL)   - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM1)  + FUNCTN(IC,JSTAP1,KSTOM1))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP2,KSTOM1) + FUNCTN(IC,JSTAP1,KSTOM1))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP3,KSTOM1) + FUNCTN(IC,JSTAP1,KSTOM1))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP4,KSTOM1) + FUNCTN(IC,JSTAP1,KSTOM1))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTOL)   - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM2)  + FUNCTN(IC,JSTAP1,KSTOM2))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP2,KSTOM2) + FUNCTN(IC,JSTAP1,KSTOM2))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP3,KSTOM2) + FUNCTN(IC,JSTAP1,KSTOM2))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP4,KSTOM2) + FUNCTN(IC,JSTAP1,KSTOM2))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTOL)   - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM3)  + FUNCTN(IC,JSTAP1,KSTOM3))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP2,KSTOM3) + FUNCTN(IC,JSTAP1,KSTOM3))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP3,KSTOM3) + FUNCTN(IC,JSTAP1,KSTOM3))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP4,KSTOM3) + FUNCTN(IC,JSTAP1,KSTOM3))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTAL,KSTOL)   - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAL,KSTOM4)  + FUNCTN(IC,JSTAP1,KSTOM4))
     +    + BCF2YZ*(FUNCTN(IC,JSTAP2,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP2,KSTOM4) + FUNCTN(IC,JSTAP1,KSTOM4))
     +    + CCF2YZ*(FUNCTN(IC,JSTAP3,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP3,KSTOM4) + FUNCTN(IC,JSTAP1,KSTOM4))
     +    + DCF2YZ*(FUNCTN(IC,JSTAP4,KSTOL)  - FUNCTN(IC,JSTAP1,KSTOL)
     +            - FUNCTN(IC,JSTAP4,KSTOM4) + FUNCTN(IC,JSTAP1,KSTOM4))
            FDERIV(IC,JSTAP1,KSTOL) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           RH EDGE IN Z
            DO JC = JSTAP2,JSTAP4

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOL)  + FUNCTN(IC,JCM1,KSTOL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOL)  + FUNCTN(IC,JCM2,KSTOL))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM2) + FUNCTN(IC,JCM1,KSTOM2))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM2) + FUNCTN(IC,JCM2,KSTOM2))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM3) + FUNCTN(IC,JCM1,KSTOM3))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM3) + FUNCTN(IC,JCM2,KSTOM3))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM4) + FUNCTN(IC,JCM1,KSTOM4))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM4) + FUNCTN(IC,JCM2,KSTOM4))
              FDERIV(IC,JC,KSTOM1) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM1) + FUNCTN(IC,JCM1,KSTOM1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM1) + FUNCTN(IC,JCM2,KSTOM1))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM2) + FUNCTN(IC,JCM1,KSTOM2))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM2) + FUNCTN(IC,JCM2,KSTOM2))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM3) + FUNCTN(IC,JCM1,KSTOM3))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM3) + FUNCTN(IC,JCM2,KSTOM3))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM4) + FUNCTN(IC,JCM1,KSTOM4))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM4) + FUNCTN(IC,JCM2,KSTOM4))
              FDERIV(IC,JC,KSTOL) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD

            ENDDO

C           LH EDGE IN Y
            DO KC = KSTOM4,KSTOM2

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTAP1,KCP1) - FUNCTN(IC,JSTAP1,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP1,KCP2) - FUNCTN(IC,JSTAP1,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTAP2,KCP1) - FUNCTN(IC,JSTAP2,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP2,KCP2) - FUNCTN(IC,JSTAP2,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTAP3,KCP1) - FUNCTN(IC,JSTAP3,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP3,KCP2) - FUNCTN(IC,JSTAP3,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTAP4,KCP1) - FUNCTN(IC,JSTAP4,KCM1)
     +                - FUNCTN(IC,JSTAL,KCP1)  + FUNCTN(IC,JSTAL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP4,KCP2) - FUNCTN(IC,JSTAP4,KCM2)
     +                - FUNCTN(IC,JSTAL,KCP2)  + FUNCTN(IC,JSTAL,KCM2))
              FDERIV(IC,JSTAL,KC) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTAL,KCP1)  - FUNCTN(IC,JSTAL,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAL,KCP2)  - FUNCTN(IC,JSTAL,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTAP2,KCP1) - FUNCTN(IC,JSTAP2,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP2,KCP2) - FUNCTN(IC,JSTAP2,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTAP3,KCP1) - FUNCTN(IC,JSTAP3,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP3,KCP2) - FUNCTN(IC,JSTAP3,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTAP4,KCP1) - FUNCTN(IC,JSTAP4,KCM1)
     +                - FUNCTN(IC,JSTAP1,KCP1) + FUNCTN(IC,JSTAP1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTAP4,KCP2) - FUNCTN(IC,JSTAP4,KCM2)
     +                - FUNCTN(IC,JSTAP1,KCP2) + FUNCTN(IC,JSTAP1,KCM2))
              FDERIV(IC,JSTAP1,KC) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM2

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              JS = 0
              DO JC = JSTAP2,JSTAP4

                JS = JS+1
                JCM2 = JC-2
                JCM1 = JC-1
                JCP1 = JC+1
                JCP2 = JC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(IC,JCP1,KCP1) - FUNCTN(IC,JCP1,KCM1)
     +                 - FUNCTN(IC,JCM1,KCP1) + FUNCTN(IC,JCM1,KCM1)
                FDIFFB = FUNCTN(IC,JCP2,KCP2) - FUNCTN(IC,JCP2,KCM2)
     +                 - FUNCTN(IC,JCM2,KCP2) + FUNCTN(IC,JCM2,KCM2)
                FDERIV(IC,JC,KC) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB
                FSTORA(JS,KS) = FDIFFA
                FSTORB(JS,KS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM3

              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              JS = 1
              DO JC = JSTAP3,JSTAP4

                JSM1 = JS
                JS = JS+1
                JCM3 = JC-3
                JCP3 = JC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(IC,JCP3,KCP3) - FUNCTN(IC,JCP3,KCM3)
     +                 - FUNCTN(IC,JCM3,KCP3) + FUNCTN(IC,JCM3,KCM3)
                FDERIV(IC,JC,KC) = ACF4YZ*FSTORA(JS,KS)
     +                           + BCF4YZ*FSTORB(JS,KS)
     +                           + CCF4YZ*FDIFFC
                FSTORC(JSM1,KS) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 1
            JS = 3
            JSM1 = 2
            KC = KSTOM4
            JC = JSTAP4
            KCM4 = KC-4
            KCP4 = KC+4
            JCM4 = JC-4
            JCP4 = JC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(IC,JCP4,KCP4) - FUNCTN(IC,JCP4,KCM4)
     +             - FUNCTN(IC,JCM4,KCP4) + FUNCTN(IC,JCM4,KCM4)
            FDERIV(IC,JC,KC) = ACF5YZ*FSTORA(JS,KS)
     +                       + BCF5YZ*FSTORB(JS,KS)
     +                       + CCF5YZ*FSTORC(JSM1,KS)
     +                       + DCF5YZ*FDIFFD

          ENDDO

        ENDIF

C       RH IN Y RH IN Z CORNER
C       ======================
        IF(NENDYR.EQ.NBOUND)THEN

          DO IC = ISTAL,ISTOL

C           RH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(IC,JSTOM1,KSTOM1) - FUNCTN(IC,JSTOM1,KSTOL)
     +             - FUNCTN(IC,JSTOL,KSTOM1)  + FUNCTN(IC,JSTOL,KSTOL)
            FDIFFB = FUNCTN(IC,JSTOM2,KSTOM2) - FUNCTN(IC,JSTOM2,KSTOL)
     +             - FUNCTN(IC,JSTOL,KSTOM2)  + FUNCTN(IC,JSTOL,KSTOL)
            FDIFFC = FUNCTN(IC,JSTOM3,KSTOM3) - FUNCTN(IC,JSTOM3,KSTOL)
     +             - FUNCTN(IC,JSTOL,KSTOM3)  + FUNCTN(IC,JSTOL,KSTOL)
            FDIFFD = FUNCTN(IC,JSTOM4,KSTOM4) - FUNCTN(IC,JSTOM4,KSTOL)
     +             - FUNCTN(IC,JSTOL,KSTOM4)  + FUNCTN(IC,JSTOL,KSTOL)
            FDERIV(IC,JSTOL,KSTOL) = ACC1YZ*FDIFFA
     +                             + BCC1YZ*FDIFFB
     +                             + CCC1YZ*FDIFFC
     +                             + DCC1YZ*FDIFFD

C           RH-1 RH-1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTOL,KSTOL)   - FUNCTN(IC,JSTOL,KSTOM1)
     +             - FUNCTN(IC,JSTOM1,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOM1)
            FDIFFB = FUNCTN(IC,JSTOM2,KSTOM2) - FUNCTN(IC,JSTOM2,KSTOM1)
     +             - FUNCTN(IC,JSTOM1,KSTOM2) + FUNCTN(IC,JSTOM1,KSTOM1)
            FDIFFC = FUNCTN(IC,JSTOM3,KSTOM3) - FUNCTN(IC,JSTOM3,KSTOM1)
     +             - FUNCTN(IC,JSTOM1,KSTOM3) + FUNCTN(IC,JSTOM1,KSTOM1)
            FDIFFD = FUNCTN(IC,JSTOM4,KSTOM4) - FUNCTN(IC,JSTOM4,KSTOM1)
     +             - FUNCTN(IC,JSTOM1,KSTOM4) + FUNCTN(IC,JSTOM1,KSTOM1)
            FDERIV(IC,JSTOM1,KSTOM1) = ACC2YZ*FDIFFA
     +                               + BCC2YZ*FDIFFB
     +                               + CCC2YZ*FDIFFC
     +                               + DCC2YZ*FDIFFD

C           RH RH-1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTOM1,KSTOL)  - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOL,KSTOM1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM1,KSTOM2) - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM2)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM1,KSTOM3) - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM3)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM1,KSTOM4) - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM4)  + FUNCTN(IC,JSTOL,KSTOM1))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTOM2,KSTOL)  - FUNCTN(IC,JSTOM2,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOL,KSTOM1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM2,KSTOM2) - FUNCTN(IC,JSTOM2,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM2)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM2,KSTOM3) - FUNCTN(IC,JSTOM2,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM3)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM2,KSTOM4) - FUNCTN(IC,JSTOM2,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM4)  + FUNCTN(IC,JSTOL,KSTOM1))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTOM3,KSTOL)  - FUNCTN(IC,JSTOM3,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOL,KSTOM1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM3,KSTOM2) - FUNCTN(IC,JSTOM3,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM2)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM3,KSTOM3) - FUNCTN(IC,JSTOM3,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM3)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM3,KSTOM4) - FUNCTN(IC,JSTOM3,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM4)  + FUNCTN(IC,JSTOL,KSTOM1))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTOM4,KSTOL)  - FUNCTN(IC,JSTOM4,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOL,KSTOM1))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM4,KSTOM2) - FUNCTN(IC,JSTOM4,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM2)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM4,KSTOM3) - FUNCTN(IC,JSTOM4,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM3)  + FUNCTN(IC,JSTOL,KSTOM1))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM4,KSTOM4) - FUNCTN(IC,JSTOM4,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOM4)  + FUNCTN(IC,JSTOL,KSTOM1))
            FDERIV(IC,JSTOL,KSTOM1) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           RH+1 RH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTOM1)  - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOM1,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM2,KSTOM1) - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOM2,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM3,KSTOM1) - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOM3,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM4,KSTOM1) - FUNCTN(IC,JSTOM1,KSTOM1)
     +            - FUNCTN(IC,JSTOM4,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
            FDIFFB =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTOM2)  - FUNCTN(IC,JSTOM1,KSTOM2)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOM1,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM2,KSTOM2) - FUNCTN(IC,JSTOM1,KSTOM2)
     +            - FUNCTN(IC,JSTOM2,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM3,KSTOM2) - FUNCTN(IC,JSTOM1,KSTOM2)
     +            - FUNCTN(IC,JSTOM3,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM4,KSTOM2) - FUNCTN(IC,JSTOM1,KSTOM2)
     +            - FUNCTN(IC,JSTOM4,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
            FDIFFC =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTOM3)  - FUNCTN(IC,JSTOM1,KSTOM3)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOM1,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM2,KSTOM3) - FUNCTN(IC,JSTOM1,KSTOM3)
     +            - FUNCTN(IC,JSTOM2,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM3,KSTOM3) - FUNCTN(IC,JSTOM1,KSTOM3)
     +            - FUNCTN(IC,JSTOM3,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM4,KSTOM3) - FUNCTN(IC,JSTOM1,KSTOM3)
     +            - FUNCTN(IC,JSTOM4,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
            FDIFFD =
     +      ACF2YZ*(FUNCTN(IC,JSTOL,KSTOM4)  - FUNCTN(IC,JSTOM1,KSTOM4)
     +            - FUNCTN(IC,JSTOL,KSTOL)   + FUNCTN(IC,JSTOM1,KSTOL))
     +    + BCF2YZ*(FUNCTN(IC,JSTOM2,KSTOM4) - FUNCTN(IC,JSTOM1,KSTOM4)
     +            - FUNCTN(IC,JSTOM2,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + CCF2YZ*(FUNCTN(IC,JSTOM3,KSTOM4) - FUNCTN(IC,JSTOM1,KSTOM4)
     +            - FUNCTN(IC,JSTOM3,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
     +    + DCF2YZ*(FUNCTN(IC,JSTOM4,KSTOM4) - FUNCTN(IC,JSTOM1,KSTOM4)
     +            - FUNCTN(IC,JSTOM4,KSTOL)  + FUNCTN(IC,JSTOM1,KSTOL))
            FDERIV(IC,JSTOM1,KSTOL) = ACF1YZ*FDIFFA
     +                              + BCF1YZ*FDIFFB
     +                              + CCF1YZ*FDIFFC
     +                              + DCF1YZ*FDIFFD

C           RH EDGE IN Z
            DO JC = JSTOM4,JSTOM2

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOL)  + FUNCTN(IC,JCM1,KSTOL))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOL)  + FUNCTN(IC,JCM2,KSTOL))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM2) + FUNCTN(IC,JCM1,KSTOM2))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM2) + FUNCTN(IC,JCM2,KSTOM2))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM3) + FUNCTN(IC,JCM1,KSTOM3))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM3) + FUNCTN(IC,JCM2,KSTOM3))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOM1) - FUNCTN(IC,JCM1,KSTOM1)
     +                - FUNCTN(IC,JCP1,KSTOM4) + FUNCTN(IC,JCM1,KSTOM4))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOM1) - FUNCTN(IC,JCM2,KSTOM1)
     +                - FUNCTN(IC,JCP2,KSTOM4) + FUNCTN(IC,JCM2,KSTOM4))
              FDERIV(IC,JC,KSTOM1) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM1) + FUNCTN(IC,JCM1,KSTOM1))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM1) + FUNCTN(IC,JCM2,KSTOM1))
              FDIFFB =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM2) + FUNCTN(IC,JCM1,KSTOM2))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM2) + FUNCTN(IC,JCM2,KSTOM2))
              FDIFFC =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM3) + FUNCTN(IC,JCM1,KSTOM3))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM3) + FUNCTN(IC,JCM2,KSTOM3))
              FDIFFD =
     +          ACOFY1*(FUNCTN(IC,JCP1,KSTOL)  - FUNCTN(IC,JCM1,KSTOL)
     +                - FUNCTN(IC,JCP1,KSTOM4) + FUNCTN(IC,JCM1,KSTOM4))
     +        + BCOFY1*(FUNCTN(IC,JCP2,KSTOL)  - FUNCTN(IC,JCM2,KSTOL)
     +                - FUNCTN(IC,JCP2,KSTOM4) + FUNCTN(IC,JCM2,KSTOM4))
              FDERIV(IC,JC,KSTOL) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD

            ENDDO

C           RH EDGE IN Y
            DO KC = KSTOM4,KSTOM2

              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM1,KCP1) + FUNCTN(IC,JSTOM1,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM1,KCP2) + FUNCTN(IC,JSTOM1,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM2,KCP1) + FUNCTN(IC,JSTOM2,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM2,KCP2) + FUNCTN(IC,JSTOM2,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM3,KCP1) + FUNCTN(IC,JSTOM3,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM3,KCP2) + FUNCTN(IC,JSTOM3,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTOL,KCP1)  - FUNCTN(IC,JSTOL,KCM1)
     +                - FUNCTN(IC,JSTOM4,KCP1) + FUNCTN(IC,JSTOM4,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOL,KCP2)  - FUNCTN(IC,JSTOL,KCM2)
     +                - FUNCTN(IC,JSTOM4,KCP2) + FUNCTN(IC,JSTOM4,KCM2))
              FDERIV(IC,JSTOL,KC) = ACF1YZ*FDIFFA
     +                            + BCF1YZ*FDIFFB
     +                            + CCF1YZ*FDIFFC
     +                            + DCF1YZ*FDIFFD
     
C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOL,KCP1)  + FUNCTN(IC,JSTOL,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOL,KCP2)  + FUNCTN(IC,JSTOL,KCM2))
              FDIFFB =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM2,KCP1) + FUNCTN(IC,JSTOM2,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM2,KCP2) + FUNCTN(IC,JSTOM2,KCM2))
              FDIFFC =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM3,KCP1) + FUNCTN(IC,JSTOM3,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM3,KCP2) + FUNCTN(IC,JSTOM3,KCM2))
              FDIFFD =
     +          ACOFZ1*(FUNCTN(IC,JSTOM1,KCP1) - FUNCTN(IC,JSTOM1,KCM1)
     +                - FUNCTN(IC,JSTOM4,KCP1) + FUNCTN(IC,JSTOM4,KCM1))
     +        + BCOFZ1*(FUNCTN(IC,JSTOM1,KCP2) - FUNCTN(IC,JSTOM1,KCM2)
     +                - FUNCTN(IC,JSTOM4,KCP2) + FUNCTN(IC,JSTOM4,KCM2))
              FDERIV(IC,JSTOM1,KC) = ACF2YZ*FDIFFA
     +                             + BCF2YZ*FDIFFB
     +                             + CCF2YZ*FDIFFC
     +                             + DCF2YZ*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM2

              KS = KS+1
              KCM2 = KC-2
              KCM1 = KC-1
              KCP1 = KC+1
              KCP2 = KC+2

              JS = 0
              DO JC = JSTOM4,JSTOM2

                JS = JS+1
                JCM2 = JC-2
                JCM1 = JC-1
                JCP1 = JC+1
                JCP2 = JC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(IC,JCP1,KCP1) - FUNCTN(IC,JCP1,KCM1)
     +                 - FUNCTN(IC,JCM1,KCP1) + FUNCTN(IC,JCM1,KCM1)
                FDIFFB = FUNCTN(IC,JCP2,KCP2) - FUNCTN(IC,JCP2,KCM2)
     +                 - FUNCTN(IC,JCM2,KCP2) + FUNCTN(IC,JCM2,KCM2)
                FDERIV(IC,JC,KC) = ACF3YZ*FDIFFA
     +                           + BCF3YZ*FDIFFB
                FSTORA(JS,KS) = FDIFFA
                FSTORB(JS,KS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            KS = 0
            DO KC = KSTOM4,KSTOM3

              KS = KS+1
              KCM3 = KC-3
              KCP3 = KC+3

              JS = 0
              DO JC = JSTOM4,JSTOM3

                JS = JS+1
                JCM3 = JC-3
                JCP3 = JC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(IC,JCP3,KCP3) - FUNCTN(IC,JCP3,KCM3)
     +                 - FUNCTN(IC,JCM3,KCP3) + FUNCTN(IC,JCM3,KCM3)
                FDERIV(IC,JC,KC) = ACF4YZ*FSTORA(JS,KS)
     +                           + BCF4YZ*FSTORB(JS,KS)
     +                           + CCF4YZ*FDIFFC
                FSTORC(JS,KS) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            KS = 1
            JS = 1
            KC = KSTOM4
            JC = JSTOM4
            KCM4 = KC-4
            KCP4 = KC+4
            JCM4 = JC-4
            JCP4 = JC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(IC,JCP4,KCP4) - FUNCTN(IC,JCP4,KCM4)
     +             - FUNCTN(IC,JCM4,KCP4) + FUNCTN(IC,JCM4,KCM4)
            FDERIV(IC,JC,KC) = ACF5YZ*FSTORA(JS,KS)
     +                       + BCF5YZ*FSTORB(JS,KS)
     +                       + CCF5YZ*FSTORC(JS,KS)
     +                       + DCF5YZ*FDIFFD

          ENDDO

        ENDIF

      ENDIF

C     =========================================================================

C     SCALING
C     =======
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDELY*OVDELZ

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
