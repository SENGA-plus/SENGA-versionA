      SUBROUTINE D2FDXY(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     D2FDXY
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
C     EVALUATES SECOND XY-DERIVATIVE OF SPECIFIED FUNCTION
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
      INTEGER IS,JS,ISM1,JSM1
      INTEGER ISTART,IFINIS,JSTART,JFINIS
      INTEGER ICM5,ICM4,ICM3,ICM2,ICM1,ICCC,ICP1,ICP2,ICP3,ICP4,ICP5
      INTEGER JCM5,JCM4,JCM3,JCM2,JCM1,JCCC,JCP1,JCP2,JCP3,JCP4,JCP5


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      ISTART = ISTAL
      IFINIS = ISTOL
      JSTART = JSTAL
      JFINIS = JSTOL
      IF(NENDXL.EQ.NBOUND)ISTART = ISTAP5
      IF(NENDXR.EQ.NBOUND)IFINIS = ISTOM5
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

            FDIFFA = FUNCTN(ICP1,JCP1,KC) - FUNCTN(ICP1,JCM1,KC) 
     +             - FUNCTN(ICM1,JCP1,KC) + FUNCTN(ICM1,JCM1,KC) 
            FDIFFB = FUNCTN(ICP2,JCP2,KC) - FUNCTN(ICP2,JCM2,KC) 
     +             - FUNCTN(ICM2,JCP2,KC) + FUNCTN(ICM2,JCM2,KC) 
            FDIFFC = FUNCTN(ICP3,JCP3,KC) - FUNCTN(ICP3,JCM3,KC) 
     +             - FUNCTN(ICM3,JCP3,KC) + FUNCTN(ICM3,JCM3,KC) 
            FDIFFD = FUNCTN(ICP4,JCP4,KC) - FUNCTN(ICP4,JCM4,KC) 
     +             - FUNCTN(ICM4,JCP4,KC) + FUNCTN(ICM4,JCM4,KC) 
            FDIFFE = FUNCTN(ICP5,JCP5,KC) - FUNCTN(ICP5,JCM5,KC) 
     +             - FUNCTN(ICM5,JCP5,KC) + FUNCTN(ICM5,JCM5,KC) 

            FDERIV(IC,JC,KC) = ACOFXY*FDIFFA
     +                       + BCOFXY*FDIFFB
     +                       + CCOFXY*FDIFFC
     +                       + DCOFXY*FDIFFD
     +                       + ECOFXY*FDIFFE

          ENDDO

        ENDDO

      ENDDO

C     =========================================================================

C     LH END X-DIRECTION
C     ==================
      IF(NENDXL.EQ.NBOUND)THEN

C       TAKE SECOND XY-DERIVATIVE IN X-LEFT INNER HALO
C       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL

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

C           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(ISTAP1,JCP1,KC) - FUNCTN(ISTAP1,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP1,JCP2,KC) - FUNCTN(ISTAP1,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
            FDIFFB =
     +          ACOFY1*(FUNCTN(ISTAP2,JCP1,KC) - FUNCTN(ISTAP2,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP2,JCP2,KC) - FUNCTN(ISTAP2,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
            FDIFFC =
     +          ACOFY1*(FUNCTN(ISTAP3,JCP1,KC) - FUNCTN(ISTAP3,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP3,JCP2,KC) - FUNCTN(ISTAP3,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
            FDIFFD =
     +          ACOFY1*(FUNCTN(ISTAP4,JCP1,KC) - FUNCTN(ISTAP4,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP4,JCP2,KC) - FUNCTN(ISTAP4,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
            FDERIV(ISTAL,JC,KC) = ACF1XY*FDIFFA
     +                          + BCF1XY*FDIFFB
     +                          + CCF1XY*FDIFFC
     +                          + DCF1XY*FDIFFD
     
C           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(ISTAL,JCP1,KC)  - FUNCTN(ISTAL,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAL,JCP2,KC)  - FUNCTN(ISTAL,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
            FDIFFB =
     +          ACOFY1*(FUNCTN(ISTAP2,JCP1,KC) - FUNCTN(ISTAP2,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP2,JCP2,KC) - FUNCTN(ISTAP2,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
            FDIFFC =
     +          ACOFY1*(FUNCTN(ISTAP3,JCP1,KC) - FUNCTN(ISTAP3,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP3,JCP2,KC) - FUNCTN(ISTAP3,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
            FDIFFD =
     +          ACOFY1*(FUNCTN(ISTAP4,JCP1,KC) - FUNCTN(ISTAP4,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP4,JCP2,KC) - FUNCTN(ISTAP4,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
            FDERIV(ISTAP1,JC,KC) = ACF2XY*FDIFFA
     +                           + BCF2XY*FDIFFB
     +                           + CCF2XY*FDIFFC
     +                           + DCF2XY*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP3,JCP1,KC) - FUNCTN(ISTAP3,JCM1,KC) 
     +             - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC) 
            FDIFFB = FUNCTN(ISTAP4,JCP2,KC) - FUNCTN(ISTAP4,JCM2,KC) 
     +             - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC) 
            FDERIV(ISTAP2,JC,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP4,JCP1,KC) - FUNCTN(ISTAP4,JCM1,KC) 
     +             - FUNCTN(ISTAP2,JCP1,KC) + FUNCTN(ISTAP2,JCM1,KC) 
            FDIFFB = FUNCTN(ISTAP5,JCP2,KC) - FUNCTN(ISTAP5,JCM2,KC) 
     +             - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC) 
            FDIFFC = FUNCTN(ISTAP6,JCP3,KC) - FUNCTN(ISTAP6,JCM3,KC) 
     +             - FUNCTN(ISTAL,JCP3,KC)  + FUNCTN(ISTAL,JCM3,KC) 
            FDERIV(ISTAP3,JC,KC) = ACF4XY*FDIFFA
     +                           + BCF4XY*FDIFFB
     +                           + CCF4XY*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP5,JCP1,KC) - FUNCTN(ISTAP5,JCM1,KC) 
     +             - FUNCTN(ISTAP3,JCP1,KC) + FUNCTN(ISTAP3,JCM1,KC) 
            FDIFFB = FUNCTN(ISTAP6,JCP2,KC) - FUNCTN(ISTAP6,JCM2,KC) 
     +             - FUNCTN(ISTAP2,JCP2,KC) + FUNCTN(ISTAP2,JCM2,KC) 
            FDIFFC = FUNCTN(ISTAP7,JCP3,KC) - FUNCTN(ISTAP7,JCM3,KC) 
     +             - FUNCTN(ISTAP1,JCP3,KC) + FUNCTN(ISTAP1,JCM3,KC) 
            FDIFFD = FUNCTN(ISTAP8,JCP4,KC) - FUNCTN(ISTAP8,JCM4,KC) 
     +             - FUNCTN(ISTAL,JCP4,KC)  + FUNCTN(ISTAL,JCM4,KC) 
            FDERIV(ISTAP4,JC,KC) = ACF5XY*FDIFFA
     +                           + BCF5XY*FDIFFB
     +                           + CCF5XY*FDIFFC
     +                           + DCF5XY*FDIFFD
      
          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END X-DIRECTION
C     ==================
      IF(NENDXR.EQ.NBOUND)THEN

C       TAKE SECOND XY-DERIVATIVE IN X-RIGHT INNER HALO
C       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL

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

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM3,JCP1,KC) - FUNCTN(ISTOM3,JCM1,KC) 
     +             - FUNCTN(ISTOM5,JCP1,KC) + FUNCTN(ISTOM5,JCM1,KC) 
            FDIFFB = FUNCTN(ISTOM2,JCP2,KC) - FUNCTN(ISTOM2,JCM2,KC) 
     +             - FUNCTN(ISTOM6,JCP2,KC) + FUNCTN(ISTOM6,JCM2,KC) 
            FDIFFC = FUNCTN(ISTOM1,JCP3,KC) - FUNCTN(ISTOM1,JCM3,KC) 
     +             - FUNCTN(ISTOM7,JCP3,KC) + FUNCTN(ISTOM7,JCM3,KC) 
            FDIFFD = FUNCTN(ISTOL,JCP4,KC)  - FUNCTN(ISTOL,JCM4,KC) 
     +             - FUNCTN(ISTOM8,JCP4,KC) + FUNCTN(ISTOM8,JCM4,KC) 
            FDERIV(ISTOM4,JC,KC) = ACF5XY*FDIFFA
     +                           + BCF5XY*FDIFFB
     +                           + CCF5XY*FDIFFC
     +                           + DCF5XY*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM2,JCP1,KC) - FUNCTN(ISTOM2,JCM1,KC) 
     +             - FUNCTN(ISTOM4,JCP1,KC) + FUNCTN(ISTOM4,JCM1,KC) 
            FDIFFB = FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC) 
     +             - FUNCTN(ISTOM5,JCP2,KC) + FUNCTN(ISTOM5,JCM2,KC) 
            FDIFFC = FUNCTN(ISTOL,JCP3,KC)  - FUNCTN(ISTOL,JCM3,KC) 
     +             - FUNCTN(ISTOM6,JCP3,KC) + FUNCTN(ISTOM6,JCM3,KC) 
            FDERIV(ISTOM3,JC,KC) = ACF4XY*FDIFFA
     +                           + BCF4XY*FDIFFB
     +                           + CCF4XY*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC) 
     +             - FUNCTN(ISTOM3,JCP1,KC) + FUNCTN(ISTOM3,JCM1,KC) 
            FDIFFB = FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC) 
     +             - FUNCTN(ISTOM4,JCP2,KC) + FUNCTN(ISTOM4,JCM2,KC) 
            FDERIV(ISTOM2,JC,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOL,JCP1,KC)  + FUNCTN(ISTOL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOL,JCP2,KC)  + FUNCTN(ISTOL,JCM2,KC))
            FDIFFB =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM2,JCP1,KC) + FUNCTN(ISTOM2,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM2,JCP2,KC) + FUNCTN(ISTOM2,JCM2,KC))
            FDIFFC =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM3,JCP1,KC) + FUNCTN(ISTOM3,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM3,JCP2,KC) + FUNCTN(ISTOM3,JCM2,KC))
            FDIFFD =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM4,JCP1,KC) + FUNCTN(ISTOM4,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM4,JCP2,KC) + FUNCTN(ISTOM4,JCM2,KC))
            FDERIV(ISTOM1,JC,KC) = ACF2XY*FDIFFA
     +                           + BCF2XY*FDIFFB
     +                           + CCF2XY*FDIFFC
     +                           + DCF2XY*FDIFFD

C           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM1,JCP1,KC) + FUNCTN(ISTOM1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM1,JCP2,KC) + FUNCTN(ISTOM1,JCM2,KC))
            FDIFFB =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM2,JCP1,KC) + FUNCTN(ISTOM2,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM2,JCP2,KC) + FUNCTN(ISTOM2,JCM2,KC))
            FDIFFC =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM3,JCP1,KC) + FUNCTN(ISTOM3,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM3,JCP2,KC) + FUNCTN(ISTOM3,JCM2,KC))
            FDIFFD =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM4,JCP1,KC) + FUNCTN(ISTOM4,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM4,JCP2,KC) + FUNCTN(ISTOM4,JCM2,KC))
            FDERIV(ISTOL,JC,KC) = ACF1XY*FDIFFA
     +                          + BCF1XY*FDIFFB
     +                          + CCF1XY*FDIFFC
     +                          + DCF1XY*FDIFFD

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     LH END Y-DIRECTION
C     ==================
      IF(NENDYL.EQ.NBOUND)THEN

C       TAKE SECOND XY-DERIVATIVE IN Y-LEFT INNER HALO
C       EXPLICIT 4TH,4TH,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL

          ICM3 = ISTART-4
          ICM2 = ISTART-3
          ICM1 = ISTART-2
          ICCC = ISTART-1
          ICP1 = ISTART
          ICP2 = ISTART+1
          ICP3 = ISTART+2
          ICP4 = ISTART+3

          DO IC = ISTART,IFINIS

            ICM4 = ICM3
            ICM3 = ICM2
            ICM2 = ICM1
            ICM1 = ICCC
            ICCC = ICP1
            ICP1 = ICP2
            ICP2 = ICP3
            ICP3 = ICP4
            ICP4 = IC+4

C           LH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP1,KC) - FUNCTN(ICM1,JSTAP1,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP1,KC) - FUNCTN(ICM2,JSTAP1,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP2,KC) - FUNCTN(ICM1,JSTAP2,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP2,KC) - FUNCTN(ICM2,JSTAP2,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP3,KC) - FUNCTN(ICM1,JSTAP3,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP3,KC) - FUNCTN(ICM2,JSTAP3,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP4,KC) - FUNCTN(ICM1,JSTAP4,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP4,KC) - FUNCTN(ICM2,JSTAP4,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
            FDERIV(IC,JSTAL,KC) = ACF1XY*FDIFFA
     +                          + BCF1XY*FDIFFB
     +                          + CCF1XY*FDIFFC
     +                          + DCF1XY*FDIFFD
     
C           LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTAL,KC)  - FUNCTN(ICM1,JSTAL,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAL,KC)  - FUNCTN(ICM2,JSTAL,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP2,KC) - FUNCTN(ICM1,JSTAP2,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP2,KC) - FUNCTN(ICM2,JSTAP2,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP3,KC) - FUNCTN(ICM1,JSTAP3,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP3,KC) - FUNCTN(ICM2,JSTAP3,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP4,KC) - FUNCTN(ICM1,JSTAP4,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP4,KC) - FUNCTN(ICM2,JSTAP4,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
            FDERIV(IC,JSTAP1,KC) = ACF2XY*FDIFFA
     +                           + BCF2XY*FDIFFB
     +                           + CCF2XY*FDIFFC
     +                           + DCF2XY*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JSTAP3,KC) - FUNCTN(ICM1,JSTAP3,KC) 
     +             - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC) 
            FDIFFB = FUNCTN(ICP2,JSTAP4,KC) - FUNCTN(ICM2,JSTAP4,KC) 
     +             - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC) 
            FDERIV(IC,JSTAP2,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JSTAP4,KC) - FUNCTN(ICM1,JSTAP4,KC) 
     +             - FUNCTN(ICP1,JSTAP2,KC) + FUNCTN(ICM1,JSTAP2,KC) 
            FDIFFB = FUNCTN(ICP2,JSTAP5,KC) - FUNCTN(ICM2,JSTAP5,KC) 
     +             - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC) 
            FDIFFC = FUNCTN(ICP3,JSTAP6,KC) - FUNCTN(ICM3,JSTAP6,KC) 
     +             - FUNCTN(ICP3,JSTAL,KC)  + FUNCTN(ICM3,JSTAL,KC) 
            FDERIV(IC,JSTAP3,KC) = ACF4XY*FDIFFA
     +                           + BCF4XY*FDIFFB
     +                           + CCF4XY*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JSTAP5,KC) - FUNCTN(ICM1,JSTAP5,KC) 
     +             - FUNCTN(ICP1,JSTAP3,KC) + FUNCTN(ICM1,JSTAP3,KC) 
            FDIFFB = FUNCTN(ICP2,JSTAP6,KC) - FUNCTN(ICM2,JSTAP6,KC) 
     +             - FUNCTN(ICP2,JSTAP2,KC) + FUNCTN(ICM2,JSTAP2,KC) 
            FDIFFC = FUNCTN(ICP3,JSTAP7,KC) - FUNCTN(ICM3,JSTAP7,KC) 
     +             - FUNCTN(ICP3,JSTAP1,KC) + FUNCTN(ICM3,JSTAP1,KC) 
            FDIFFD = FUNCTN(ICP4,JSTAP8,KC) - FUNCTN(ICM4,JSTAP8,KC) 
     +             - FUNCTN(ICP4,JSTAL,KC)  + FUNCTN(ICM4,JSTAL,KC) 
            FDERIV(IC,JSTAP4,KC) = ACF5XY*FDIFFA
     +                           + BCF5XY*FDIFFB
     +                           + CCF5XY*FDIFFC
     +                           + DCF5XY*FDIFFD
      
          ENDDO
        ENDDO

C       LH IN X LH IN Y CORNER
C       ======================
        IF(NENDXL.EQ.NBOUND)THEN

          DO KC = KSTAL,KSTOL

C           LH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTAP1,JSTAP1,KC) - FUNCTN(ISTAP1,JSTAL,KC)
     +             - FUNCTN(ISTAL,JSTAP1,KC)  + FUNCTN(ISTAL,JSTAL,KC)
            FDIFFB = FUNCTN(ISTAP2,JSTAP2,KC) - FUNCTN(ISTAP2,JSTAL,KC)
     +             - FUNCTN(ISTAL,JSTAP2,KC)  + FUNCTN(ISTAL,JSTAL,KC)
            FDIFFC = FUNCTN(ISTAP3,JSTAP3,KC) - FUNCTN(ISTAP3,JSTAL,KC)
     +             - FUNCTN(ISTAL,JSTAP3,KC)  + FUNCTN(ISTAL,JSTAL,KC)
            FDIFFD = FUNCTN(ISTAP4,JSTAP4,KC) - FUNCTN(ISTAP4,JSTAL,KC)
     +             - FUNCTN(ISTAL,JSTAP4,KC)  + FUNCTN(ISTAL,JSTAL,KC)
            FDERIV(ISTAL,JSTAL,KC) = ACC1XY*FDIFFA
     +                             + BCC1XY*FDIFFB
     +                             + CCC1XY*FDIFFC
     +                             + DCC1XY*FDIFFD

C           LH+1 LH+1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTAL,JSTAL,KC)   - FUNCTN(ISTAL,JSTAP1,KC)
     +             - FUNCTN(ISTAP1,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAP1,KC)
            FDIFFB = FUNCTN(ISTAP2,JSTAP2,KC) - FUNCTN(ISTAP2,JSTAP1,KC)
     +             - FUNCTN(ISTAP1,JSTAP2,KC) + FUNCTN(ISTAP1,JSTAP1,KC)
            FDIFFC = FUNCTN(ISTAP3,JSTAP3,KC) - FUNCTN(ISTAP3,JSTAP1,KC)
     +             - FUNCTN(ISTAP1,JSTAP3,KC) + FUNCTN(ISTAP1,JSTAP1,KC)
            FDIFFD = FUNCTN(ISTAP4,JSTAP4,KC) - FUNCTN(ISTAP4,JSTAP1,KC)
     +             - FUNCTN(ISTAP1,JSTAP4,KC) + FUNCTN(ISTAP1,JSTAP1,KC)
            FDERIV(ISTAP1,JSTAP1,KC) = ACC2XY*FDIFFA
     +                               + BCC2XY*FDIFFB
     +                               + CCC2XY*FDIFFC
     +                               + DCC2XY*FDIFFD

C           LH LH+1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTAP1,JSTAL,KC)  - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAL,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTAP1,JSTAP2,KC) - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP2,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTAP1,JSTAP3,KC) - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP3,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTAP1,JSTAP4,KC) - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP4,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTAP2,JSTAL,KC)  - FUNCTN(ISTAP2,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAL,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTAP2,KC) - FUNCTN(ISTAP2,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP2,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTAP2,JSTAP3,KC) - FUNCTN(ISTAP2,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP3,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTAP2,JSTAP4,KC) - FUNCTN(ISTAP2,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP4,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTAP3,JSTAL,KC)  - FUNCTN(ISTAP3,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAL,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTAP3,JSTAP2,KC) - FUNCTN(ISTAP3,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP2,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTAP3,KC) - FUNCTN(ISTAP3,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP3,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTAP3,JSTAP4,KC) - FUNCTN(ISTAP3,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP4,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTAP4,JSTAL,KC)  - FUNCTN(ISTAP4,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAL,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTAP4,JSTAP2,KC) - FUNCTN(ISTAP4,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP2,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTAP4,JSTAP3,KC) - FUNCTN(ISTAP4,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP3,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTAP4,KC) - FUNCTN(ISTAP4,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAP4,KC)  + FUNCTN(ISTAL,JSTAP1,KC))
            FDERIV(ISTAL,JSTAP1,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           LH+1 LH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTAL,JSTAP1,KC)  - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAP1,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTAP1,KC) - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAP2,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTAP1,KC) - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAP3,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTAP1,KC) - FUNCTN(ISTAP1,JSTAP1,KC)
     +            - FUNCTN(ISTAP4,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTAL,JSTAP2,KC)  - FUNCTN(ISTAP1,JSTAP2,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAP1,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTAP2,KC) - FUNCTN(ISTAP1,JSTAP2,KC)
     +            - FUNCTN(ISTAP2,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTAP2,KC) - FUNCTN(ISTAP1,JSTAP2,KC)
     +            - FUNCTN(ISTAP3,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTAP2,KC) - FUNCTN(ISTAP1,JSTAP2,KC)
     +            - FUNCTN(ISTAP4,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTAL,JSTAP3,KC)  - FUNCTN(ISTAP1,JSTAP3,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAP1,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTAP3,KC) - FUNCTN(ISTAP1,JSTAP3,KC)
     +            - FUNCTN(ISTAP2,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTAP3,KC) - FUNCTN(ISTAP1,JSTAP3,KC)
     +            - FUNCTN(ISTAP3,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTAP3,KC) - FUNCTN(ISTAP1,JSTAP3,KC)
     +            - FUNCTN(ISTAP4,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTAL,JSTAP4,KC)  - FUNCTN(ISTAP1,JSTAP4,KC)
     +            - FUNCTN(ISTAL,JSTAL,KC)   + FUNCTN(ISTAP1,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTAP4,KC) - FUNCTN(ISTAP1,JSTAP4,KC)
     +            - FUNCTN(ISTAP2,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTAP4,KC) - FUNCTN(ISTAP1,JSTAP4,KC)
     +            - FUNCTN(ISTAP3,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTAP4,KC) - FUNCTN(ISTAP1,JSTAP4,KC)
     +            - FUNCTN(ISTAP4,JSTAL,KC)  + FUNCTN(ISTAP1,JSTAL,KC))
            FDERIV(ISTAP1,JSTAL,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           LH EDGE IN Y
            DO IC = ISTAP2,ISTAP4

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP1,KC) - FUNCTN(ICM1,JSTAP1,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP1,KC) - FUNCTN(ICM2,JSTAP1,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP2,KC) - FUNCTN(ICM1,JSTAP2,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP2,KC) - FUNCTN(ICM2,JSTAP2,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP3,KC) - FUNCTN(ICM1,JSTAP3,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP3,KC) - FUNCTN(ICM2,JSTAP3,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP4,KC) - FUNCTN(ICM1,JSTAP4,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP4,KC) - FUNCTN(ICM2,JSTAP4,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDERIV(IC,JSTAL,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTAL,KC)  - FUNCTN(ICM1,JSTAL,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAL,KC)  - FUNCTN(ICM2,JSTAL,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP2,KC) - FUNCTN(ICM1,JSTAP2,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP2,KC) - FUNCTN(ICM2,JSTAP2,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP3,KC) - FUNCTN(ICM1,JSTAP3,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP3,KC) - FUNCTN(ICM2,JSTAP3,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP4,KC) - FUNCTN(ICM1,JSTAP4,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP4,KC) - FUNCTN(ICM2,JSTAP4,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDERIV(IC,JSTAP1,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

            ENDDO

C           LH EDGE IN X
            DO JC = JSTAP2,JSTAP4

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTAP1,JCP1,KC) - FUNCTN(ISTAP1,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP1,JCP2,KC) - FUNCTN(ISTAP1,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTAP2,JCP1,KC) - FUNCTN(ISTAP2,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP2,JCP2,KC) - FUNCTN(ISTAP2,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTAP3,JCP1,KC) - FUNCTN(ISTAP3,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP3,JCP2,KC) - FUNCTN(ISTAP3,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTAP4,JCP1,KC) - FUNCTN(ISTAP4,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP4,JCP2,KC) - FUNCTN(ISTAP4,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDERIV(ISTAL,JC,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTAL,JCP1,KC)  - FUNCTN(ISTAL,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAL,JCP2,KC)  - FUNCTN(ISTAL,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTAP2,JCP1,KC) - FUNCTN(ISTAP2,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP2,JCP2,KC) - FUNCTN(ISTAP2,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTAP3,JCP1,KC) - FUNCTN(ISTAP3,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP3,JCP2,KC) - FUNCTN(ISTAP3,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTAP4,JCP1,KC) - FUNCTN(ISTAP4,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP4,JCP2,KC) - FUNCTN(ISTAP4,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDERIV(ISTAP1,JC,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            JS = 0
            DO JC = JSTAP2,JSTAP4

              JS = JS+1
              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

              IS = 0
              DO IC = ISTAP2,ISTAP4

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JCP1,KC) - FUNCTN(ICP1,JCM1,KC)
     +                 - FUNCTN(ICM1,JCP1,KC) + FUNCTN(ICM1,JCM1,KC)
                FDIFFB = FUNCTN(ICP2,JCP2,KC) - FUNCTN(ICP2,JCM2,KC)
     +                 - FUNCTN(ICM2,JCP2,KC) + FUNCTN(ICM2,JCM2,KC)
                FDERIV(IC,JC,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB
                FSTORA(IS,JS) = FDIFFA
                FSTORB(IS,JS) = FDIFFB
 
              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            JS = 1
            DO JC = JSTAP3,JSTAP4

              JSM1 = JS
              JS = JS+1
              JCM3 = JC-3
              JCP3 = JC+3

              IS = 1
              DO IC = ISTAP3,ISTAP4

                ISM1 = IS
                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JCP3,KC) - FUNCTN(ICP3,JCM3,KC)
     +                 - FUNCTN(ICM3,JCP3,KC) + FUNCTN(ICM3,JCM3,KC)
                FDERIV(IC,JC,KC) = ACF4XY*FSTORA(IS,JS)
     +                           + BCF4XY*FSTORB(IS,JS)
     +                           + CCF4XY*FDIFFC
                FSTORC(ISM1,JSM1) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            JS = 3
            IS = 3
            JSM1 = 2
            ISM1 = 2
            JC = JSTAP4
            IC = ISTAP4
            JCM4 = JC-4
            JCP4 = JC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JCP4,KC) - FUNCTN(ICP4,JCM4,KC)
     +             - FUNCTN(ICM4,JCP4,KC) + FUNCTN(ICM4,JCM4,KC)
            FDERIV(IC,JC,KC) = ACF5XY*FSTORA(IS,JS)
     +                       + BCF5XY*FSTORB(IS,JS)
     +                       + CCF5XY*FSTORC(ISM1,JSM1)
     +                       + DCF5XY*FDIFFD

          ENDDO

        ENDIF

C       RH IN X LH IN Y CORNER
C       ======================
        IF(NENDXR.EQ.NBOUND)THEN

          DO KC = KSTAL,KSTOL

C           RH LH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTOL,JSTAP1,KC)  - FUNCTN(ISTOL,JSTAL,KC)
     +             - FUNCTN(ISTOM1,JSTAP1,KC) + FUNCTN(ISTOM1,JSTAL,KC)
            FDIFFB = FUNCTN(ISTOL,JSTAP2,KC)  - FUNCTN(ISTOL,JSTAL,KC)
     +             - FUNCTN(ISTOM2,JSTAP2,KC) + FUNCTN(ISTOM2,JSTAL,KC)
            FDIFFC = FUNCTN(ISTOL,JSTAP3,KC)  - FUNCTN(ISTOL,JSTAL,KC)
     +             - FUNCTN(ISTOM3,JSTAP3,KC) + FUNCTN(ISTOM3,JSTAL,KC)
            FDIFFD = FUNCTN(ISTOL,JSTAP4,KC)  - FUNCTN(ISTOL,JSTAL,KC)
     +             - FUNCTN(ISTOM4,JSTAP4,KC) + FUNCTN(ISTOM4,JSTAL,KC)
            FDERIV(ISTOL,JSTAL,KC) = ACC1XY*FDIFFA
     +                             + BCC1XY*FDIFFB
     +                             + CCC1XY*FDIFFC
     +                             + DCC1XY*FDIFFD

C           RH-1 LH+1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTOM1,JSTAL,KC)  - FUNCTN(ISTOM1,JSTAP1,KC)
     +             - FUNCTN(ISTOL,JSTAL,KC)   + FUNCTN(ISTOL,JSTAP1,KC)
            FDIFFB = FUNCTN(ISTOM1,JSTAP2,KC) - FUNCTN(ISTOM1,JSTAP1,KC)
     +             - FUNCTN(ISTOM2,JSTAP2,KC) + FUNCTN(ISTOM2,JSTAP1,KC)
            FDIFFC = FUNCTN(ISTOM1,JSTAP3,KC) - FUNCTN(ISTOM1,JSTAP1,KC)
     +             - FUNCTN(ISTOM3,JSTAP3,KC) + FUNCTN(ISTOM3,JSTAP1,KC)
            FDIFFD = FUNCTN(ISTOM1,JSTAP4,KC) - FUNCTN(ISTOM1,JSTAP1,KC)
     +             - FUNCTN(ISTOM4,JSTAP4,KC) + FUNCTN(ISTOM4,JSTAP1,KC)
            FDERIV(ISTOM1,JSTAP1,KC) = ACC2XY*FDIFFA
     +                               + BCC2XY*FDIFFB
     +                               + CCC2XY*FDIFFC
     +                               + DCC2XY*FDIFFD

C           RH LH+1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTOL,JSTAL,KC)   - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM1,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTOL,JSTAP2,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAP2,KC) + FUNCTN(ISTOM1,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTOL,JSTAP3,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAP3,KC) + FUNCTN(ISTOM1,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTOL,JSTAP4,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAP4,KC) + FUNCTN(ISTOM1,JSTAP1,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTOL,JSTAL,KC)   - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM2,JSTAL,KC)  + FUNCTN(ISTOM2,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTOL,JSTAP2,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM2,JSTAP2,KC) + FUNCTN(ISTOM2,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTOL,JSTAP3,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM2,JSTAP3,KC) + FUNCTN(ISTOM2,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTOL,JSTAP4,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM2,JSTAP4,KC) + FUNCTN(ISTOM2,JSTAP1,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTOL,JSTAL,KC)   - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM3,JSTAL,KC)  + FUNCTN(ISTOM3,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTOL,JSTAP2,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM3,JSTAP2,KC) + FUNCTN(ISTOM3,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTOL,JSTAP3,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM3,JSTAP3,KC) + FUNCTN(ISTOM3,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTOL,JSTAP4,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM3,JSTAP4,KC) + FUNCTN(ISTOM3,JSTAP1,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTOL,JSTAL,KC)   - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM4,JSTAL,KC)  + FUNCTN(ISTOM4,JSTAP1,KC))
     +    + BCF2XY*(FUNCTN(ISTOL,JSTAP2,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM4,JSTAP2,KC) + FUNCTN(ISTOM4,JSTAP1,KC))
     +    + CCF2XY*(FUNCTN(ISTOL,JSTAP3,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM4,JSTAP3,KC) + FUNCTN(ISTOM4,JSTAP1,KC))
     +    + DCF2XY*(FUNCTN(ISTOL,JSTAP4,KC)  - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM4,JSTAP4,KC) + FUNCTN(ISTOM4,JSTAP1,KC))
            FDERIV(ISTOL,JSTAP1,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           RH-1 LH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTOM1,JSTAP1,KC) - FUNCTN(ISTOL,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOL,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM1,JSTAP1,KC) - FUNCTN(ISTOM2,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM2,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM1,JSTAP1,KC) - FUNCTN(ISTOM3,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM3,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM1,JSTAP1,KC) - FUNCTN(ISTOM4,JSTAP1,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM4,JSTAL,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTOM1,JSTAP2,KC) - FUNCTN(ISTOL,JSTAP2,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOL,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM1,JSTAP2,KC) - FUNCTN(ISTOM2,JSTAP2,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM2,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM1,JSTAP2,KC) - FUNCTN(ISTOM3,JSTAP2,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM3,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM1,JSTAP2,KC) - FUNCTN(ISTOM4,JSTAP2,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM4,JSTAL,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTOM1,JSTAP3,KC) - FUNCTN(ISTOL,JSTAP3,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOL,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM1,JSTAP3,KC) - FUNCTN(ISTOM2,JSTAP3,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM2,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM1,JSTAP3,KC) - FUNCTN(ISTOM3,JSTAP3,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM3,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM1,JSTAP3,KC) - FUNCTN(ISTOM4,JSTAP3,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM4,JSTAL,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTOM1,JSTAP4,KC) - FUNCTN(ISTOL,JSTAP4,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOL,JSTAL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM1,JSTAP4,KC) - FUNCTN(ISTOM2,JSTAP4,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM2,JSTAL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM1,JSTAP4,KC) - FUNCTN(ISTOM3,JSTAP4,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM3,JSTAL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM1,JSTAP4,KC) - FUNCTN(ISTOM4,JSTAP4,KC)
     +            - FUNCTN(ISTOM1,JSTAL,KC)  + FUNCTN(ISTOM4,JSTAL,KC))
            FDERIV(ISTOM1,JSTAL,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           LH EDGE IN Y
            DO IC = ISTOM4,ISTOM2

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP1,KC) - FUNCTN(ICM1,JSTAP1,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP1,KC) - FUNCTN(ICM2,JSTAP1,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP2,KC) - FUNCTN(ICM1,JSTAP2,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP2,KC) - FUNCTN(ICM2,JSTAP2,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP3,KC) - FUNCTN(ICM1,JSTAP3,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP3,KC) - FUNCTN(ICM2,JSTAP3,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP4,KC) - FUNCTN(ICM1,JSTAP4,KC)
     +                - FUNCTN(ICP1,JSTAL,KC)  + FUNCTN(ICM1,JSTAL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP4,KC) - FUNCTN(ICM2,JSTAP4,KC)
     +                - FUNCTN(ICP2,JSTAL,KC)  + FUNCTN(ICM2,JSTAL,KC))
              FDERIV(IC,JSTAL,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTAL,KC)  - FUNCTN(ICM1,JSTAL,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAL,KC)  - FUNCTN(ICM2,JSTAL,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP2,KC) - FUNCTN(ICM1,JSTAP2,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP2,KC) - FUNCTN(ICM2,JSTAP2,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP3,KC) - FUNCTN(ICM1,JSTAP3,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP3,KC) - FUNCTN(ICM2,JSTAP3,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTAP4,KC) - FUNCTN(ICM1,JSTAP4,KC)
     +                - FUNCTN(ICP1,JSTAP1,KC) + FUNCTN(ICM1,JSTAP1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTAP4,KC) - FUNCTN(ICM2,JSTAP4,KC)
     +                - FUNCTN(ICP2,JSTAP1,KC) + FUNCTN(ICM2,JSTAP1,KC))
              FDERIV(IC,JSTAP1,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

            ENDDO

C           RH EDGE IN X
            DO JC = JSTAP2,JSTAP4

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM1,JCP1,KC) + FUNCTN(ISTOM1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM1,JCP2,KC) + FUNCTN(ISTOM1,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM2,JCP1,KC) + FUNCTN(ISTOM2,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM2,JCP2,KC) + FUNCTN(ISTOM2,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM3,JCP1,KC) + FUNCTN(ISTOM3,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM3,JCP2,KC) + FUNCTN(ISTOM3,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM4,JCP1,KC) + FUNCTN(ISTOM4,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM4,JCP2,KC) + FUNCTN(ISTOM4,JCM2,KC))
              FDERIV(ISTOL,JC,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD
    
C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOL,JCP1,KC)  + FUNCTN(ISTOL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOL,JCP2,KC)  + FUNCTN(ISTOL,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM2,JCP1,KC) + FUNCTN(ISTOM2,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM2,JCP2,KC) + FUNCTN(ISTOM2,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM3,JCP1,KC) + FUNCTN(ISTOM3,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM3,JCP2,KC) + FUNCTN(ISTOM3,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM4,JCP1,KC) + FUNCTN(ISTOM4,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM4,JCP2,KC) + FUNCTN(ISTOM4,JCM2,KC))
              FDERIV(ISTOM1,JC,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            JS = 0
            DO JC = JSTAP2,JSTAP4

              JS = JS+1
              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

              IS = 0
              DO IC = ISTOM4,ISTOM2

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JCP1,KC) - FUNCTN(ICP1,JCM1,KC)
     +                 - FUNCTN(ICM1,JCP1,KC) + FUNCTN(ICM1,JCM1,KC)
                FDIFFB = FUNCTN(ICP2,JCP2,KC) - FUNCTN(ICP2,JCM2,KC)
     +                 - FUNCTN(ICM2,JCP2,KC) + FUNCTN(ICM2,JCM2,KC)
                FDERIV(IC,JC,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB
                FSTORA(IS,JS) = FDIFFA
                FSTORB(IS,JS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            JS = 1
            DO JC = JSTAP3,JSTAP4

              JSM1 = JS
              JS = JS+1
              JCM3 = JC-3
              JCP3 = JC+3

              IS = 0
              DO IC = ISTOM4,ISTOM3

                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JCP3,KC) - FUNCTN(ICP3,JCM3,KC)
     +                 - FUNCTN(ICM3,JCP3,KC) + FUNCTN(ICM3,JCM3,KC)
                FDERIV(IC,JC,KC) = ACF4XY*FSTORA(IS,JS)
     +                           + BCF4XY*FSTORB(IS,JS)
     +                           + CCF4XY*FDIFFC
                FSTORC(IS,JSM1) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            JS = 3
            IS = 1
            JSM1 = 2
            JC = JSTAP4
            IC = ISTOM4
            JCM4 = JC-4
            JCP4 = JC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JCP4,KC) - FUNCTN(ICP4,JCM4,KC)
     +             - FUNCTN(ICM4,JCP4,KC) + FUNCTN(ICM4,JCM4,KC)
            FDERIV(IC,JC,KC) = ACF5XY*FSTORA(IS,JS)
     +                       + BCF5XY*FSTORB(IS,JS)
     +                       + CCF5XY*FSTORC(IS,JSM1)
     +                       + DCF5XY*FDIFFD

          ENDDO

        ENDIF

      ENDIF 

C     =========================================================================

C     RH END Y-DIRECTION
C     ==================
      IF(NENDYR.EQ.NBOUND)THEN

C       TAKE SECOND XY-DERIVATIVE IN Y-RIGHT INNER HALO
C       EXPLICIT 2ND,2ND,4TH,6TH,8TH COMPATIBLE ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL

          ICM3 = ISTART-4
          ICM2 = ISTART-3
          ICM1 = ISTART-2
          ICCC = ISTART-1
          ICP1 = ISTART
          ICP2 = ISTART+1
          ICP3 = ISTART+2
          ICP4 = ISTART+3

          DO IC = ISTART,IFINIS

            ICM4 = ICM3
            ICM3 = ICM2
            ICM2 = ICM1
            ICM1 = ICCC
            ICCC = ICP1
            ICP1 = ICP2
            ICP2 = ICP3
            ICP3 = ICP4
            ICP4 = IC+4

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JSTOM3,KC) - FUNCTN(ICM1,JSTOM3,KC) 
     +             - FUNCTN(ICP1,JSTOM5,KC) + FUNCTN(ICM1,JSTOM5,KC) 
            FDIFFB = FUNCTN(ICP2,JSTOM2,KC) - FUNCTN(ICM2,JSTOM2,KC) 
     +             - FUNCTN(ICP2,JSTOM6,KC) + FUNCTN(ICM2,JSTOM6,KC) 
            FDIFFC = FUNCTN(ICP3,JSTOM1,KC) - FUNCTN(ICM3,JSTOM1,KC) 
     +             - FUNCTN(ICP3,JSTOM7,KC) + FUNCTN(ICM3,JSTOM7,KC) 
            FDIFFD = FUNCTN(ICP4,JSTOL,KC)  - FUNCTN(ICM4,JSTOL,KC) 
     +             - FUNCTN(ICP4,JSTOM8,KC) + FUNCTN(ICM4,JSTOM8,KC) 
            FDERIV(IC,JSTOM4,KC) = ACF5XY*FDIFFA
     +                           + BCF5XY*FDIFFB
     +                           + CCF5XY*FDIFFC
     +                           + DCF5XY*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JSTOM2,KC) - FUNCTN(ICM1,JSTOM2,KC) 
     +             - FUNCTN(ICP1,JSTOM4,KC) + FUNCTN(ICM1,JSTOM4,KC) 
            FDIFFB = FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC) 
     +             - FUNCTN(ICP2,JSTOM5,KC) + FUNCTN(ICM2,JSTOM5,KC) 
            FDIFFC = FUNCTN(ICP3,JSTOL,KC)  - FUNCTN(ICM3,JSTOL,KC) 
     +             - FUNCTN(ICP3,JSTOM6,KC) + FUNCTN(ICM3,JSTOM6,KC) 
            FDERIV(IC,JSTOM3,KC) = ACF4XY*FDIFFA
     +                           + BCF4XY*FDIFFB
     +                           + CCF4XY*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC) 
     +             - FUNCTN(ICP1,JSTOM3,KC) + FUNCTN(ICM1,JSTOM3,KC) 
            FDIFFB = FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC) 
     +             - FUNCTN(ICP2,JSTOM4,KC) + FUNCTN(ICM2,JSTOM4,KC) 
            FDERIV(IC,JSTOM2,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOL,KC)  + FUNCTN(ICM1,JSTOL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOL,KC)  + FUNCTN(ICM2,JSTOL,KC))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM2,KC) + FUNCTN(ICM1,JSTOM2,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM2,KC) + FUNCTN(ICM2,JSTOM2,KC))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM3,KC) + FUNCTN(ICM1,JSTOM3,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM3,KC) + FUNCTN(ICM2,JSTOM3,KC))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM4,KC) + FUNCTN(ICM1,JSTOM4,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM4,KC) + FUNCTN(ICM2,JSTOM4,KC))
            FDERIV(IC,JSTOM1,KC) = ACF2XY*FDIFFA
     +                           + BCF2XY*FDIFFB
     +                           + CCF2XY*FDIFFC
     +                           + DCF2XY*FDIFFD

C           RH POINT: 4TH ORDER ONE-SIDED/CENTRED
            FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM1,KC) + FUNCTN(ICM1,JSTOM1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM1,KC) + FUNCTN(ICM2,JSTOM1,KC))
            FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM2,KC) + FUNCTN(ICM1,JSTOM2,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM2,KC) + FUNCTN(ICM2,JSTOM2,KC))
            FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM3,KC) + FUNCTN(ICM1,JSTOM3,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM3,KC) + FUNCTN(ICM2,JSTOM3,KC))
            FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM4,KC) + FUNCTN(ICM1,JSTOM4,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM4,KC) + FUNCTN(ICM2,JSTOM4,KC))
            FDERIV(IC,JSTOL,KC) = ACF1XY*FDIFFA
     +                          + BCF1XY*FDIFFB
     +                          + CCF1XY*FDIFFC
     +                          + DCF1XY*FDIFFD

          ENDDO
        ENDDO

C       LH IN X RH IN Y CORNER
C       ======================
        IF(NENDXL.EQ.NBOUND)THEN

          DO KC = KSTAL,KSTOL

C           LH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTAP1,JSTOL,KC) - FUNCTN(ISTAP1,JSTOM1,KC)
     +             - FUNCTN(ISTAL,JSTOL,KC)  + FUNCTN(ISTAL,JSTOM1,KC) 
            FDIFFB = FUNCTN(ISTAP2,JSTOL,KC) - FUNCTN(ISTAP2,JSTOM2,KC)
     +             - FUNCTN(ISTAL,JSTOL,KC)  + FUNCTN(ISTAL,JSTOM2,KC) 
            FDIFFC = FUNCTN(ISTAP3,JSTOL,KC) - FUNCTN(ISTAP3,JSTOM3,KC)
     +             - FUNCTN(ISTAL,JSTOL,KC)  + FUNCTN(ISTAL,JSTOM3,KC) 
            FDIFFD = FUNCTN(ISTAP4,JSTOL,KC) - FUNCTN(ISTAP4,JSTOM4,KC)
     +             - FUNCTN(ISTAL,JSTOL,KC)  + FUNCTN(ISTAL,JSTOM4,KC) 
            FDERIV(ISTAL,JSTOL,KC) = ACC1XY*FDIFFA
     +                             + BCC1XY*FDIFFB
     +                             + CCC1XY*FDIFFC
     +                             + DCC1XY*FDIFFD

C           LH+1 RH-1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTAL,JSTOM1,KC)  - FUNCTN(ISTAL,JSTOL,KC)
     +             - FUNCTN(ISTAP1,JSTOM1,KC) + FUNCTN(ISTAP1,JSTOL,KC)
            FDIFFB = FUNCTN(ISTAP2,JSTOM1,KC) - FUNCTN(ISTAP2,JSTOM2,KC)
     +             - FUNCTN(ISTAP1,JSTOM1,KC) + FUNCTN(ISTAP1,JSTOM2,KC)
            FDIFFC = FUNCTN(ISTAP3,JSTOM1,KC) - FUNCTN(ISTAP3,JSTOM3,KC)
     +             - FUNCTN(ISTAP1,JSTOM1,KC) + FUNCTN(ISTAP1,JSTOM3,KC)
            FDIFFD = FUNCTN(ISTAP4,JSTOM1,KC) - FUNCTN(ISTAP4,JSTOM4,KC)
     +             - FUNCTN(ISTAP1,JSTOM1,KC) + FUNCTN(ISTAP1,JSTOM4,KC)
            FDERIV(ISTAP1,JSTOM1,KC) = ACC2XY*FDIFFA
     +                               + BCC2XY*FDIFFB
     +                               + CCC2XY*FDIFFC
     +                               + DCC2XY*FDIFFD

C           LH RH-1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTAP1,JSTOM1,KC) - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP1,JSTOM1,KC) - FUNCTN(ISTAP1,JSTOM2,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM2,KC))
     +    + CCF2XY*(FUNCTN(ISTAP1,JSTOM1,KC) - FUNCTN(ISTAP1,JSTOM3,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM3,KC))
     +    + DCF2XY*(FUNCTN(ISTAP1,JSTOM1,KC) - FUNCTN(ISTAP1,JSTOM4,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM4,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTAP2,JSTOM1,KC) - FUNCTN(ISTAP2,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTOM1,KC) - FUNCTN(ISTAP2,JSTOM2,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM2,KC))
     +    + CCF2XY*(FUNCTN(ISTAP2,JSTOM1,KC) - FUNCTN(ISTAP2,JSTOM3,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM3,KC))
     +    + DCF2XY*(FUNCTN(ISTAP2,JSTOM1,KC) - FUNCTN(ISTAP2,JSTOM4,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM4,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTAP3,JSTOM1,KC) - FUNCTN(ISTAP3,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP3,JSTOM1,KC) - FUNCTN(ISTAP3,JSTOM2,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM2,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTOM1,KC) - FUNCTN(ISTAP3,JSTOM3,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM3,KC))
     +    + DCF2XY*(FUNCTN(ISTAP3,JSTOM1,KC) - FUNCTN(ISTAP3,JSTOM4,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM4,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTAP4,JSTOM1,KC) - FUNCTN(ISTAP4,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTAP4,JSTOM1,KC) - FUNCTN(ISTAP4,JSTOM2,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM2,KC))
     +    + CCF2XY*(FUNCTN(ISTAP4,JSTOM1,KC) - FUNCTN(ISTAP4,JSTOM3,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM3,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTOM1,KC) - FUNCTN(ISTAP4,JSTOM4,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAL,JSTOM4,KC))
            FDERIV(ISTAL,JSTOM1,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           LH+1 RH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTAL,JSTOL,KC)   - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM1,KC)  + FUNCTN(ISTAP1,JSTOM1,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP2,JSTOM1,KC) + FUNCTN(ISTAP1,JSTOM1,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP3,JSTOM1,KC) + FUNCTN(ISTAP1,JSTOM1,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP4,JSTOM1,KC) + FUNCTN(ISTAP1,JSTOM1,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTAL,JSTOL,KC)   - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM2,KC)  + FUNCTN(ISTAP1,JSTOM2,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP2,JSTOM2,KC) + FUNCTN(ISTAP1,JSTOM2,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP3,JSTOM2,KC) + FUNCTN(ISTAP1,JSTOM2,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP4,JSTOM2,KC) + FUNCTN(ISTAP1,JSTOM2,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTAL,JSTOL,KC)   - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM3,KC)  + FUNCTN(ISTAP1,JSTOM3,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP2,JSTOM3,KC) + FUNCTN(ISTAP1,JSTOM3,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP3,JSTOM3,KC) + FUNCTN(ISTAP1,JSTOM3,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP4,JSTOM3,KC) + FUNCTN(ISTAP1,JSTOM3,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTAL,JSTOL,KC)   - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAL,JSTOM4,KC)  + FUNCTN(ISTAP1,JSTOM4,KC))
     +    + BCF2XY*(FUNCTN(ISTAP2,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP2,JSTOM4,KC) + FUNCTN(ISTAP1,JSTOM4,KC))
     +    + CCF2XY*(FUNCTN(ISTAP3,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP3,JSTOM4,KC) + FUNCTN(ISTAP1,JSTOM4,KC))
     +    + DCF2XY*(FUNCTN(ISTAP4,JSTOL,KC)  - FUNCTN(ISTAP1,JSTOL,KC)
     +            - FUNCTN(ISTAP4,JSTOM4,KC) + FUNCTN(ISTAP1,JSTOM4,KC))
            FDERIV(ISTAP1,JSTOL,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           RH EDGE IN Y
            DO IC = ISTAP2,ISTAP4

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOL,KC)  + FUNCTN(ICM1,JSTOL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOL,KC)  + FUNCTN(ICM2,JSTOL,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM2,KC) + FUNCTN(ICM1,JSTOM2,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM2,KC) + FUNCTN(ICM2,JSTOM2,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM3,KC) + FUNCTN(ICM1,JSTOM3,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM3,KC) + FUNCTN(ICM2,JSTOM3,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM4,KC) + FUNCTN(ICM1,JSTOM4,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM4,KC) + FUNCTN(ICM2,JSTOM4,KC))
              FDERIV(IC,JSTOM1,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM1,KC) + FUNCTN(ICM1,JSTOM1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM1,KC) + FUNCTN(ICM2,JSTOM1,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM2,KC) + FUNCTN(ICM1,JSTOM2,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM2,KC) + FUNCTN(ICM2,JSTOM2,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM3,KC) + FUNCTN(ICM1,JSTOM3,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM3,KC) + FUNCTN(ICM2,JSTOM3,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM4,KC) + FUNCTN(ICM1,JSTOM4,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM4,KC) + FUNCTN(ICM2,JSTOM4,KC))
              FDERIV(IC,JSTOL,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD

            ENDDO

C           LH EDGE IN X
            DO JC = JSTOM4,JSTOM2

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             LH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTAP1,JCP1,KC) - FUNCTN(ISTAP1,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP1,JCP2,KC) - FUNCTN(ISTAP1,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTAP2,JCP1,KC) - FUNCTN(ISTAP2,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP2,JCP2,KC) - FUNCTN(ISTAP2,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTAP3,JCP1,KC) - FUNCTN(ISTAP3,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP3,JCP2,KC) - FUNCTN(ISTAP3,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTAP4,JCP1,KC) - FUNCTN(ISTAP4,JCM1,KC)
     +                - FUNCTN(ISTAL,JCP1,KC)  + FUNCTN(ISTAL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP4,JCP2,KC) - FUNCTN(ISTAP4,JCM2,KC)
     +                - FUNCTN(ISTAL,JCP2,KC)  + FUNCTN(ISTAL,JCM2,KC))
              FDERIV(ISTAL,JC,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD
     
C             LH POINT PLUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTAL,JCP1,KC)  - FUNCTN(ISTAL,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAL,JCP2,KC)  - FUNCTN(ISTAL,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTAP2,JCP1,KC) - FUNCTN(ISTAP2,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP2,JCP2,KC) - FUNCTN(ISTAP2,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTAP3,JCP1,KC) - FUNCTN(ISTAP3,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP3,JCP2,KC) - FUNCTN(ISTAP3,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTAP4,JCP1,KC) - FUNCTN(ISTAP4,JCM1,KC)
     +                - FUNCTN(ISTAP1,JCP1,KC) + FUNCTN(ISTAP1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTAP4,JCP2,KC) - FUNCTN(ISTAP4,JCM2,KC)
     +                - FUNCTN(ISTAP1,JCP2,KC) + FUNCTN(ISTAP1,JCM2,KC))
              FDERIV(ISTAP1,JC,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            JS = 0
            DO JC = JSTOM4,JSTOM2

              JS = JS+1
              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

              IS = 0
              DO IC = ISTAP2,ISTAP4

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JCP1,KC) - FUNCTN(ICP1,JCM1,KC)
     +                 - FUNCTN(ICM1,JCP1,KC) + FUNCTN(ICM1,JCM1,KC)
                FDIFFB = FUNCTN(ICP2,JCP2,KC) - FUNCTN(ICP2,JCM2,KC)
     +                 - FUNCTN(ICM2,JCP2,KC) + FUNCTN(ICM2,JCM2,KC)
                FDERIV(IC,JC,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB
                FSTORA(IS,JS) = FDIFFA
                FSTORB(IS,JS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            JS = 0
            DO JC = JSTOM4,JSTOM3

              JS = JS+1
              JCM3 = JC-3
              JCP3 = JC+3

              IS = 1
              DO IC = ISTAP3,ISTAP4

                ISM1 = IS
                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JCP3,KC) - FUNCTN(ICP3,JCM3,KC)
     +                 - FUNCTN(ICM3,JCP3,KC) + FUNCTN(ICM3,JCM3,KC)
                FDERIV(IC,JC,KC) = ACF4XY*FSTORA(IS,JS)
     +                           + BCF4XY*FSTORB(IS,JS)
     +                           + CCF4XY*FDIFFC
                FSTORC(ISM1,JS) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            JS = 1
            IS = 3
            ISM1 = 2
            JC = JSTOM4
            IC = ISTAP4
            JCM4 = JC-4
            JCP4 = JC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JCP4,KC) - FUNCTN(ICP4,JCM4,KC)
     +             - FUNCTN(ICM4,JCP4,KC) + FUNCTN(ICM4,JCM4,KC)
            FDERIV(IC,JC,KC) = ACF5XY*FSTORA(IS,JS)
     +                       + BCF5XY*FSTORB(IS,JS)
     +                       + CCF5XY*FSTORC(ISM1,JS)
     +                       + DCF5XY*FDIFFD

          ENDDO

        ENDIF

C       RH IN X RH IN Y CORNER
C       ======================
        IF(NENDXR.EQ.NBOUND)THEN

          DO KC = KSTAL,KSTOL

C           RH RH CORNER POINT: 4TH ORDER ONE SIDED/ONE SIDED
            FDIFFA = FUNCTN(ISTOM1,JSTOM1,KC) - FUNCTN(ISTOM1,JSTOL,KC)
     +             - FUNCTN(ISTOL,JSTOM1,KC)  + FUNCTN(ISTOL,JSTOL,KC)
            FDIFFB = FUNCTN(ISTOM2,JSTOM2,KC) - FUNCTN(ISTOM2,JSTOL,KC)
     +             - FUNCTN(ISTOL,JSTOM2,KC)  + FUNCTN(ISTOL,JSTOL,KC)
            FDIFFC = FUNCTN(ISTOM3,JSTOM3,KC) - FUNCTN(ISTOM3,JSTOL,KC)
     +             - FUNCTN(ISTOL,JSTOM3,KC)  + FUNCTN(ISTOL,JSTOL,KC)
            FDIFFD = FUNCTN(ISTOM4,JSTOM4,KC) - FUNCTN(ISTOM4,JSTOL,KC)
     +             - FUNCTN(ISTOL,JSTOM4,KC)  + FUNCTN(ISTOL,JSTOL,KC)
            FDERIV(ISTOL,JSTOL,KC) = ACC1XY*FDIFFA
     +                             + BCC1XY*FDIFFB
     +                             + CCC1XY*FDIFFC
     +                             + DCC1XY*FDIFFD

C           RH-1 RH-1 CORNER POINT: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTOL,JSTOL,KC)   - FUNCTN(ISTOL,JSTOM1,KC)
     +             - FUNCTN(ISTOM1,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOM1,KC)
            FDIFFB = FUNCTN(ISTOM2,JSTOM2,KC) - FUNCTN(ISTOM2,JSTOM1,KC)
     +             - FUNCTN(ISTOM1,JSTOM2,KC) + FUNCTN(ISTOM1,JSTOM1,KC)
            FDIFFC = FUNCTN(ISTOM3,JSTOM3,KC) - FUNCTN(ISTOM3,JSTOM1,KC)
     +             - FUNCTN(ISTOM1,JSTOM3,KC) + FUNCTN(ISTOM1,JSTOM1,KC)
            FDIFFD = FUNCTN(ISTOM4,JSTOM4,KC) - FUNCTN(ISTOM4,JSTOM1,KC)
     +             - FUNCTN(ISTOM1,JSTOM4,KC) + FUNCTN(ISTOM1,JSTOM1,KC)
            FDERIV(ISTOM1,JSTOM1,KC) = ACC2XY*FDIFFA
     +                               + BCC2XY*FDIFFB
     +                               + CCC2XY*FDIFFC
     +                               + DCC2XY*FDIFFD

C           RH RH-1 EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTOM1,JSTOL,KC)  - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOL,JSTOM1,KC))
     +    + BCF2XY*(FUNCTN(ISTOM1,JSTOM2,KC) - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM2,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + CCF2XY*(FUNCTN(ISTOM1,JSTOM3,KC) - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM3,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + DCF2XY*(FUNCTN(ISTOM1,JSTOM4,KC) - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM4,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTOM2,JSTOL,KC)  - FUNCTN(ISTOM2,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOL,JSTOM1,KC))
     +    + BCF2XY*(FUNCTN(ISTOM2,JSTOM2,KC) - FUNCTN(ISTOM2,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM2,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + CCF2XY*(FUNCTN(ISTOM2,JSTOM3,KC) - FUNCTN(ISTOM2,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM3,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + DCF2XY*(FUNCTN(ISTOM2,JSTOM4,KC) - FUNCTN(ISTOM2,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM4,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTOM3,JSTOL,KC)  - FUNCTN(ISTOM3,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOL,JSTOM1,KC))
     +    + BCF2XY*(FUNCTN(ISTOM3,JSTOM2,KC) - FUNCTN(ISTOM3,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM2,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + CCF2XY*(FUNCTN(ISTOM3,JSTOM3,KC) - FUNCTN(ISTOM3,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM3,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + DCF2XY*(FUNCTN(ISTOM3,JSTOM4,KC) - FUNCTN(ISTOM3,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM4,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTOM4,JSTOL,KC)  - FUNCTN(ISTOM4,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOL,JSTOM1,KC))
     +    + BCF2XY*(FUNCTN(ISTOM4,JSTOM2,KC) - FUNCTN(ISTOM4,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM2,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + CCF2XY*(FUNCTN(ISTOM4,JSTOM3,KC) - FUNCTN(ISTOM4,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM3,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
     +    + DCF2XY*(FUNCTN(ISTOM4,JSTOM4,KC) - FUNCTN(ISTOM4,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOM4,KC)  + FUNCTN(ISTOL,JSTOM1,KC))
            FDERIV(ISTOL,JSTOM1,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           RH+1 RH EDGE POINT: 4TH ORDER MIXED
            FDIFFA =
     +      ACF2XY*(FUNCTN(ISTOL,JSTOM1,KC)  - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOM1,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM2,JSTOM1,KC) - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOM2,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM3,JSTOM1,KC) - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOM3,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM4,JSTOM1,KC) - FUNCTN(ISTOM1,JSTOM1,KC)
     +            - FUNCTN(ISTOM4,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
            FDIFFB =
     +      ACF2XY*(FUNCTN(ISTOL,JSTOM2,KC)  - FUNCTN(ISTOM1,JSTOM2,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOM1,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM2,JSTOM2,KC) - FUNCTN(ISTOM1,JSTOM2,KC)
     +            - FUNCTN(ISTOM2,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM3,JSTOM2,KC) - FUNCTN(ISTOM1,JSTOM2,KC)
     +            - FUNCTN(ISTOM3,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM4,JSTOM2,KC) - FUNCTN(ISTOM1,JSTOM2,KC)
     +            - FUNCTN(ISTOM4,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
            FDIFFC =
     +      ACF2XY*(FUNCTN(ISTOL,JSTOM3,KC)  - FUNCTN(ISTOM1,JSTOM3,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOM1,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM2,JSTOM3,KC) - FUNCTN(ISTOM1,JSTOM3,KC)
     +            - FUNCTN(ISTOM2,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM3,JSTOM3,KC) - FUNCTN(ISTOM1,JSTOM3,KC)
     +            - FUNCTN(ISTOM3,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM4,JSTOM3,KC) - FUNCTN(ISTOM1,JSTOM3,KC)
     +            - FUNCTN(ISTOM4,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
            FDIFFD =
     +      ACF2XY*(FUNCTN(ISTOL,JSTOM4,KC)  - FUNCTN(ISTOM1,JSTOM4,KC)
     +            - FUNCTN(ISTOL,JSTOL,KC)   + FUNCTN(ISTOM1,JSTOL,KC))
     +    + BCF2XY*(FUNCTN(ISTOM2,JSTOM4,KC) - FUNCTN(ISTOM1,JSTOM4,KC)
     +            - FUNCTN(ISTOM2,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + CCF2XY*(FUNCTN(ISTOM3,JSTOM4,KC) - FUNCTN(ISTOM1,JSTOM4,KC)
     +            - FUNCTN(ISTOM3,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
     +    + DCF2XY*(FUNCTN(ISTOM4,JSTOM4,KC) - FUNCTN(ISTOM1,JSTOM4,KC)
     +            - FUNCTN(ISTOM4,JSTOL,KC)  + FUNCTN(ISTOM1,JSTOL,KC))
            FDERIV(ISTOM1,JSTOL,KC) = ACF1XY*FDIFFA
     +                              + BCF1XY*FDIFFB
     +                              + CCF1XY*FDIFFC
     +                              + DCF1XY*FDIFFD

C           RH EDGE IN Y
            DO IC = ISTOM4,ISTOM2

              ICM2 = IC-2
              ICM1 = IC-1
              ICP1 = IC+1
              ICP2 = IC+2

C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOL,KC)  + FUNCTN(ICM1,JSTOL,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOL,KC)  + FUNCTN(ICM2,JSTOL,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM2,KC) + FUNCTN(ICM1,JSTOM2,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM2,KC) + FUNCTN(ICM2,JSTOM2,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM3,KC) + FUNCTN(ICM1,JSTOM3,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM3,KC) + FUNCTN(ICM2,JSTOM3,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTOM1,KC) - FUNCTN(ICM1,JSTOM1,KC)
     +                - FUNCTN(ICP1,JSTOM4,KC) + FUNCTN(ICM1,JSTOM4,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOM1,KC) - FUNCTN(ICM2,JSTOM1,KC)
     +                - FUNCTN(ICP2,JSTOM4,KC) + FUNCTN(ICM2,JSTOM4,KC))
              FDERIV(IC,JSTOM1,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM1,KC) + FUNCTN(ICM1,JSTOM1,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM1,KC) + FUNCTN(ICM2,JSTOM1,KC))
              FDIFFB =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM2,KC) + FUNCTN(ICM1,JSTOM2,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM2,KC) + FUNCTN(ICM2,JSTOM2,KC))
              FDIFFC =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM3,KC) + FUNCTN(ICM1,JSTOM3,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM3,KC) + FUNCTN(ICM2,JSTOM3,KC))
              FDIFFD =
     +          ACOFX1*(FUNCTN(ICP1,JSTOL,KC)  - FUNCTN(ICM1,JSTOL,KC)
     +                - FUNCTN(ICP1,JSTOM4,KC) + FUNCTN(ICM1,JSTOM4,KC))
     +        + BCOFX1*(FUNCTN(ICP2,JSTOL,KC)  - FUNCTN(ICM2,JSTOL,KC)
     +                - FUNCTN(ICP2,JSTOM4,KC) + FUNCTN(ICM2,JSTOM4,KC))
              FDERIV(IC,JSTOL,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD

            ENDDO

C           RH EDGE IN X
            DO JC = JSTOM4,JSTOM2

              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

C             RH POINT: 4TH ORDER ONE-SIDED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM1,JCP1,KC) + FUNCTN(ISTOM1,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM1,JCP2,KC) + FUNCTN(ISTOM1,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM2,JCP1,KC) + FUNCTN(ISTOM2,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM2,JCP2,KC) + FUNCTN(ISTOM2,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM3,JCP1,KC) + FUNCTN(ISTOM3,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM3,JCP2,KC) + FUNCTN(ISTOM3,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTOL,JCP1,KC)  - FUNCTN(ISTOL,JCM1,KC)
     +                - FUNCTN(ISTOM4,JCP1,KC) + FUNCTN(ISTOM4,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOL,JCP2,KC)  - FUNCTN(ISTOL,JCM2,KC)
     +                - FUNCTN(ISTOM4,JCP2,KC) + FUNCTN(ISTOM4,JCM2,KC))
              FDERIV(ISTOL,JC,KC) = ACF1XY*FDIFFA
     +                            + BCF1XY*FDIFFB
     +                            + CCF1XY*FDIFFC
     +                            + DCF1XY*FDIFFD
     
C             RH POINT MINUS 1: 4TH ORDER MIXED/CENTRED
              FDIFFA =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOL,JCP1,KC)  + FUNCTN(ISTOL,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOL,JCP2,KC)  + FUNCTN(ISTOL,JCM2,KC))
              FDIFFB =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM2,JCP1,KC) + FUNCTN(ISTOM2,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM2,JCP2,KC) + FUNCTN(ISTOM2,JCM2,KC))
              FDIFFC =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM3,JCP1,KC) + FUNCTN(ISTOM3,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM3,JCP2,KC) + FUNCTN(ISTOM3,JCM2,KC))
              FDIFFD =
     +          ACOFY1*(FUNCTN(ISTOM1,JCP1,KC) - FUNCTN(ISTOM1,JCM1,KC)
     +                - FUNCTN(ISTOM4,JCP1,KC) + FUNCTN(ISTOM4,JCM1,KC))
     +        + BCOFY1*(FUNCTN(ISTOM1,JCP2,KC) - FUNCTN(ISTOM1,JCM2,KC)
     +                - FUNCTN(ISTOM4,JCP2,KC) + FUNCTN(ISTOM4,JCM2,KC))
              FDERIV(ISTOM1,JC,KC) = ACF2XY*FDIFFA
     +                             + BCF2XY*FDIFFB
     +                             + CCF2XY*FDIFFC
     +                             + DCF2XY*FDIFFD

            ENDDO

C           INTERIOR POINTS 4TH ORDER
            JS = 0
            DO JC = JSTOM4,JSTOM2

              JS = JS+1
              JCM2 = JC-2
              JCM1 = JC-1
              JCP1 = JC+1
              JCP2 = JC+2

              IS = 0
              DO IC = ISTOM4,ISTOM2

                IS = IS+1
                ICM2 = IC-2
                ICM1 = IC-1
                ICP1 = IC+1
                ICP2 = IC+2

C               4TH ORDER CENTRED
                FDIFFA = FUNCTN(ICP1,JCP1,KC) - FUNCTN(ICP1,JCM1,KC)
     +                 - FUNCTN(ICM1,JCP1,KC) + FUNCTN(ICM1,JCM1,KC)
                FDIFFB = FUNCTN(ICP2,JCP2,KC) - FUNCTN(ICP2,JCM2,KC)
     +                 - FUNCTN(ICM2,JCP2,KC) + FUNCTN(ICM2,JCM2,KC)
                FDERIV(IC,JC,KC) = ACF3XY*FDIFFA
     +                           + BCF3XY*FDIFFB
                FSTORA(IS,JS) = FDIFFA
                FSTORB(IS,JS) = FDIFFB

              ENDDO
            ENDDO

C           INTERIOR POINTS 6TH ORDER
            JS = 0
            DO JC = JSTOM4,JSTOM3

              JS = JS+1
              JCM3 = JC-3
              JCP3 = JC+3

              IS = 0
              DO IC = ISTOM4,ISTOM3

                IS = IS+1
                ICM3 = IC-3
                ICP3 = IC+3

C               6TH ORDER CENTRED
                FDIFFC = FUNCTN(ICP3,JCP3,KC) - FUNCTN(ICP3,JCM3,KC)
     +                 - FUNCTN(ICM3,JCP3,KC) + FUNCTN(ICM3,JCM3,KC)
                FDERIV(IC,JC,KC) = ACF4XY*FSTORA(IS,JS)
     +                           + BCF4XY*FSTORB(IS,JS)
     +                           + CCF4XY*FDIFFC
                FSTORC(IS,JS) = FDIFFC

              ENDDO
            ENDDO

C           INTERIOR POINT 8TH ORDER
            JS = 1
            IS = 1
            JC = JSTOM4
            IC = ISTOM4
            JCM4 = JC-4
            JCP4 = JC+4
            ICM4 = IC-4
            ICP4 = IC+4

C           8TH ORDER CENTRED
            FDIFFD = FUNCTN(ICP4,JCP4,KC) - FUNCTN(ICP4,JCM4,KC)
     +             - FUNCTN(ICM4,JCP4,KC) + FUNCTN(ICM4,JCM4,KC)
            FDERIV(IC,JC,KC) = ACF5XY*FSTORA(IS,JS)
     +                       + BCF5XY*FSTORB(IS,JS)
     +                       + CCF5XY*FSTORC(IS,JS)
     +                       + DCF5XY*FDIFFD

          ENDDO

        ENDIF

      ENDIF

C     =========================================================================

C     SCALING
C     =======
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDELX*OVDELY

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
