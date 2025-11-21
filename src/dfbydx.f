      SUBROUTINE DFBYDX(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     DFBYDX
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     28-MAR-2003:  RSC MODIFIED FOR SENGA2
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES FIRST X-DERIVATIVE OF SPECIFIED FUNCTION
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
      DOUBLE PRECISION FDIFFA,FDIFFB,FDIFFC,FDIFFD,FDIFFE
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

            FDIFFA = FUNCTN(ICP1,JC,KC) - FUNCTN(ICM1,JC,KC) 
            FDIFFB = FUNCTN(ICP2,JC,KC) - FUNCTN(ICM2,JC,KC) 
            FDIFFC = FUNCTN(ICP3,JC,KC) - FUNCTN(ICM3,JC,KC) 
            FDIFFD = FUNCTN(ICP4,JC,KC) - FUNCTN(ICM4,JC,KC) 
            FDIFFE = FUNCTN(ICP5,JC,KC) - FUNCTN(ICM5,JC,KC) 

            FDERIV(IC,JC,KC) = ACOFFX*FDIFFA
     +                       + BCOFFX*FDIFFB
     +                       + CCOFFX*FDIFFC
     +                       + DCOFFX*FDIFFD
     +                       + ECOFFX*FDIFFE

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
            FDIFFA = FUNCTN(ISTAP1,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFFB = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFFC = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDIFFD = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAL,JC,KC) = ACOF1X*FDIFFA
     +                          + BCOF1X*FDIFFB
     +                          + CCOF1X*FDIFFC
     +                          + DCOF1X*FDIFFD

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTAL,JC,KC)  - FUNCTN(ISTAP1,JC,KC) 
            FDIFFB = FUNCTN(ISTAP2,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFFC = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFFD = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDERIV(ISTAP1,JC,KC) = ACOF2X*FDIFFA
     +                           + BCOF2X*FDIFFB
     +                           + CCOF2X*FDIFFC
     +                           + DCOF2X*FDIFFD

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP3,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFFB = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP2,JC,KC) = ACOF3X*FDIFFA
     +                           + BCOF3X*FDIFFB

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP4,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFFB = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFFC = FUNCTN(ISTAP6,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP3,JC,KC) = ACOF4X*FDIFFA
     +                           + BCOF4X*FDIFFB
     +                           + CCOF4X*FDIFFC

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTAP5,JC,KC) - FUNCTN(ISTAP3,JC,KC) 
            FDIFFB = FUNCTN(ISTAP6,JC,KC) - FUNCTN(ISTAP2,JC,KC) 
            FDIFFC = FUNCTN(ISTAP7,JC,KC) - FUNCTN(ISTAP1,JC,KC) 
            FDIFFD = FUNCTN(ISTAP8,JC,KC) - FUNCTN(ISTAL,JC,KC) 
            FDERIV(ISTAP4,JC,KC) = ACOF5X*FDIFFA
     +                           + BCOF5X*FDIFFB
     +                           + CCOF5X*FDIFFC
     +                           + DCOF5X*FDIFFD
      
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
            FDIFFA = FUNCTN(ISTOM3,JC,KC) - FUNCTN(ISTOM5,JC,KC) 
            FDIFFB = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM6,JC,KC) 
            FDIFFC = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM7,JC,KC) 
            FDIFFD = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM8,JC,KC) 
            FDERIV(ISTOM4,JC,KC) = ACOF5X*FDIFFA
     +                           + BCOF5X*FDIFFB
     +                           + CCOF5X*FDIFFC
     +                           + DCOF5X*FDIFFD
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM2,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDIFFB = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM5,JC,KC) 
            FDIFFC = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM6,JC,KC) 
            FDERIV(ISTOM3,JC,KC) = ACOF4X*FDIFFA
     +                           + BCOF4X*FDIFFB
     +                           + CCOF4X*FDIFFC
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFFB = FUNCTN(ISTOL,JC,KC)  - FUNCTN(ISTOM4,JC,KC) 
            FDERIV(ISTOM2,JC,KC) = ACOF3X*FDIFFA
     +                           + BCOF3X*FDIFFB
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOL,JC,KC) 
            FDIFFB = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM2,JC,KC) 
            FDIFFC = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFFD = FUNCTN(ISTOM1,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDERIV(ISTOM1,JC,KC) = ACOF2X*FDIFFA
     +                           + BCOF2X*FDIFFB
     +                           + CCOF2X*FDIFFC
     +                           + DCOF2X*FDIFFD
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFFA = FUNCTN(ISTOL,JC,KC) - FUNCTN(ISTOM1,JC,KC) 
            FDIFFB = FUNCTN(ISTOL,JC,KC) - FUNCTN(ISTOM2,JC,KC) 
            FDIFFC = FUNCTN(ISTOL,JC,KC) - FUNCTN(ISTOM3,JC,KC) 
            FDIFFD = FUNCTN(ISTOL,JC,KC) - FUNCTN(ISTOM4,JC,KC) 
            FDERIV(ISTOL,JC,KC) = ACOF1X*FDIFFA
     +                          + BCOF1X*FDIFFB
     +                          + CCOF1X*FDIFFC
     +                          + DCOF1X*FDIFFD

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     SCALING
C     =======
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDELX

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
