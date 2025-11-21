      SUBROUTINE FILTRY(FUNCTN)
 
C     *************************************************************************
C
C     FILTRY
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     30-AUG-2009:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EXPLICIT 12TH ORDER FINITE DIFFERENCE FILTER
C     WITH EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER END CONDITIONS
C     Y DIRECTION
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


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION FILTER(NXSIZE,NYSIZE,NZSIZE)
      DOUBLE PRECISION FDIFFA,FDIFFB,FDIFFC,FDIFFD,FDIFFE,FDIFFF
      INTEGER IC,JC,KC
      INTEGER JSTART,JFINIS
      INTEGER JCM6,JCM5,JCM4,JCM3,JCM2,JCM1,JCCC
      INTEGER JCP1,JCP2,JCP3,JCP4,JCP5,JCP6


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      JSTART = JSTAL
      JFINIS = JSTOL
      IF(NENDYL.EQ.NBOUND)JSTART = JSTAP6
      IF(NENDYR.EQ.NBOUND)JFINIS = JSTOM6

C     =========================================================================

C     INTERIOR SCHEME
C     ===============

C     TWELFTH ORDER EXPLICIT FILTER
      DO KC = KSTAL,KSTOL

        JCM5 = JSTART-6
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
        JCP6 = JSTART+5

        DO JC = JSTART,JFINIS

          JCM6 = JCM5
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
          JCP5 = JCP6
          JCP6 = JC+6

          DO IC = ISTAL,ISTOL

            FDIFFA = FUNCTN(IC,JCP1,KC) + FUNCTN(IC,JCM1,KC) 
            FDIFFB = FUNCTN(IC,JCP2,KC) + FUNCTN(IC,JCM2,KC) 
            FDIFFC = FUNCTN(IC,JCP3,KC) + FUNCTN(IC,JCM3,KC) 
            FDIFFD = FUNCTN(IC,JCP4,KC) + FUNCTN(IC,JCM4,KC) 
            FDIFFE = FUNCTN(IC,JCP5,KC) + FUNCTN(IC,JCM5,KC) 
            FDIFFF = FUNCTN(IC,JCP6,KC) + FUNCTN(IC,JCM6,KC) 

            FILTER(IC,JC,KC) = FACOFY*FDIFFA
     +                       + FBCOFY*FDIFFB
     +                       + FCCOFY*FDIFFC
     +                       + FDCOFY*FDIFFD
     +                       + FECOFY*FDIFFE
     +                       + FFCOFY*FDIFFF
     +                       + FGCOFY*FUNCTN(IC,JC,KC)

          ENDDO

        ENDDO
      ENDDO

C     =========================================================================

C     LH END
C     ======
      IF(NENDYL.EQ.NBOUND)THEN

C       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C           LH POINT: 6TH ORDER ONE-SIDED
            FILTER(IC,JSTAL,KC) = FACF1Y*FUNCTN(IC,JSTAL,KC)
     +                          + FBCF1Y*FUNCTN(IC,JSTAP1,KC)
     +                          + FCCF1Y*FUNCTN(IC,JSTAP2,KC)
     +                          + FDCF1Y*FUNCTN(IC,JSTAP3,KC)
     +                          + FECF1Y*FUNCTN(IC,JSTAP4,KC)
     +                          + FFCF1Y*FUNCTN(IC,JSTAP5,KC)
     +                          + FGCF1Y*FUNCTN(IC,JSTAP6,KC)

C           LH POINT PLUS 1: 7TH ORDER MIXED
            FILTER(IC,JSTAP1,KC) = FACF2Y*FUNCTN(IC,JSTAL,KC)
     +                           + FBCF2Y*FUNCTN(IC,JSTAP1,KC)
     +                           + FCCF2Y*FUNCTN(IC,JSTAP2,KC)
     +                           + FDCF2Y*FUNCTN(IC,JSTAP3,KC)
     +                           + FECF2Y*FUNCTN(IC,JSTAP4,KC)
     +                           + FFCF2Y*FUNCTN(IC,JSTAP5,KC)
     +                           + FGCF2Y*FUNCTN(IC,JSTAP6,KC)
     +                           + FHCF2Y*FUNCTN(IC,JSTAP7,KC)

C           LH POINT PLUS 2: 8TH ORDER MIXED
            FILTER(IC,JSTAP2,KC) = FACF3Y*FUNCTN(IC,JSTAL,KC)
     +                           + FBCF3Y*FUNCTN(IC,JSTAP1,KC)
     +                           + FCCF3Y*FUNCTN(IC,JSTAP2,KC)
     +                           + FDCF3Y*FUNCTN(IC,JSTAP3,KC)
     +                           + FECF3Y*FUNCTN(IC,JSTAP4,KC)
     +                           + FFCF3Y*FUNCTN(IC,JSTAP5,KC)
     +                           + FGCF3Y*FUNCTN(IC,JSTAP6,KC)
     +                           + FHCF3Y*FUNCTN(IC,JSTAP7,KC)
     +                           + FICF3Y*FUNCTN(IC,JSTAP8,KC)

C           LH POINT PLUS 3: 9TH ORDER MIXED
            FILTER(IC,JSTAP3,KC) = FACF4Y*FUNCTN(IC,JSTAL,KC)
     +                           + FBCF4Y*FUNCTN(IC,JSTAP1,KC)
     +                           + FCCF4Y*FUNCTN(IC,JSTAP2,KC)
     +                           + FDCF4Y*FUNCTN(IC,JSTAP3,KC)
     +                           + FECF4Y*FUNCTN(IC,JSTAP4,KC)
     +                           + FFCF4Y*FUNCTN(IC,JSTAP5,KC)
     +                           + FGCF4Y*FUNCTN(IC,JSTAP6,KC)
     +                           + FHCF4Y*FUNCTN(IC,JSTAP7,KC)
     +                           + FICF4Y*FUNCTN(IC,JSTAP8,KC)
     +                           + FJCF4Y*FUNCTN(IC,JSTAP9,KC)

C           LH POINT PLUS 4: 10TH ORDER MIXED
            FILTER(IC,JSTAP4,KC) = FACF5Y*FUNCTN(IC,JSTAL,KC)
     +                           + FBCF5Y*FUNCTN(IC,JSTAP1,KC)
     +                           + FCCF5Y*FUNCTN(IC,JSTAP2,KC)
     +                           + FDCF5Y*FUNCTN(IC,JSTAP3,KC)
     +                           + FECF5Y*FUNCTN(IC,JSTAP4,KC)
     +                           + FFCF5Y*FUNCTN(IC,JSTAP5,KC)
     +                           + FGCF5Y*FUNCTN(IC,JSTAP6,KC)
     +                           + FHCF5Y*FUNCTN(IC,JSTAP7,KC)
     +                           + FICF5Y*FUNCTN(IC,JSTAP8,KC)
     +                           + FJCF5Y*FUNCTN(IC,JSTAP9,KC)
     +                           + FKCF5Y*FUNCTN(IC,JSTAPA,KC)
      
C           LH POINT PLUS 5: 11TH ORDER MIXED
            FILTER(IC,JSTAP5,KC) = FACF6Y*FUNCTN(IC,JSTAL,KC)
     +                           + FBCF6Y*FUNCTN(IC,JSTAP1,KC)
     +                           + FCCF6Y*FUNCTN(IC,JSTAP2,KC)
     +                           + FDCF6Y*FUNCTN(IC,JSTAP3,KC)
     +                           + FECF6Y*FUNCTN(IC,JSTAP4,KC)
     +                           + FFCF6Y*FUNCTN(IC,JSTAP5,KC)
     +                           + FGCF6Y*FUNCTN(IC,JSTAP6,KC)
     +                           + FHCF6Y*FUNCTN(IC,JSTAP7,KC)
     +                           + FICF6Y*FUNCTN(IC,JSTAP8,KC)
     +                           + FJCF6Y*FUNCTN(IC,JSTAP9,KC)
     +                           + FKCF6Y*FUNCTN(IC,JSTAPA,KC)
     +                           + FLCF6Y*FUNCTN(IC,JSTAPB,KC)
      
          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END
C     ======
      IF(NENDYR.EQ.NBOUND)THEN

C       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 5: 11TH ORDER MIXED
            FILTER(IC,JSTOM5,KC) = FACF6Y*FUNCTN(IC,JSTOL,KC)
     +                           + FBCF6Y*FUNCTN(IC,JSTOM1,KC)
     +                           + FCCF6Y*FUNCTN(IC,JSTOM2,KC)
     +                           + FDCF6Y*FUNCTN(IC,JSTOM3,KC)
     +                           + FECF6Y*FUNCTN(IC,JSTOM4,KC)
     +                           + FFCF6Y*FUNCTN(IC,JSTOM5,KC)
     +                           + FGCF6Y*FUNCTN(IC,JSTOM6,KC)
     +                           + FHCF6Y*FUNCTN(IC,JSTOM7,KC)
     +                           + FICF6Y*FUNCTN(IC,JSTOM8,KC)
     +                           + FJCF6Y*FUNCTN(IC,JSTOM9,KC)
     +                           + FKCF6Y*FUNCTN(IC,JSTOMA,KC)
     +                           + FLCF6Y*FUNCTN(IC,JSTOMB,KC)
      
C           RH POINT MINUS 4: 10TH ORDER MIXED
            FILTER(IC,JSTOM4,KC) = FACF5Y*FUNCTN(IC,JSTOL,KC)
     +                           + FBCF5Y*FUNCTN(IC,JSTOM1,KC)
     +                           + FCCF5Y*FUNCTN(IC,JSTOM2,KC)
     +                           + FDCF5Y*FUNCTN(IC,JSTOM3,KC)
     +                           + FECF5Y*FUNCTN(IC,JSTOM4,KC)
     +                           + FFCF5Y*FUNCTN(IC,JSTOM5,KC)
     +                           + FGCF5Y*FUNCTN(IC,JSTOM6,KC)
     +                           + FHCF5Y*FUNCTN(IC,JSTOM7,KC)
     +                           + FICF5Y*FUNCTN(IC,JSTOM8,KC)
     +                           + FJCF5Y*FUNCTN(IC,JSTOM9,KC)
     +                           + FKCF5Y*FUNCTN(IC,JSTOMA,KC)
      
C           RH POINT MINUS 3: 9TH ORDER MIXED
            FILTER(IC,JSTOM3,KC) = FACF4Y*FUNCTN(IC,JSTOL,KC)
     +                           + FBCF4Y*FUNCTN(IC,JSTOM1,KC)
     +                           + FCCF4Y*FUNCTN(IC,JSTOM2,KC)
     +                           + FDCF4Y*FUNCTN(IC,JSTOM3,KC)
     +                           + FECF4Y*FUNCTN(IC,JSTOM4,KC)
     +                           + FFCF4Y*FUNCTN(IC,JSTOM5,KC)
     +                           + FGCF4Y*FUNCTN(IC,JSTOM6,KC)
     +                           + FHCF4Y*FUNCTN(IC,JSTOM7,KC)
     +                           + FICF4Y*FUNCTN(IC,JSTOM8,KC)
     +                           + FJCF4Y*FUNCTN(IC,JSTOM9,KC)

C           RH POINT MINUS 2: 8TH ORDER MIXED
            FILTER(IC,JSTOM2,KC) = FACF3Y*FUNCTN(IC,JSTOL,KC)
     +                           + FBCF3Y*FUNCTN(IC,JSTOM1,KC)
     +                           + FCCF3Y*FUNCTN(IC,JSTOM2,KC)
     +                           + FDCF3Y*FUNCTN(IC,JSTOM3,KC)
     +                           + FECF3Y*FUNCTN(IC,JSTOM4,KC)
     +                           + FFCF3Y*FUNCTN(IC,JSTOM5,KC)
     +                           + FGCF3Y*FUNCTN(IC,JSTOM6,KC)
     +                           + FHCF3Y*FUNCTN(IC,JSTOM7,KC)
     +                           + FICF3Y*FUNCTN(IC,JSTOM8,KC)

C           RH POINT PLUS 1: 7TH ORDER MIXED
            FILTER(IC,JSTOM1,KC) = FACF2Y*FUNCTN(IC,JSTOL,KC)
     +                           + FBCF2Y*FUNCTN(IC,JSTOM1,KC)
     +                           + FCCF2Y*FUNCTN(IC,JSTOM2,KC)
     +                           + FDCF2Y*FUNCTN(IC,JSTOM3,KC)
     +                           + FECF2Y*FUNCTN(IC,JSTOM4,KC)
     +                           + FFCF2Y*FUNCTN(IC,JSTOM5,KC)
     +                           + FGCF2Y*FUNCTN(IC,JSTOM6,KC)
     +                           + FHCF2Y*FUNCTN(IC,JSTOM7,KC)

C           RH POINT: 6TH ORDER ONE-SIDED
            FILTER(IC,JSTOL,KC) = FACF1Y*FUNCTN(IC,JSTOL,KC)
     +                          + FBCF1Y*FUNCTN(IC,JSTOM1,KC)
     +                          + FCCF1Y*FUNCTN(IC,JSTOM2,KC)
     +                          + FDCF1Y*FUNCTN(IC,JSTOM3,KC)
     +                          + FECF1Y*FUNCTN(IC,JSTOM4,KC)
     +                          + FFCF1Y*FUNCTN(IC,JSTOM5,KC)
     +                          + FGCF1Y*FUNCTN(IC,JSTOM6,KC)

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     COPY BACK
C     =========
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FUNCTN(IC,JC,KC) = FILTER(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
