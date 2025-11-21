      SUBROUTINE FILTRX(FUNCTN)
 
C     *************************************************************************
C
C     FILTRX
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
C     X DIRECTION
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
      INTEGER ISTART,IFINIS
      INTEGER ICM6,ICM5,ICM4,ICM3,ICM2,ICM1,ICCC
      INTEGER ICP1,ICP2,ICP3,ICP4,ICP5,ICP6


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      ISTART = ISTAL
      IFINIS = ISTOL
      IF(NENDXL.EQ.NBOUND)ISTART = ISTAP6
      IF(NENDXR.EQ.NBOUND)IFINIS = ISTOM6

C     =========================================================================

C     INTERIOR SCHEME
C     ===============

C     TWELFTH ORDER EXPLICIT FILTER
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL

          ICM5 = ISTART-6
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
          ICP6 = ISTART+5

          DO IC = ISTART,IFINIS

            ICM6 = ICM5
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
            ICP5 = ICP6
            ICP6 = IC+6

            FDIFFA = FUNCTN(ICP1,JC,KC) + FUNCTN(ICM1,JC,KC) 
            FDIFFB = FUNCTN(ICP2,JC,KC) + FUNCTN(ICM2,JC,KC) 
            FDIFFC = FUNCTN(ICP3,JC,KC) + FUNCTN(ICM3,JC,KC) 
            FDIFFD = FUNCTN(ICP4,JC,KC) + FUNCTN(ICM4,JC,KC) 
            FDIFFE = FUNCTN(ICP5,JC,KC) + FUNCTN(ICM5,JC,KC) 
            FDIFFF = FUNCTN(ICP6,JC,KC) + FUNCTN(ICM6,JC,KC) 

            FILTER(IC,JC,KC) = FACOFX*FDIFFA
     +                       + FBCOFX*FDIFFB
     +                       + FCCOFX*FDIFFC
     +                       + FDCOFX*FDIFFD
     +                       + FECOFX*FDIFFE
     +                       + FFCOFX*FDIFFF
     +                       + FGCOFX*FUNCTN(IC,JC,KC)

          ENDDO

        ENDDO
      ENDDO

C     =========================================================================

C     LH END
C     ======
      IF(NENDXL.EQ.NBOUND)THEN

C       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C           LH POINT: 6TH ORDER ONE-SIDED
            FILTER(ISTAL,JC,KC) = FACF1X*FUNCTN(ISTAL,JC,KC)
     +                          + FBCF1X*FUNCTN(ISTAP1,JC,KC)
     +                          + FCCF1X*FUNCTN(ISTAP2,JC,KC)
     +                          + FDCF1X*FUNCTN(ISTAP3,JC,KC)
     +                          + FECF1X*FUNCTN(ISTAP4,JC,KC)
     +                          + FFCF1X*FUNCTN(ISTAP5,JC,KC)
     +                          + FGCF1X*FUNCTN(ISTAP6,JC,KC)

C           LH POINT PLUS 1: 7TH ORDER MIXED
            FILTER(ISTAP1,JC,KC) = FACF2X*FUNCTN(ISTAL,JC,KC)
     +                           + FBCF2X*FUNCTN(ISTAP1,JC,KC)
     +                           + FCCF2X*FUNCTN(ISTAP2,JC,KC)
     +                           + FDCF2X*FUNCTN(ISTAP3,JC,KC)
     +                           + FECF2X*FUNCTN(ISTAP4,JC,KC)
     +                           + FFCF2X*FUNCTN(ISTAP5,JC,KC)
     +                           + FGCF2X*FUNCTN(ISTAP6,JC,KC)
     +                           + FHCF2X*FUNCTN(ISTAP7,JC,KC)

C           LH POINT PLUS 2: 8TH ORDER MIXED
            FILTER(ISTAP2,JC,KC) = FACF3X*FUNCTN(ISTAL,JC,KC)
     +                           + FBCF3X*FUNCTN(ISTAP1,JC,KC)
     +                           + FCCF3X*FUNCTN(ISTAP2,JC,KC)
     +                           + FDCF3X*FUNCTN(ISTAP3,JC,KC)
     +                           + FECF3X*FUNCTN(ISTAP4,JC,KC)
     +                           + FFCF3X*FUNCTN(ISTAP5,JC,KC)
     +                           + FGCF3X*FUNCTN(ISTAP6,JC,KC)
     +                           + FHCF3X*FUNCTN(ISTAP7,JC,KC)
     +                           + FICF3X*FUNCTN(ISTAP8,JC,KC)

C           LH POINT PLUS 3: 9TH ORDER MIXED
            FILTER(ISTAP3,JC,KC) = FACF4X*FUNCTN(ISTAL,JC,KC)
     +                           + FBCF4X*FUNCTN(ISTAP1,JC,KC)
     +                           + FCCF4X*FUNCTN(ISTAP2,JC,KC)
     +                           + FDCF4X*FUNCTN(ISTAP3,JC,KC)
     +                           + FECF4X*FUNCTN(ISTAP4,JC,KC)
     +                           + FFCF4X*FUNCTN(ISTAP5,JC,KC)
     +                           + FGCF4X*FUNCTN(ISTAP6,JC,KC)
     +                           + FHCF4X*FUNCTN(ISTAP7,JC,KC)
     +                           + FICF4X*FUNCTN(ISTAP8,JC,KC)
     +                           + FJCF4X*FUNCTN(ISTAP9,JC,KC)

C           LH POINT PLUS 4: 10TH ORDER MIXED
            FILTER(ISTAP4,JC,KC) = FACF5X*FUNCTN(ISTAL,JC,KC)
     +                           + FBCF5X*FUNCTN(ISTAP1,JC,KC)
     +                           + FCCF5X*FUNCTN(ISTAP2,JC,KC)
     +                           + FDCF5X*FUNCTN(ISTAP3,JC,KC)
     +                           + FECF5X*FUNCTN(ISTAP4,JC,KC)
     +                           + FFCF5X*FUNCTN(ISTAP5,JC,KC)
     +                           + FGCF5X*FUNCTN(ISTAP6,JC,KC)
     +                           + FHCF5X*FUNCTN(ISTAP7,JC,KC)
     +                           + FICF5X*FUNCTN(ISTAP8,JC,KC)
     +                           + FJCF5X*FUNCTN(ISTAP9,JC,KC)
     +                           + FKCF5X*FUNCTN(ISTAPA,JC,KC)
      
C           LH POINT PLUS 5: 11TH ORDER MIXED
            FILTER(ISTAP5,JC,KC) = FACF6X*FUNCTN(ISTAL,JC,KC)
     +                           + FBCF6X*FUNCTN(ISTAP1,JC,KC)
     +                           + FCCF6X*FUNCTN(ISTAP2,JC,KC)
     +                           + FDCF6X*FUNCTN(ISTAP3,JC,KC)
     +                           + FECF6X*FUNCTN(ISTAP4,JC,KC)
     +                           + FFCF6X*FUNCTN(ISTAP5,JC,KC)
     +                           + FGCF6X*FUNCTN(ISTAP6,JC,KC)
     +                           + FHCF6X*FUNCTN(ISTAP7,JC,KC)
     +                           + FICF6X*FUNCTN(ISTAP8,JC,KC)
     +                           + FJCF6X*FUNCTN(ISTAP9,JC,KC)
     +                           + FKCF6X*FUNCTN(ISTAPA,JC,KC)
     +                           + FLCF6X*FUNCTN(ISTAPB,JC,KC)
      
          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END
C     ======
      IF(NENDXR.EQ.NBOUND)THEN

C       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

C           RH POINT MINUS 5: 11TH ORDER MIXED
            FILTER(ISTOM5,JC,KC) = FACF6X*FUNCTN(ISTOL,JC,KC)
     +                           + FBCF6X*FUNCTN(ISTOM1,JC,KC)
     +                           + FCCF6X*FUNCTN(ISTOM2,JC,KC)
     +                           + FDCF6X*FUNCTN(ISTOM3,JC,KC)
     +                           + FECF6X*FUNCTN(ISTOM4,JC,KC)
     +                           + FFCF6X*FUNCTN(ISTOM5,JC,KC)
     +                           + FGCF6X*FUNCTN(ISTOM6,JC,KC)
     +                           + FHCF6X*FUNCTN(ISTOM7,JC,KC)
     +                           + FICF6X*FUNCTN(ISTOM8,JC,KC)
     +                           + FJCF6X*FUNCTN(ISTOM9,JC,KC)
     +                           + FKCF6X*FUNCTN(ISTOMA,JC,KC)
     +                           + FLCF6X*FUNCTN(ISTOMB,JC,KC)
      
C           RH POINT MINUS 4: 10TH ORDER MIXED
            FILTER(ISTOM4,JC,KC) = FACF5X*FUNCTN(ISTOL,JC,KC)
     +                           + FBCF5X*FUNCTN(ISTOM1,JC,KC)
     +                           + FCCF5X*FUNCTN(ISTOM2,JC,KC)
     +                           + FDCF5X*FUNCTN(ISTOM3,JC,KC)
     +                           + FECF5X*FUNCTN(ISTOM4,JC,KC)
     +                           + FFCF5X*FUNCTN(ISTOM5,JC,KC)
     +                           + FGCF5X*FUNCTN(ISTOM6,JC,KC)
     +                           + FHCF5X*FUNCTN(ISTOM7,JC,KC)
     +                           + FICF5X*FUNCTN(ISTOM8,JC,KC)
     +                           + FJCF5X*FUNCTN(ISTOM9,JC,KC)
     +                           + FKCF5X*FUNCTN(ISTOMA,JC,KC)
      
C           RH POINT MINUS 3: 9TH ORDER MIXED
            FILTER(ISTOM3,JC,KC) = FACF4X*FUNCTN(ISTOL,JC,KC)
     +                           + FBCF4X*FUNCTN(ISTOM1,JC,KC)
     +                           + FCCF4X*FUNCTN(ISTOM2,JC,KC)
     +                           + FDCF4X*FUNCTN(ISTOM3,JC,KC)
     +                           + FECF4X*FUNCTN(ISTOM4,JC,KC)
     +                           + FFCF4X*FUNCTN(ISTOM5,JC,KC)
     +                           + FGCF4X*FUNCTN(ISTOM6,JC,KC)
     +                           + FHCF4X*FUNCTN(ISTOM7,JC,KC)
     +                           + FICF4X*FUNCTN(ISTOM8,JC,KC)
     +                           + FJCF4X*FUNCTN(ISTOM9,JC,KC)

C           RH POINT MINUS 2: 8TH ORDER MIXED
            FILTER(ISTOM2,JC,KC) = FACF3X*FUNCTN(ISTOL,JC,KC)
     +                           + FBCF3X*FUNCTN(ISTOM1,JC,KC)
     +                           + FCCF3X*FUNCTN(ISTOM2,JC,KC)
     +                           + FDCF3X*FUNCTN(ISTOM3,JC,KC)
     +                           + FECF3X*FUNCTN(ISTOM4,JC,KC)
     +                           + FFCF3X*FUNCTN(ISTOM5,JC,KC)
     +                           + FGCF3X*FUNCTN(ISTOM6,JC,KC)
     +                           + FHCF3X*FUNCTN(ISTOM7,JC,KC)
     +                           + FICF3X*FUNCTN(ISTOM8,JC,KC)

C           RH POINT PLUS 1: 7TH ORDER MIXED
            FILTER(ISTOM1,JC,KC) = FACF2X*FUNCTN(ISTOL,JC,KC)
     +                           + FBCF2X*FUNCTN(ISTOM1,JC,KC)
     +                           + FCCF2X*FUNCTN(ISTOM2,JC,KC)
     +                           + FDCF2X*FUNCTN(ISTOM3,JC,KC)
     +                           + FECF2X*FUNCTN(ISTOM4,JC,KC)
     +                           + FFCF2X*FUNCTN(ISTOM5,JC,KC)
     +                           + FGCF2X*FUNCTN(ISTOM6,JC,KC)
     +                           + FHCF2X*FUNCTN(ISTOM7,JC,KC)

C           RH POINT: 6TH ORDER ONE-SIDED
            FILTER(ISTOL,JC,KC) = FACF1X*FUNCTN(ISTOL,JC,KC)
     +                          + FBCF1X*FUNCTN(ISTOM1,JC,KC)
     +                          + FCCF1X*FUNCTN(ISTOM2,JC,KC)
     +                          + FDCF1X*FUNCTN(ISTOM3,JC,KC)
     +                          + FECF1X*FUNCTN(ISTOM4,JC,KC)
     +                          + FFCF1X*FUNCTN(ISTOM5,JC,KC)
     +                          + FGCF1X*FUNCTN(ISTOM6,JC,KC)

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
