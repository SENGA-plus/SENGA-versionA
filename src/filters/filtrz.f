      SUBROUTINE FILTRZ(FUNCTN)
 
C     *************************************************************************
C
C     FILTRZ
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
C     Z DIRECTION
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
      INTEGER KSTART,KFINIS
      INTEGER KCM6,KCM5,KCM4,KCM3,KCM2,KCM1,KCCC
      INTEGER KCP1,KCP2,KCP3,KCP4,KCP5,KCP6


C     BEGIN
C     =====

C     =========================================================================

C     END CONDITIONS
C     ==============

      KSTART = KSTAL
      KFINIS = KSTOL
      IF(NENDZL.EQ.NBOUND)KSTART = KSTAP6
      IF(NENDZR.EQ.NBOUND)KFINIS = KSTOM6

C     =========================================================================

C     INTERIOR SCHEME
C     ===============

C     TWELFTH ORDER EXPLICIT FILTER
      KCM5 = KSTART-6
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
      KCP6 = KSTART+5

      DO KC = KSTART,KFINIS

        KCM6 = KCM5
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
        KCP5 = KCP6
        KCP6 = KC+6

        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FDIFFA = FUNCTN(IC,JC,KCP1) + FUNCTN(IC,JC,KCM1) 
            FDIFFB = FUNCTN(IC,JC,KCP2) + FUNCTN(IC,JC,KCM2) 
            FDIFFC = FUNCTN(IC,JC,KCP3) + FUNCTN(IC,JC,KCM3) 
            FDIFFD = FUNCTN(IC,JC,KCP4) + FUNCTN(IC,JC,KCM4) 
            FDIFFE = FUNCTN(IC,JC,KCP5) + FUNCTN(IC,JC,KCM5) 
            FDIFFF = FUNCTN(IC,JC,KCP6) + FUNCTN(IC,JC,KCM6) 

            FILTER(IC,JC,KC) = FACOFZ*FDIFFA
     +                       + FBCOFZ*FDIFFB
     +                       + FCCOFZ*FDIFFC
     +                       + FDCOFZ*FDIFFD
     +                       + FECOFZ*FDIFFE
     +                       + FFCOFZ*FDIFFF
     +                       + FGCOFZ*FUNCTN(IC,JC,KC)

          ENDDO

        ENDDO
      ENDDO

C     =========================================================================

C     LH END
C     ======
      IF(NENDZL.EQ.NBOUND)THEN

C       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

C           LH POINT: 6TH ORDER ONE-SIDED
            FILTER(IC,JC,KSTAL) = FACF1Z*FUNCTN(IC,JC,KSTAL)
     +                          + FBCF1Z*FUNCTN(IC,JC,KSTAP1)
     +                          + FCCF1Z*FUNCTN(IC,JC,KSTAP2)
     +                          + FDCF1Z*FUNCTN(IC,JC,KSTAP3)
     +                          + FECF1Z*FUNCTN(IC,JC,KSTAP4)
     +                          + FFCF1Z*FUNCTN(IC,JC,KSTAP5)
     +                          + FGCF1Z*FUNCTN(IC,JC,KSTAP6)

C           LH POINT PLUS 1: 7TH ORDER MIXED
            FILTER(IC,JC,KSTAP1) = FACF2Z*FUNCTN(IC,JC,KSTAL)
     +                           + FBCF2Z*FUNCTN(IC,JC,KSTAP1)
     +                           + FCCF2Z*FUNCTN(IC,JC,KSTAP2)
     +                           + FDCF2Z*FUNCTN(IC,JC,KSTAP3)
     +                           + FECF2Z*FUNCTN(IC,JC,KSTAP4)
     +                           + FFCF2Z*FUNCTN(IC,JC,KSTAP5)
     +                           + FGCF2Z*FUNCTN(IC,JC,KSTAP6)
     +                           + FHCF2Z*FUNCTN(IC,JC,KSTAP7)

C           LH POINT PLUS 2: 8TH ORDER MIXED
            FILTER(IC,JC,KSTAP2) = FACF3Z*FUNCTN(IC,JC,KSTAL)
     +                           + FBCF3Z*FUNCTN(IC,JC,KSTAP1)
     +                           + FCCF3Z*FUNCTN(IC,JC,KSTAP2)
     +                           + FDCF3Z*FUNCTN(IC,JC,KSTAP3)
     +                           + FECF3Z*FUNCTN(IC,JC,KSTAP4)
     +                           + FFCF3Z*FUNCTN(IC,JC,KSTAP5)
     +                           + FGCF3Z*FUNCTN(IC,JC,KSTAP6)
     +                           + FHCF3Z*FUNCTN(IC,JC,KSTAP7)
     +                           + FICF3Z*FUNCTN(IC,JC,KSTAP8)

C           LH POINT PLUS 3: 9TH ORDER MIXED
            FILTER(IC,JC,KSTAP3) = FACF4Z*FUNCTN(IC,JC,KSTAL)
     +                           + FBCF4Z*FUNCTN(IC,JC,KSTAP1)
     +                           + FCCF4Z*FUNCTN(IC,JC,KSTAP2)
     +                           + FDCF4Z*FUNCTN(IC,JC,KSTAP3)
     +                           + FECF4Z*FUNCTN(IC,JC,KSTAP4)
     +                           + FFCF4Z*FUNCTN(IC,JC,KSTAP5)
     +                           + FGCF4Z*FUNCTN(IC,JC,KSTAP6)
     +                           + FHCF4Z*FUNCTN(IC,JC,KSTAP7)
     +                           + FICF4Z*FUNCTN(IC,JC,KSTAP8)
     +                           + FJCF4Z*FUNCTN(IC,JC,KSTAP9)

C           LH POINT PLUS 4: 10TH ORDER MIXED
            FILTER(IC,JC,KSTAP4) = FACF5Z*FUNCTN(IC,JC,KSTAL)
     +                           + FBCF5Z*FUNCTN(IC,JC,KSTAP1)
     +                           + FCCF5Z*FUNCTN(IC,JC,KSTAP2)
     +                           + FDCF5Z*FUNCTN(IC,JC,KSTAP3)
     +                           + FECF5Z*FUNCTN(IC,JC,KSTAP4)
     +                           + FFCF5Z*FUNCTN(IC,JC,KSTAP5)
     +                           + FGCF5Z*FUNCTN(IC,JC,KSTAP6)
     +                           + FHCF5Z*FUNCTN(IC,JC,KSTAP7)
     +                           + FICF5Z*FUNCTN(IC,JC,KSTAP8)
     +                           + FJCF5Z*FUNCTN(IC,JC,KSTAP9)
     +                           + FKCF5Z*FUNCTN(IC,JC,KSTAPA)
      
C           LH POINT PLUS 5: 11TH ORDER MIXED
            FILTER(IC,JC,KSTAP5) = FACF6Z*FUNCTN(IC,JC,KSTAL)
     +                           + FBCF6Z*FUNCTN(IC,JC,KSTAP1)
     +                           + FCCF6Z*FUNCTN(IC,JC,KSTAP2)
     +                           + FDCF6Z*FUNCTN(IC,JC,KSTAP3)
     +                           + FECF6Z*FUNCTN(IC,JC,KSTAP4)
     +                           + FFCF6Z*FUNCTN(IC,JC,KSTAP5)
     +                           + FGCF6Z*FUNCTN(IC,JC,KSTAP6)
     +                           + FHCF6Z*FUNCTN(IC,JC,KSTAP7)
     +                           + FICF6Z*FUNCTN(IC,JC,KSTAP8)
     +                           + FJCF6Z*FUNCTN(IC,JC,KSTAP9)
     +                           + FKCF6Z*FUNCTN(IC,JC,KSTAPA)
     +                           + FLCF6Z*FUNCTN(IC,JC,KSTAPB)
      
          ENDDO
        ENDDO

      ENDIF 

C     =========================================================================

C     RH END
C     ======
      IF(NENDZR.EQ.NBOUND)THEN

C       EXPLICIT 6TH,7TH,8TH,9TH,10TH,11TH ORDER BOUNDARY TREATMENT
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 5: 11TH ORDER MIXED
            FILTER(IC,JC,KSTOM5) = FACF6Z*FUNCTN(IC,JC,KSTOL)
     +                           + FBCF6Z*FUNCTN(IC,JC,KSTOM1)
     +                           + FCCF6Z*FUNCTN(IC,JC,KSTOM2)
     +                           + FDCF6Z*FUNCTN(IC,JC,KSTOM3)
     +                           + FECF6Z*FUNCTN(IC,JC,KSTOM4)
     +                           + FFCF6Z*FUNCTN(IC,JC,KSTOM5)
     +                           + FGCF6Z*FUNCTN(IC,JC,KSTOM6)
     +                           + FHCF6Z*FUNCTN(IC,JC,KSTOM7)
     +                           + FICF6Z*FUNCTN(IC,JC,KSTOM8)
     +                           + FJCF6Z*FUNCTN(IC,JC,KSTOM9)
     +                           + FKCF6Z*FUNCTN(IC,JC,KSTOMA)
     +                           + FLCF6Z*FUNCTN(IC,JC,KSTOMB)
      
C           RH POINT MINUS 4: 10TH ORDER MIXED
            FILTER(IC,JC,KSTOM4) = FACF5Z*FUNCTN(IC,JC,KSTOL)
     +                           + FBCF5Z*FUNCTN(IC,JC,KSTOM1)
     +                           + FCCF5Z*FUNCTN(IC,JC,KSTOM2)
     +                           + FDCF5Z*FUNCTN(IC,JC,KSTOM3)
     +                           + FECF5Z*FUNCTN(IC,JC,KSTOM4)
     +                           + FFCF5Z*FUNCTN(IC,JC,KSTOM5)
     +                           + FGCF5Z*FUNCTN(IC,JC,KSTOM6)
     +                           + FHCF5Z*FUNCTN(IC,JC,KSTOM7)
     +                           + FICF5Z*FUNCTN(IC,JC,KSTOM8)
     +                           + FJCF5Z*FUNCTN(IC,JC,KSTOM9)
     +                           + FKCF5Z*FUNCTN(IC,JC,KSTOMA)
      
C           RH POINT MINUS 3: 9TH ORDER MIXED
            FILTER(IC,JC,KSTOM3) = FACF4Z*FUNCTN(IC,JC,KSTOL)
     +                           + FBCF4Z*FUNCTN(IC,JC,KSTOM1)
     +                           + FCCF4Z*FUNCTN(IC,JC,KSTOM2)
     +                           + FDCF4Z*FUNCTN(IC,JC,KSTOM3)
     +                           + FECF4Z*FUNCTN(IC,JC,KSTOM4)
     +                           + FFCF4Z*FUNCTN(IC,JC,KSTOM5)
     +                           + FGCF4Z*FUNCTN(IC,JC,KSTOM6)
     +                           + FHCF4Z*FUNCTN(IC,JC,KSTOM7)
     +                           + FICF4Z*FUNCTN(IC,JC,KSTOM8)
     +                           + FJCF4Z*FUNCTN(IC,JC,KSTOM9)

C           RH POINT MINUS 2: 8TH ORDER MIXED
            FILTER(IC,JC,KSTOM2) = FACF3Z*FUNCTN(IC,JC,KSTOL)
     +                           + FBCF3Z*FUNCTN(IC,JC,KSTOM1)
     +                           + FCCF3Z*FUNCTN(IC,JC,KSTOM2)
     +                           + FDCF3Z*FUNCTN(IC,JC,KSTOM3)
     +                           + FECF3Z*FUNCTN(IC,JC,KSTOM4)
     +                           + FFCF3Z*FUNCTN(IC,JC,KSTOM5)
     +                           + FGCF3Z*FUNCTN(IC,JC,KSTOM6)
     +                           + FHCF3Z*FUNCTN(IC,JC,KSTOM7)
     +                           + FICF3Z*FUNCTN(IC,JC,KSTOM8)

C           RH POINT PLUS 1: 7TH ORDER MIXED
            FILTER(IC,JC,KSTOM1) = FACF2Z*FUNCTN(IC,JC,KSTOL)
     +                           + FBCF2Z*FUNCTN(IC,JC,KSTOM1)
     +                           + FCCF2Z*FUNCTN(IC,JC,KSTOM2)
     +                           + FDCF2Z*FUNCTN(IC,JC,KSTOM3)
     +                           + FECF2Z*FUNCTN(IC,JC,KSTOM4)
     +                           + FFCF2Z*FUNCTN(IC,JC,KSTOM5)
     +                           + FGCF2Z*FUNCTN(IC,JC,KSTOM6)
     +                           + FHCF2Z*FUNCTN(IC,JC,KSTOM7)

C           RH POINT: 6TH ORDER ONE-SIDED
            FILTER(IC,JC,KSTOL) = FACF1Z*FUNCTN(IC,JC,KSTOL)
     +                          + FBCF1Z*FUNCTN(IC,JC,KSTOM1)
     +                          + FCCF1Z*FUNCTN(IC,JC,KSTOM2)
     +                          + FDCF1Z*FUNCTN(IC,JC,KSTOM3)
     +                          + FECF1Z*FUNCTN(IC,JC,KSTOM4)
     +                          + FFCF1Z*FUNCTN(IC,JC,KSTOM5)
     +                          + FGCF1Z*FUNCTN(IC,JC,KSTOM6)

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
