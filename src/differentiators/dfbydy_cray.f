      SUBROUTINE DFBYDY(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     DFBYDY
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
C     EVALUATES FIRST Y-DERIVATIVE OF SPECIFIED FUNCTION
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
      DOUBLE PRECISION FDIFFA,FDIFFB,FDIFFC,FDIFFD,FDIFFE
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

            FDIFFA = FUNCTN(IC,JC+1,KC) - FUNCTN(IC,JC-1,KC) 
            FDIFFB = FUNCTN(IC,JC+2,KC) - FUNCTN(IC,JC-2,KC) 
            FDIFFC = FUNCTN(IC,JC+3,KC) - FUNCTN(IC,JC-3,KC) 
            FDIFFD = FUNCTN(IC,JC+4,KC) - FUNCTN(IC,JC-4,KC) 
            FDIFFE = FUNCTN(IC,JC+5,KC) - FUNCTN(IC,JC-5,KC) 

            FDERIV(IC,JC,KC) = ACOFFY*FDIFFA
     +                       + BCOFFY*FDIFFB
     +                       + CCOFFY*FDIFFC
     +                       + DCOFFY*FDIFFD
     +                       + ECOFFY*FDIFFE

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDELY

          ENDDO
        ENDDO

C  `    =======================================================================

C       LH END
C       ======
        IF(NENDYL.EQ.NBOUND)THEN

C         EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFFA = FUNCTN(IC,JSTAP1,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFFB = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFFC = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAL,KC) 
            FDIFFD = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC) 
            FDERIV(IC,JSTAL,KC) = (ACOF1Y*FDIFFA
     +                          +  BCOF1Y*FDIFFB
     +                          +  CCOF1Y*FDIFFC
     +                          +  DCOF1Y*FDIFFD)*OVDELY

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTAL,KC)  - FUNCTN(IC,JSTAP1,KC) 
            FDIFFB = FUNCTN(IC,JSTAP2,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFC = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFD = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDERIV(IC,JSTAP1,KC) = (ACOF2Y*FDIFFA
     +                           +  BCOF2Y*FDIFFB
     +                           +  CCOF2Y*FDIFFC
     +                           +  DCOF2Y*FDIFFD)*OVDELY

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP3,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFB = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAL,KC) 
            FDERIV(IC,JSTAP2,KC) = (ACOF3Y*FDIFFA
     +                           +  BCOF3Y*FDIFFB)*OVDELY

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP4,KC) - FUNCTN(IC,JSTAP2,KC)
            FDIFFB = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFC = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAL,KC) 
            FDERIV(IC,JSTAP3,KC) = (ACOF4Y*FDIFFA
     +                           +  BCOF4Y*FDIFFB
     +                           +  CCOF4Y*FDIFFC)*OVDELY

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTAP5,KC) - FUNCTN(IC,JSTAP3,KC) 
            FDIFFB = FUNCTN(IC,JSTAP6,KC) - FUNCTN(IC,JSTAP2,KC) 
            FDIFFC = FUNCTN(IC,JSTAP7,KC) - FUNCTN(IC,JSTAP1,KC) 
            FDIFFD = FUNCTN(IC,JSTAP8,KC) - FUNCTN(IC,JSTAL,KC) 
            FDERIV(IC,JSTAP4,KC) = (ACOF5Y*FDIFFA
     +                           +  BCOF5Y*FDIFFB
     +                           +  CCOF5Y*FDIFFC
     +                           +  DCOF5Y*FDIFFD)*OVDELY
      
          ENDDO

        ENDIF 

C       =======================================================================

C       RH END
C       ======
        IF(NENDYR.EQ.NBOUND)THEN

C         EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM3,KC) - FUNCTN(IC,JSTOM5,KC) 
            FDIFFB = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM6,KC) 
            FDIFFC = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM7,KC) 
            FDIFFD = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM8,KC) 
            FDERIV(IC,JSTOM4,KC) = (ACOF5Y*FDIFFA
     +                           +  BCOF5Y*FDIFFB
     +                           +  CCOF5Y*FDIFFC
     +                           +  DCOF5Y*FDIFFD)*OVDELY
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM2,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDIFFB = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM5,KC) 
            FDIFFC = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM6,KC) 
            FDERIV(IC,JSTOM3,KC) = (ACOF4Y*FDIFFA
     +                           +  BCOF4Y*FDIFFB
     +                           +  CCOF4Y*FDIFFC)*OVDELY
      
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFFB = FUNCTN(IC,JSTOL,KC)  - FUNCTN(IC,JSTOM4,KC) 
            FDERIV(IC,JSTOM2,KC) = (ACOF3Y*FDIFFA
     +                           +  BCOF3Y*FDIFFB)*OVDELY
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOL,KC) 
            FDIFFB = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM2,KC) 
            FDIFFC = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFFD = FUNCTN(IC,JSTOM1,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDERIV(IC,JSTOM1,KC) = (ACOF2Y*FDIFFA
     +                           +  BCOF2Y*FDIFFB
     +                           +  CCOF2Y*FDIFFC
     +                           +  DCOF2Y*FDIFFD)*OVDELY
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFFA = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM1,KC) 
            FDIFFB = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM2,KC) 
            FDIFFC = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM3,KC) 
            FDIFFD = FUNCTN(IC,JSTOL,KC) - FUNCTN(IC,JSTOM4,KC) 
            FDERIV(IC,JSTOL,KC) = (ACOF1Y*FDIFFA
     +                          +  BCOF1Y*FDIFFB
     +                          +  CCOF1Y*FDIFFC
     +                          +  DCOF1Y*FDIFFD)*OVDELY
      
          ENDDO

        ENDIF

C       =======================================================================

      ENDDO

C     =========================================================================


      RETURN
      END
