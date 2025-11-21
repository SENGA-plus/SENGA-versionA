      SUBROUTINE DFBYDZ(FUNCTN,FDERIV)
 
C     *************************************************************************
C
C     DFBYDZ
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
C     EVALUATES FIRST Z-DERIVATIVE OF SPECIFIED FUNCTION
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

            FDIFFA = FUNCTN(IC,JC,KC+1) - FUNCTN(IC,JC,KC-1) 
            FDIFFB = FUNCTN(IC,JC,KC+2) - FUNCTN(IC,JC,KC-2) 
            FDIFFC = FUNCTN(IC,JC,KC+3) - FUNCTN(IC,JC,KC-3) 
            FDIFFD = FUNCTN(IC,JC,KC+4) - FUNCTN(IC,JC,KC-4) 
            FDIFFE = FUNCTN(IC,JC,KC+5) - FUNCTN(IC,JC,KC-5) 

            FDERIV(IC,JC,KC) = ACOFFZ*FDIFFA
     +                       + BCOFFZ*FDIFFB
     +                       + CCOFFZ*FDIFFC
     +                       + DCOFFZ*FDIFFD
     +                       + ECOFFZ*FDIFFE

            FDERIV(IC,JC,KC) = FDERIV(IC,JC,KC)*OVDELZ

          ENDDO
        ENDDO

C       =======================================================================

C       LH END
C       ======
        IF(NENDZL.EQ.NBOUND)THEN

C         EXPLICIT 4TH,4TH,4TH,6TH,8TH ORDER BOUNDARY TREATMENT
          DO IC = ISTAL,ISTOL

C           LH POINT: 4TH ORDER ONE-SIDED
            FDIFFA = FUNCTN(IC,JC,KSTAP1) - FUNCTN(IC,JC,KSTAL) 
            FDIFFB = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAL) 
            FDIFFC = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAL) 
            FDIFFD = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAL) = (ACOF1Z*FDIFFA
     +                          +  BCOF1Z*FDIFFB
     +                          +  CCOF1Z*FDIFFC
     +                          +  DCOF1Z*FDIFFD)*OVDELZ

C           LH POINT PLUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JC,KSTAL)  - FUNCTN(IC,JC,KSTAP1) 
            FDIFFB = FUNCTN(IC,JC,KSTAP2) - FUNCTN(IC,JC,KSTAP1) 
            FDIFFC = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP1) 
            FDIFFD = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP1) 
            FDERIV(IC,JC,KSTAP1) = (ACOF2Z*FDIFFA
     +                           +  BCOF2Z*FDIFFB
     +                           +  CCOF2Z*FDIFFC
     +                           +  DCOF2Z*FDIFFD)*OVDELZ

C           LH POINT PLUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JC,KSTAP3) - FUNCTN(IC,JC,KSTAP1) 
            FDIFFB = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP2) = (ACOF3Z*FDIFFA
     +                           +  BCOF3Z*FDIFFB)*OVDELZ

C           LH POINT PLUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JC,KSTAP4) - FUNCTN(IC,JC,KSTAP2) 
            FDIFFB = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP1) 
            FDIFFC = FUNCTN(IC,JC,KSTAP6) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP3) = (ACOF4Z*FDIFFA
     +                           +  BCOF4Z*FDIFFB
     +                           +  CCOF4Z*FDIFFC)*OVDELZ

C           LH POINT PLUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JC,KSTAP5) - FUNCTN(IC,JC,KSTAP3) 
            FDIFFB = FUNCTN(IC,JC,KSTAP6) - FUNCTN(IC,JC,KSTAP2) 
            FDIFFC = FUNCTN(IC,JC,KSTAP7) - FUNCTN(IC,JC,KSTAP1) 
            FDIFFD = FUNCTN(IC,JC,KSTAP8) - FUNCTN(IC,JC,KSTAL) 
            FDERIV(IC,JC,KSTAP4) = (ACOF5Z*FDIFFA
     +                           +  BCOF5Z*FDIFFB
     +                           +  CCOF5Z*FDIFFC
     +                           +  DCOF5Z*FDIFFD)*OVDELZ
      
          ENDDO

        ENDIF 

C       =======================================================================

C       RH END
C       ======
        IF(NENDZR.EQ.NBOUND)THEN

          DO IC = ISTAL,ISTOL

C           RH POINT MINUS 4: 8TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JC,KSTOM3) - FUNCTN(IC,JC,KSTOM5) 
            FDIFFB = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM6) 
            FDIFFC = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM7) 
            FDIFFD = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM8) 
            FDERIV(IC,JC,KSTOM4) = (ACOF5Z*FDIFFA
     +                           +  BCOF5Z*FDIFFB
     +                           +  CCOF5Z*FDIFFC
     +                           +  DCOF5Z*FDIFFD)*OVDELZ
      
C           RH POINT MINUS 3: 6TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JC,KSTOM2) - FUNCTN(IC,JC,KSTOM4) 
            FDIFFB = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM5) 
            FDIFFC = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM6) 
            FDERIV(IC,JC,KSTOM3) = (ACOF4Z*FDIFFA
     +                           +  BCOF4Z*FDIFFB
     +                           +  CCOF4Z*FDIFFC)*OVDELZ
 
C           RH POINT MINUS 2: 4TH ORDER CENTRED
            FDIFFA = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM3) 
            FDIFFB = FUNCTN(IC,JC,KSTOL)  - FUNCTN(IC,JC,KSTOM4) 
            FDERIV(IC,JC,KSTOM2) = (ACOF3Z*FDIFFA
     +                           +  BCOF3Z*FDIFFB)*OVDELZ
      
C           RH POINT MINUS 1: 4TH ORDER MIXED
            FDIFFA = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOL) 
            FDIFFB = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM2) 
            FDIFFC = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM3) 
            FDIFFD = FUNCTN(IC,JC,KSTOM1) - FUNCTN(IC,JC,KSTOM4) 
            FDERIV(IC,JC,KSTOM1) = (ACOF2Z*FDIFFA
     +                           +  BCOF2Z*FDIFFB
     +                           +  CCOF2Z*FDIFFC
     +                           +  DCOF2Z*FDIFFD)*OVDELZ
      
C           RH POINT: 4TH ORDER ONE-SIDED
            FDIFFA = FUNCTN(IC,JC,KSTOL) - FUNCTN(IC,JC,KSTOM1) 
            FDIFFB = FUNCTN(IC,JC,KSTOL) - FUNCTN(IC,JC,KSTOM2) 
            FDIFFC = FUNCTN(IC,JC,KSTOL) - FUNCTN(IC,JC,KSTOM3) 
            FDIFFD = FUNCTN(IC,JC,KSTOL) - FUNCTN(IC,JC,KSTOM4) 
            FDERIV(IC,JC,KSTOL) = (ACOF1Z*FDIFFA
     +                          +  BCOF1Z*FDIFFB
     +                          +  CCOF1Z*FDIFFC
     +                          +  DCOF1Z*FDIFFD)*OVDELZ

          ENDDO

        ENDIF

C       =======================================================================

      ENDDO

C     =========================================================================


      RETURN
      END
