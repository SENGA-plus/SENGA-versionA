      SUBROUTINE RADCAL
 
C     *************************************************************************
C
C     RADCAL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     14-JUL-2013:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     RADIATION TREATMENT
C     USING OPTICALLY THIN ASSUMPTION: Ju et al: JFM 342, 315-334, 1997.
C     AFTER TOM DUNSTAN 2012
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION PLSPEC,FORNOW
      INTEGER IC,JC,KC,ISPEC,JSPEC,ICP


C     BEGIN
C     =====

C     =========================================================================

C     BUILD THE PLANCK MEAN ABSORPTION COEFFICIENT OF THE MIXTURE
C     -----------------------------------------------------------

C     INITIALISE THE ACCUMULATOR
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            STORE1(IC,JC,KC) = ZERO

          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------

C     RUN THROUGH ALL RADIATING SPECIES
      DO JSPEC = 1, NSPRAD

C       PLANCK MEAN ABSORPTION COEFFICIENT OF EACH SPECIES
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = TRUN(IC,JC,KC)
              PLSPEC = AKPRAD(NKPRAD(JSPEC),JSPEC)
              DO ICP = NKPRM1(JSPEC),1,-1
                PLSPEC = PLSPEC*FORNOW + AKPRAD(ICP,JSPEC)
              ENDDO
              STORE2(IC,JC,KC) = PLSPEC

            ENDDO
          ENDDO
        ENDDO

C       SPECIES ID
        ISPEC = NSPRID(JSPEC)

C       ADD THE SPECIES CONTRIBUTION
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              FORNOW = YRHS(IC,JC,KC,ISPEC)*RGSPEC(ISPEC)*TRUN(IC,JC,KC)
              STORE1(IC,JC,KC) = STORE1(IC,JC,KC)
     +                         + STORE2(IC,JC,KC)*FORNOW

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     =========================================================================

C     INCLUDE THE RADIATION TERM IN THE ENERGY EQUATION

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FORNOW = TRUN(IC,JC,KC)
            FORNOW = FORNOW*FORNOW*FORNOW*FORNOW

            ERHS(IC,JC,KC) = ERHS(IC,JC,KC)
     +                     - FOURSB*STORE1(IC,JC,KC)*(FORNOW - TRFRTH)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
