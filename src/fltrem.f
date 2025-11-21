      SUBROUTINE FLTREM
 
C     *************************************************************************
C
C     FLTREM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-AUG-2009:  CREATED
C     08-AUG-2012:  RSC EVALUATE ALL SPECIES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     CARRIES OUT SPATIAL FILTERING
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC,ISPEC


C     BEGIN
C     =====
      
C     =========================================================================

C     DENSITY
      CALL FILTRX(DRHS)
      CALL FILTRY(DRHS)
      CALL FILTRZ(DRHS)

C     U-VELOCITY
      CALL FILTRX(URHS)
      CALL FILTRY(URHS)
      CALL FILTRZ(URHS)

C     V-VELOCITY
      CALL FILTRX(VRHS)
      CALL FILTRY(VRHS)
      CALL FILTRZ(VRHS)

C     W-VELOCITY
      CALL FILTRX(WRHS)
      CALL FILTRY(WRHS)
      CALL FILTRZ(WRHS)

C     INTERNAL ENERGY
      CALL FILTRX(ERHS)
      CALL FILTRY(ERHS)
      CALL FILTRZ(ERHS)

C     SPECIES MASS FRACTION
C     RSC 08-AUG-2012 EVALUATE ALL SPECIES
C      DO ISPEC = 1,NSPM1
      DO ISPEC = 1,NSPEC

        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              STORE7(IC,JC,KC) = YRHS(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO

        CALL FILTRX(STORE7)
        CALL FILTRY(STORE7)
        CALL FILTRZ(STORE7)

        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB

              YRHS(IC,JC,KC,ISPEC) = STORE7(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     =========================================================================


      RETURN
      END
