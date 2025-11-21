      SUBROUTINE ARDUMS(ARRAY,NXPHYS,NYPHYS,NZPHYS,NSPEC,ISPEC)
 
C     *************************************************************************
C
C     ARDUMS
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     11-SEP-2005:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     DIAGNOSTIC ROUTINE
C     PRINTS THE CONTENTS OF THE SPECIFIED SPECIES ARRAY
C
C     *************************************************************************
 

C     ARGUMENTS
C     =========
      INTEGER NXPHYS,NYPHYS,NZPHYS,NSPEC,ISPEC
      DOUBLE PRECISION ARRAY(NXPHYS,NYPHYS,NZPHYS,NSPEC)


C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC


C     BEGIN
C     =====

      WRITE(6,*)
      DO KC = NZPHYS,1,-1
        WRITE(6,'(" Z-PLANE ",I10)')KC
        DO JC = NYPHYS,1,-1
          WRITE(6,'(I5,21(1PE10.3))')
     +      JC,(ARRAY(IC,JC,KC,ISPEC),IC = 1,NXPHYS)
        ENDDO
      ENDDO


      RETURN
      END
