      SUBROUTINE ARDUMP(ARRAY,NXPHYS,NYPHYS,NZPHYS)
 
C     *************************************************************************
C
C     ARDUMP
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     23-JUL-1992:  CREATED
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     DIAGNOSTIC ROUTINE
C     PRINTS THE CONTENTS OF THE SPECIFIED ARRAY
C
C     *************************************************************************
 

C     ARGUMENTS
C     =========
      INTEGER NXPHYS,NYPHYS,NZPHYS
      DOUBLE PRECISION ARRAY(NXPHYS,NYPHYS,NZPHYS)


C     LOCAL DATA
C     ==========
      INTEGER IC,JC,KC


C     BEGIN
C     =====

      WRITE(6,*)
      DO KC = NZPHYS,1,-1
        WRITE(6,'(" Z-PLANE ",I10)')KC
        DO JC = NYPHYS,1,-1
          WRITE(6,'(I5,21(1PE10.3))')JC,(ARRAY(IC,JC,KC),IC = 1,NXPHYS)
        ENDDO
      ENDDO


      RETURN
      END
