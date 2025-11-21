      SUBROUTINE LISTEM
 
C     *************************************************************************
C
C     LISTEM
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     23-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     CREATES THE LISTS FOR PPCHEM
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_chemin.h'
      INCLUDE 'com_ppccom.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      INTEGER ISPEC,JSPEC,ISTEP
      INTEGER IC
      INTEGER MUDIFF
 
      CHARACTER*50 ERRSTR


C     BEGIN
C     =====

C     =========================================================================

C     CREATE STEP SPECIES-LIST
C     ------------------------
      DO ISTEP = 1,NSTEP
        NSSLEN(ISTEP) = 0
        DO ISPEC = 1, NSPEC
        
          IF((NRTABL(ISPEC,ISTEP).NE.0)
     +      .OR.(NPTABL(ISPEC,ISTEP).NE.0))THEN
            NSSLEN(ISTEP) = NSSLEN(ISTEP) + 1
            IF(NSSLEN(ISTEP).LE.NSSMAX)THEN
              NSSPEC(NSSLEN(ISTEP),ISTEP) = NUMSPC(ISPEC)  
            ELSE
              ERRSTR = 'step species-list length exceeded'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
          ENDIF

        ENDDO
      ENDDO

C     =========================================================================

C     CREATE STEP REACTANT-LIST
C     -------------------------
C     WITH MULTIPLE ENTRIES FOR ANY SPECIES WITH MU > 1
      DO ISTEP = 1,NSTEP
        NRSLEN(ISTEP) = 0
        DO ISPEC = 1, NSPEC
        
          IF(NRTABL(ISPEC,ISTEP).NE.0)THEN
            DO IC = 1, NRTABL(ISPEC,ISTEP)
              NRSLEN(ISTEP) = NRSLEN(ISTEP) + 1
              IF(NRSLEN(ISTEP).LE.NRSMAX)THEN
                NRSPEC(NRSLEN(ISTEP),ISTEP) = NUMSPC(ISPEC)  
              ELSE
                ERRSTR = 'step reactant-list length exceeded'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
            ENDDO
          ENDIF

        ENDDO
      ENDDO

C     =========================================================================

C     CREATE STEP PRODUCT-LIST
C     ------------------------
C     WITH MULTIPLE ENTRIES FOR ANY SPECIES WITH MU > 1
      DO ISTEP = 1,NSTEP
        NPSLEN(ISTEP) = 0
        DO ISPEC = 1, NSPEC
        
          IF(NPTABL(ISPEC,ISTEP).NE.0)THEN
            DO IC = 1, NPTABL(ISPEC,ISTEP)
              NPSLEN(ISTEP) = NPSLEN(ISTEP) + 1
              IF(NPSLEN(ISTEP).LE.NRSMAX)THEN
                NPSPEC(NPSLEN(ISTEP),ISTEP) = NUMSPC(ISPEC)  
              ELSE
                ERRSTR = 'step product-list length exceeded'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
            ENDDO
          ENDIF

        ENDDO
      ENDDO

C     =========================================================================

C     CREATE STEP REACTANT COEFFICIENT-LIST
C     -------------------------------------
      DO ISTEP = 1,NSTEP
        NRCLEN(ISTEP) = 0
        DO ISPEC = 1, NSPEC
        
          IF(ABS(CRTABL(ISPEC,ISTEP)).GT.CSCTOL)THEN
            NRCLEN(ISTEP) = NRCLEN(ISTEP) + 1
            IF(NRCLEN(ISTEP).LE.NRSMAX)THEN
              NRCPEC(NRCLEN(ISTEP),ISTEP) = NUMSPC(ISPEC)
              CRSPEC(NRCLEN(ISTEP),ISTEP) = CRTABL(ISPEC,ISTEP)
            ELSE
              ERRSTR = 'step reactant coefficient-list length exceeded'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
          ENDIF

        ENDDO
      ENDDO

C     =========================================================================

C     CREATE STEP PRODUCT COEFFICIENT-LIST
C     ------------------------------------
      DO ISTEP = 1,NSTEP
        NPCLEN(ISTEP) = 0
        DO ISPEC = 1, NSPEC
        
          IF(ABS(CPTABL(ISPEC,ISTEP)).GT.CSCTOL)THEN
            NPCLEN(ISTEP) = NPCLEN(ISTEP) + 1
            IF(NPCLEN(ISTEP).LE.NRSMAX)THEN
              NPCPEC(NPCLEN(ISTEP),ISTEP) = NUMSPC(ISPEC)
              CPSPEC(NPCLEN(ISTEP),ISTEP) = CPTABL(ISPEC,ISTEP)
            ELSE
              ERRSTR = 'step product coefficient-list length exceeded'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
          ENDIF

        ENDDO
      ENDDO

C     =========================================================================

C     CREATE SPECIES DELTA-LIST
C     -------------------------

      DO ISTEP = 1, NSTEP
        DO ISPEC = 1, NSSLEN(ISTEP)

          JSPEC = NSSPEC(ISPEC,ISTEP) 
          MUDIFF = NPTABL(JSPEC,ISTEP)-NRTABL(JSPEC,ISTEP)

          IF(MUDIFF.EQ.0)THEN
            WRITE(6,'(2I5)')ISTEP,JSPEC
            ERRSTR = 'species delta-list entry is zero'
            CALL ERHAND(ERRSTR,IEWARN)
          ENDIF

          DIFFMU(ISPEC,ISTEP) = REAL(MUDIFF)

        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
