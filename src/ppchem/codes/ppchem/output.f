      SUBROUTINE OUTPUT
 
C     *************************************************************************
C
C     OUTPUT
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     22-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     WRITES OUT THE CHEMICAL DATA FOR PPCHEM
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
      INTEGER ISPEC
      INTEGER ISTEP
      INTEGER IBODY
      INTEGER ILIND,ITROE,ISRIF
      INTEGER IC,JC
 
      CHARACTER*50 ERRSTR


C     BEGIN
C     =====

C     =========================================================================

      WRITE(NCREPT,*)'Please enter the name of the output file:'
      READ(NCINPT,'(A)')FNOUTP
      WRITE(NCREPT,'(A)')FNOUTP
      WRITE(NCREPT,*)
      OPEN(UNIT=NCOUTP,FILE=FNOUTP,STATUS='UNKNOWN',FORM='FORMATTED',
     +     ERR=1010)

C     =========================================================================

C     HEADER
      WRITE(NCOUTP,*)'*****************************'
      WRITE(NCOUTP,*)'*                           *'
      WRITE(NCOUTP,*)'*  Output file from ppchem  *'
      WRITE(NCOUTP,*)'*                           *'
      WRITE(NCOUTP,*)'*****************************'

C     =========================================================================

C     SPECIES LIST
      WRITE(NCOUTP,*)'Species list:'
      WRITE(NCOUTP,'(I5)')NSPEC
      DO ISPEC = 1, NSPEC
        WRITE(NCOUTP,'(I5,3X,A)')ISPEC,SPCSYM(ISPEC)
      ENDDO

C     SPECIES DATA
      WRITE(NCOUTP,*)'Species data:'
      WRITE(NCOUTP,'(1PE12.4)')PREFGB
      DO ISPEC = 1,NSPEC
        WRITE(NCOUTP,'(I5,1PE12.4)')ISPEC,WMOLAR(ISPEC)
        WRITE(NCOUTP,'(1PE12.4)')CLEWIS(ISPEC)
        WRITE(NCOUTP,'(I5)')NTINT(ISPEC)
        DO IC = 1,NTINT(ISPEC)
          WRITE(NCOUTP,'(2(1PE12.4),I5)')TINTLO(IC,ISPEC),
     +                  TINTHI(IC,ISPEC),NCOFCP(IC,ISPEC)
          DO JC = 1,NCOFCP(IC,ISPEC)
            WRITE(NCOUTP,'(1PE15.7)')ACOFCP(JC,IC,ISPEC)
          ENDDO
        ENDDO
      ENDDO

C     =========================================================================

C     STEP RATE DATA
      WRITE(NCOUTP,*)'Step rate data:'
      WRITE(NCOUTP,'(I5)')NSTEP
      DO ISTEP = 1, NSTEP
        WRITE(NCOUTP,'(I5,3(1PE12.4))')ISTEP,(RPARAM(IC,ISTEP),IC=1,3)
      ENDDO

C     =========================================================================

C     STEP SPECIES-LIST
      WRITE(NCOUTP,*)'Step species-list:'
      DO ISTEP = 1,NSTEP
        WRITE(NCOUTP,'(2I5)')ISTEP,NSSLEN(ISTEP)
        DO ISPEC = 1, NSSLEN(ISTEP)
          WRITE(NCOUTP,'(2I5,I8)')ISTEP,ISPEC,NSSPEC(ISPEC,ISTEP)
        ENDDO
      ENDDO

C     STEP REACTANT-LIST
      WRITE(NCOUTP,*)'Step reactant-list:'
      DO ISTEP = 1,NSTEP
        WRITE(NCOUTP,'(2I5)')ISTEP,NRSLEN(ISTEP)
        DO ISPEC = 1, NRSLEN(ISTEP)
          WRITE(NCOUTP,'(2I5,I8)')ISTEP,ISPEC,NRSPEC(ISPEC,ISTEP)
        ENDDO
      ENDDO

C     STEP PRODUCT-LIST
      WRITE(NCOUTP,*)'Step product-list:'
      DO ISTEP = 1,NSTEP
        WRITE(NCOUTP,'(2I5)')ISTEP,NPSLEN(ISTEP)
        DO ISPEC = 1, NPSLEN(ISTEP)
          WRITE(NCOUTP,'(2I5,I8)')ISTEP,ISPEC,NPSPEC(ISPEC,ISTEP)
        ENDDO
      ENDDO

C     STEP REACTANT COEFFICIENT-LIST
      WRITE(NCOUTP,*)'Step reactant coefficient-list:'
      DO ISTEP = 1,NSTEP
        WRITE(NCOUTP,'(2I5)')ISTEP,NRCLEN(ISTEP)
        DO ISPEC = 1, NRCLEN(ISTEP)
          WRITE(NCOUTP,'(2I5,I8,1PE12.4)')ISTEP,ISPEC,
     +                     NRCPEC(ISPEC,ISTEP),CRSPEC(ISPEC,ISTEP)
        ENDDO
      ENDDO

C     STEP PRODUCT COEFFICIENT-LIST
      WRITE(NCOUTP,*)'Step product coefficient-list:'
      DO ISTEP = 1,NSTEP
        WRITE(NCOUTP,'(2I5)')ISTEP,NPCLEN(ISTEP)
        DO ISPEC = 1, NPCLEN(ISTEP)
          WRITE(NCOUTP,'(2I5,I8,1PE12.4)')ISTEP,ISPEC,
     +                     NPCPEC(ISPEC,ISTEP),CPSPEC(ISPEC,ISTEP)
        ENDDO
      ENDDO

C     SPECIES DELTA-LIST
      WRITE(NCOUTP,*)'Species delta-list:'
      DO ISTEP = 1,NSTEP
        DO ISPEC = 1,NSSLEN(ISTEP)
          WRITE(NCOUTP,'(2I5,F5.1)')ISTEP,ISPEC,DIFFMU(ISPEC,ISTEP)
        ENDDO
      ENDDO

C     =========================================================================

C     THIRD-BODY LIST
      WRITE(NCOUTP,*)'Third-body list:'
      WRITE(NCOUTP,'(I5)')NBODY
      DO IBODY = 1, NBODY
        WRITE(NCOUTP,'(I5,3X,A)')IBODY,BDYSYM(IBODY)
      ENDDO

C     THIRD-BODY STEP-LIST
      WRITE(NCOUTP,*)'Third-body step-list:'
      IF(NBODY.GT.0)THEN
        DO ISTEP = 1, NSTEP
          WRITE(NCOUTP,'(I5,I8)')ISTEP,MBLIST(ISTEP)
        ENDDO
      ENDIF

C     THIRD-BODY EFFICIENCIES
      WRITE(NCOUTP,*)'Third-body efficiencies:'
      DO IBODY = 1, NBODY
        DO ISPEC = 1, NSPEC
          WRITE(NCOUTP,'(2I5,1PE12.4)')IBODY,ISPEC,EFFY3B(ISPEC,IBODY)
        ENDDO
      ENDDO

C     =========================================================================

C     GIBBS STEP-LIST
      WRITE(NCOUTP,*)'Gibbs step-list:'
      WRITE(NCOUTP,'(I5)')NGIBB
      IF(NGIBB.GT.0)THEN
        DO ISTEP = 1, NSTEP
          WRITE(NCOUTP,'(I5,I8)')ISTEP,MGLIST(ISTEP)
        ENDDO
      ENDIF

C     =========================================================================

C     LINDEMANN STEP-LIST
      WRITE(NCOUTP,*)'Lindemann step-list:'
      WRITE(NCOUTP,'(I5)')NLIND
      IF(NLIND.GT.0)THEN
        DO ISTEP = 1, NSTEP
          WRITE(NCOUTP,'(I5,I8)')ISTEP,MLLIST(ISTEP)
        ENDDO
      ENDIF

C     LINDEMANN STEP RATE DATA
      WRITE(NCOUTP,*)'Lindemann step rate data:'
      DO ILIND = 1, NLIND
        WRITE(NCOUTP,'(I5,4(1PE12.4))')ILIND,(RCLIND(IC,ILIND),IC=1,4)
      ENDDO

C     =========================================================================

C     TROE STEP-LIST
      WRITE(NCOUTP,*)'Troe step-list:'
      WRITE(NCOUTP,'(I5)')NTROE
      IF(NTROE.GT.0)THEN
        DO ISTEP = 1, NSTEP
          WRITE(NCOUTP,'(I5,I8)')ISTEP,MTLIST(ISTEP)
        ENDDO
      ENDIF

C     TROE STEP RATE DATA
      WRITE(NCOUTP,*)'Troe step rate data:'
      DO ITROE = 1, NTROE
        WRITE(NCOUTP,'(I5,6(1PE12.4))')ITROE,(RCTROE(IC,ITROE),IC=1,6)
        WRITE(NCOUTP,'(I5,6(1PE12.4))')ITROE,(RCTROE(IC,ITROE),IC=7,12)
      ENDDO

C     =========================================================================

C     SRI STEP-LIST
      WRITE(NCOUTP,*)'SRI step-list:'
      WRITE(NCOUTP,'(I5)')NSRIF
      IF(NSRIF.GT.0)THEN
        DO ISTEP = 1, NSTEP
          WRITE(NCOUTP,'(I5,I8)')ISTEP,MSLIST(ISTEP)
        ENDDO
      ENDIF

C     SRI STEP RATE DATA
      WRITE(NCOUTP,*)'SRI step rate data:'
      DO ISRIF = 1, NSRIF
        WRITE(NCOUTP,'(I5,3(1PE12.4))')ISRIF,(RCSRIF(IC,ISRIF),IC=1,3)
        WRITE(NCOUTP,'(I5,5(1PE12.4))')ISRIF,(RCSRIF(IC,ISRIF),IC=4,8)
      ENDDO

C     =========================================================================

      WRITE(NCOUTP,*)'End of file'
      CLOSE(NCOUTP)

C     =========================================================================

C     FILENAME ERROR HANDLER
      GOTO 1020
1010  ERRSTR = 'opening output file:'
      CALL ERHAND(ERRSTR,IEFATL)
1020  CONTINUE
C     END OF FILENAME ERROR HANDLER

C     =========================================================================


      RETURN
      END
