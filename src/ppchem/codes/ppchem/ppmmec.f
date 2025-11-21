      SUBROUTINE PPMMEC
 
C     *************************************************************************
C
C     PPMMEC
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     25-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     PREPROCESSES CHEMICAL MECHANISM MECHANISM DATA FOR PPCHEM
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_chemin.h'
      INCLUDE 'com_ppccom.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
      INTEGER NUMSTR
      PARAMETER(NUMSTR=12)


C     EXTERNAL FUNCTION
C     =================
      LOGICAL CHKSTR
      EXTERNAL CHKSTR


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION RSTOIC

      INTEGER LENSTR(NUMSTR)
      INTEGER ISPEC
      INTEGER ISTEP
      INTEGER NPARSE
      INTEGER ISTR,JSTR,KSTR,ISTOTL
      INTEGER LENSPC,LENSTO
      INTEGER NSTOIC
      INTEGER ITREAD,ISREAD

      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR
      CHARACTER*10 STOSTR
      CHARACTER*10 SPCSTR

      LOGICAL SFOUND
      LOGICAL RSFLAG,RDFLAG
      LOGICAL RFLAG
      LOGICAL READLN,READTR,READSR
      LOGICAL VALIDS


C     BEGIN
C     =====

C     =========================================================================

C     READ THE MECHANISM LIST
C     -----------------------
C     READ AND IGNORE HEADER LINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE

C     READ AND PARSE EACH LINE OF THE MECHANISM LIST
      NSTEP = 0
      READLN = .FALSE.
      READTR = .FALSE.
      ITREAD = 0
      READSR = .FALSE.
      ISREAD = 0
1000  CONTINUE

C       INCREMENT STEP COUNTER AND CHECK
        NSTEP = NSTEP + 1
        IF(READLN)NSTEP = NSTEP - 1
        IF(READTR)NSTEP = NSTEP - 1
        IF(READSR)NSTEP = NSTEP - 1
        IF(NSTEP.GT.NSTPMX)THEN
          ERRSTR = 'max no. of steps exceeded'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF

C       READ LINE AND INITIALISE PARSER
        READ(NCMECH,'(A)',END=1010)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE

        CALL PARMEC(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)

C       DIAGNOSTICS
        WRITE(6,*)'After PARMEC:'
        WRITE(6,'(4I5)')NPARSE,NUMSTR,ISTR,ISTOTL
        DO KSTR = 1, ISTR
          WRITE(6,'(I5,3X,A)')LENSTR(KSTR),PARSTR(KSTR)(1:LENSTR(KSTR))
        ENDDO

C       CONDITION FOR END OF MECHANISM LIST
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
        IF(VALIDS)THEN

C         MECHANISM LIST FINISHED
          NSTEP = NSTEP - 1

        ELSE

          IF(ISTOTL.EQ.0)THEN

            ERRSTR = 'blank line in mechanism list'
            CALL ERHAND(ERRSTR,IEWARN)

          ELSE

C           PROCESS THE STRINGS

CC           ZERO REACTION ORDER COUNTERS
C            NORDER(NSTEP) = 0
C            RORDER(NSTEP) = 0

C           CONVERT FIRST STRING TO STEP NUMBER
C           CHECK FOR VALID NUMBERS
            VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHINT)
            IF(.NOT.VALIDS)THEN
              ERRSTR = 'step number incorrect'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
            READ(PARSTR(1)(1:LENSTR(1)),*)ISTEP
            IF(ISTEP.NE.NSTEP)THEN
              ERRSTR = 'step number out of order'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
            NUMSTP(NSTEP) = ISTEP

C           CHECK LAST STRING FOR INDICATOR FLAG
            IF(READLN)THEN
              IF(ISTR.NE.5)THEN
                ERRSTR = 'Lindemann data format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
            ELSEIF(READTR)THEN
              IF(ITREAD.EQ.1)THEN
                IF(ISTR.NE.8)THEN
                  ERRSTR = 'Troe data format incorrect'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
              ELSEIF(ITREAD.EQ.2)THEN
                IF(ISTR.NE.6)THEN
                  ERRSTR = 'Troe data format incorrect'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
              ENDIF
            ELSEIF(READSR)THEN
              IF(ISREAD.EQ.1)THEN
                IF(ISTR.NE.4)THEN
                  ERRSTR = 'SRI data format incorrect'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
              ELSEIF(ISREAD.EQ.2)THEN
                IF(ISTR.NE.6)THEN
                  ERRSTR = 'SRI data format incorrect'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
              ENDIF
            ELSE
C             INDICATOR = L FOR LINDEMANN FORM
              IF(PARSTR(ISTR)(1:1).EQ.'L')THEN
                NLIND = NLIND + 1
                IF(NLIND.GT.NLNDMX)THEN
                  ERRSTR = 'too many Lindemann steps'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
                MLLIST(NSTEP) = NLIND
                ISTR = ISTR - 1
              ENDIF 
C             INDICATOR = T FOR TROE FORM
              IF(PARSTR(ISTR)(1:1).EQ.'T')THEN
                NTROE = NTROE + 1
                IF(NTROE.GT.NLNDMX)THEN
                  ERRSTR = 'too many Troe steps'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
                MTLIST(NSTEP) = NTROE
                ISTR = ISTR - 1
              ENDIF 
C             INDICATOR = S FOR SRI FORM
              IF(PARSTR(ISTR)(1:1).EQ.'S')THEN
                NSRIF = NSRIF + 1
                IF(NSRIF.GT.NLNDMX)THEN
                  ERRSTR = 'too many SRI steps'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
                MSLIST(NSTEP) = NSRIF
                ISTR = ISTR - 1
              ENDIF 
            ENDIF 

C           CONVERT SPECIES STRINGS
            IF(.NOT.(READLN.OR.READTR.OR.READSR))THEN

C             ZERO REACTION ORDER COUNTERS
              NORDER(NSTEP) = 0
              RORDER(NSTEP) = 0

              RFLAG = .TRUE.
              RSFLAG = .FALSE.
              DO KSTR = 2, ISTR-3

                IF(PARSTR(KSTR)(1:LENSTR(KSTR)).EQ.'=>')THEN

C                 REACTANT LIST COMPLETE: SKIP ARROW
                  STOSTR = '0'
                  LENSTO = 1
                  SPCSTR = '=>'
                  LENSPC = 2
                  RFLAG = .FALSE.

                ELSE IF(PARSTR(KSTR)(1:LENSTR(KSTR)).EQ.'==')THEN

C                 REACTANT LIST COMPLETE: SKIP ARROW
                  STOSTR = '0'
                  LENSTO = 1
                  SPCSTR = '=='
                  LENSPC = 2
                  RFLAG = .FALSE.

C                 SET GIBBS INDICATOR FOR THIS STEP
                  NGIBB = NGIBB + 1
                  IF(NGIBB.GT.NGIBMX)THEN
                    ERRSTR = 'too many Gibbs steps'
                    CALL ERHAND(ERRSTR,IEFATL)
                  ENDIF
                  MGLIST(NSTEP) = NSTEP

                ELSE

C                 PROCESS THE SPECIES STRING
                  JSTR = 1
                  LENSTO = 0
                  RSFLAG = .FALSE.

C                 CHECK FOR STOICHIOMETRIC COEFFICIENT
                  IF((PARSTR(KSTR)(1:1).GE.'0')
     +              .AND.(PARSTR(KSTR)(1:1).LE.'9'))THEN

C                   STOICHIOMETRIC COEFFICIENT IS PRESENT
1100                CONTINUE
                      RDFLAG = (PARSTR(KSTR)(JSTR:JSTR).EQ.'.')
                      IF(((PARSTR(KSTR)(JSTR:JSTR).GE.'0')
     +                  .AND.(PARSTR(KSTR)(JSTR:JSTR).LE.'9'))
     +                     .OR.RDFLAG)THEN
                        IF(RSFLAG.OR.RDFLAG)RSFLAG = .TRUE.
                        STOSTR(JSTR:JSTR) = PARSTR(KSTR)(JSTR:JSTR) 
                        LENSTO = JSTR
                        IF(JSTR.LT.LENSTR(KSTR))THEN
                          JSTR = JSTR + 1
                          GOTO 1100
                        ELSE
                          ERRSTR = 'bad species string'
                          CALL ERHAND(ERRSTR,IEFATL)
                        ENDIF

                      ENDIF
C                   END OF LOOP 1100
                  ELSE
                    STOSTR = '1  '
                    LENSTO = 1
                  ENDIF

C                 DIAGNOSTICS
                  WRITE(6,'(2I5,2X,A)')JSTR,LENSTO,STOSTR(1:LENSTO)

C                 CONVERT STOICHIOMETRIC COEFFICIENT
                  IF(RSFLAG)THEN
                    VALIDS = CHKSTR(STOSTR(1:LENSTO),LENSTO,ICHRNO)
                    IF(VALIDS)THEN
                      READ(STOSTR(1:LENSTO),*)RSTOIC
                      WRITE(6,'(F6.3)')RSTOIC
                    ELSE
                      ERRSTR = 'bad stoic coeff'
                      CALL ERHAND(ERRSTR,IEFATL)
                    ENDIF
                  ELSE
                    VALIDS = CHKSTR(STOSTR(1:LENSTO),LENSTO,ICHINT)
                    IF(VALIDS)THEN
                      READ(STOSTR(1:LENSTO),*)NSTOIC
                      WRITE(6,'(I6)')NSTOIC
                    ELSE
                      ERRSTR = 'bad stoic coeff'
                      CALL ERHAND(ERRSTR,IEFATL)
                    ENDIF
                  ENDIF

C                 CAPTURE SPECIES STRING
                  LENSPC = LENSTR(KSTR) - JSTR + 1
                  SPCSTR(1:LENSPC) = PARSTR(KSTR)(JSTR:LENSTR(KSTR)) 

                ENDIF

C               CHECK SPECIES STRING IS VALID AND ENUMERATE
                IF(SPCSTR(1:LENSPC).EQ.'=>')THEN
C                 IT'S THE FORWARD ARROW - IGNORE
                ELSE IF(SPCSTR(1:LENSPC).EQ.'==')THEN
C                 IT'S THE EQUILIBRIUM ARROW - IGNORE
                ELSE
                  SFOUND = .FALSE.
                  ISPEC = 1
1200              CONTINUE
              IF(SPCSTR(1:LENSPC).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN
                    SFOUND = .TRUE.

                    IF(RFLAG)THEN
                      IF(RSFLAG)THEN
                      CRTABL(ISPEC,NSTEP) = CRTABL(ISPEC,NSTEP) + RSTOIC
                        RORDER(NSTEP) = RORDER(NSTEP) + RSTOIC
                      ELSE
                      NRTABL(ISPEC,NSTEP) = NRTABL(ISPEC,NSTEP) + NSTOIC
                        NORDER(NSTEP) = NORDER(NSTEP) + NSTOIC
                      ENDIF
                    ELSE
                      IF(RSFLAG)THEN
                    CPTABL(ISPEC,NSTEP) = CPTABL(ISPEC,NSTEP) + RSTOIC
                      ELSE
                    NPTABL(ISPEC,NSTEP) = NPTABL(ISPEC,NSTEP) + NSTOIC
                      ENDIF
                    ENDIF

                  ELSE

                    IF(ISPEC.LT.NSPEC)THEN
                      ISPEC = ISPEC + 1
                      GOTO 1200
                    ELSE
C                     DROP THROUGH TO THIRD BODY CHECK
                      CONTINUE
                    ENDIF

                  ENDIF
C                 END OF LOOP 1200
C                 IF NECESSARY CHECK THIRD BODY STRINGS
                  IF(.NOT.SFOUND)THEN
                    ISPEC = 1
1300                CONTINUE
              IF(SPCSTR(1:LENSPC).EQ.BDYSYM(ISPEC)(1:LENBDY(ISPEC)))THEN
                        MBLIST(NSTEP) = ISPEC
                        IF(RFLAG)NORDER(NSTEP) = NORDER(NSTEP) + 1
                        IF(NSTOIC.NE.1)THEN
                        ERRSTR = 'third body stoic coeff not equal to 1'
                          CALL ERHAND(ERRSTR,IEFATL)
                        ENDIF
                      ELSE
                        IF(ISPEC.LT.NBODY)THEN
                          ISPEC = ISPEC + 1
                          GOTO 1300
                        ELSE
                          ERRSTR = 'species or third body not found'
                          CALL ERHAND(ERRSTR,IEFATL)
                        ENDIF
                      ENDIF
C                   END OF LOOP 1300
                  ENDIF
                ENDIF

              ENDDO

C             DIAGNOSTICS
              WRITE(6,'("Reaction order is:",I5,2F7.2)')
     +                   NORDER(NSTEP),RORDER(NSTEP),
     +                   REAL(NORDER(NSTEP))+RORDER(NSTEP)
            ENDIF
C           END OF CONVERT SPECIES STRINGS

C           CONVERT REMAINING STRINGS TO RATE DATA
C           CHECK FOR VALID NUMBERS
            IF(READLN)THEN
              DO KSTR = ISTR-3,ISTR
                VALIDS = CHKSTR(PARSTR(KSTR),LENSTR(KSTR),ICHRNO)
                IF(.NOT.VALIDS)THEN
                  ERRSTR = 'Lindemann rate data format incorrect'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
              ENDDO
              READ(PARSTR(ISTR-3)(1:LENSTR(ISTR-3)),*)RCLIND(4,NLIND)
              READ(PARSTR(ISTR-2)(1:LENSTR(ISTR-2)),*)RCLIND(1,NLIND)
              READ(PARSTR(ISTR-1)(1:LENSTR(ISTR-1)),*)RCLIND(2,NLIND)
              READ(PARSTR(ISTR)(1:LENSTR(ISTR)),*)RCLIND(3,NLIND)
              READLN = .FALSE.
            ELSEIF(READTR)THEN
              IF(ITREAD.EQ.1)THEN
                DO KSTR = ISTR-6,ISTR
                  VALIDS = CHKSTR(PARSTR(KSTR),LENSTR(KSTR),ICHRNO)
                  IF(.NOT.VALIDS)THEN
                    ERRSTR = 'Troe rate data format incorrect'
                    CALL ERHAND(ERRSTR,IEFATL)
                  ENDIF
                ENDDO
                READ(PARSTR(ISTR-6)(1:LENSTR(ISTR-6)),*)RCTROE(4,NTROE)
                READ(PARSTR(ISTR-5)(1:LENSTR(ISTR-5)),*)RCTROE(5,NTROE)
                READ(PARSTR(ISTR-4)(1:LENSTR(ISTR-4)),*)RCTROE(6,NTROE)
                READ(PARSTR(ISTR-3)(1:LENSTR(ISTR-3)),*)RCTROE(7,NTROE)
                READ(PARSTR(ISTR-2)(1:LENSTR(ISTR-2)),*)RCTROE(1,NTROE)
                READ(PARSTR(ISTR-1)(1:LENSTR(ISTR-1)),*)RCTROE(2,NTROE)
                READ(PARSTR(ISTR)(1:LENSTR(ISTR)),*)RCTROE(3,NTROE)
                ITREAD = 2
              ELSEIF(ITREAD.EQ.2)THEN
                DO KSTR = ISTR-4,ISTR
                  VALIDS = CHKSTR(PARSTR(KSTR),LENSTR(KSTR),ICHRNO)
                  IF(.NOT.VALIDS)THEN
                    ERRSTR = 'Troe rate data format incorrect'
                    CALL ERHAND(ERRSTR,IEFATL)
                  ENDIF
                ENDDO
                READ(PARSTR(ISTR-4)(1:LENSTR(ISTR-4)),*)RCTROE(8,NTROE)
                READ(PARSTR(ISTR-3)(1:LENSTR(ISTR-3)),*)RCTROE(9,NTROE)
                READ(PARSTR(ISTR-2)(1:LENSTR(ISTR-2)),*)RCTROE(10,NTROE)
                READ(PARSTR(ISTR-1)(1:LENSTR(ISTR-1)),*)RCTROE(11,NTROE)
                READ(PARSTR(ISTR)(1:LENSTR(ISTR)),*)RCTROE(12,NTROE)
                READTR = .FALSE.
                ITREAD = 0
              ENDIF
            ELSEIF(READSR)THEN
              IF(ISREAD.EQ.1)THEN
                DO KSTR = ISTR-2,ISTR
                  VALIDS = CHKSTR(PARSTR(KSTR),LENSTR(KSTR),ICHRNO)
                  IF(.NOT.VALIDS)THEN
                    ERRSTR = 'SRI rate data format incorrect'
                    CALL ERHAND(ERRSTR,IEFATL)
                  ENDIF
                ENDDO
                READ(PARSTR(ISTR-2)(1:LENSTR(ISTR-2)),*)RCSRIF(1,NSRIF)
                READ(PARSTR(ISTR-1)(1:LENSTR(ISTR-1)),*)RCSRIF(2,NSRIF)
                READ(PARSTR(ISTR)(1:LENSTR(ISTR)),*)RCSRIF(3,NSRIF)
                ISREAD = 2
              ELSEIF(ISREAD.EQ.2)THEN
                DO KSTR = ISTR-4,ISTR
                  VALIDS = CHKSTR(PARSTR(KSTR),LENSTR(KSTR),ICHRNO)
                  IF(.NOT.VALIDS)THEN
                    ERRSTR = 'SRI rate data format incorrect'
                    CALL ERHAND(ERRSTR,IEFATL)
                  ENDIF
                ENDDO
                READ(PARSTR(ISTR-4)(1:LENSTR(ISTR-4)),*)RCSRIF(4,NSRIF)
                READ(PARSTR(ISTR-3)(1:LENSTR(ISTR-3)),*)RCSRIF(5,NSRIF)
                READ(PARSTR(ISTR-2)(1:LENSTR(ISTR-2)),*)RCSRIF(6,NSRIF)
                READ(PARSTR(ISTR-1)(1:LENSTR(ISTR-1)),*)RCSRIF(7,NSRIF)
                READ(PARSTR(ISTR)(1:LENSTR(ISTR)),*)RCSRIF(8,NSRIF)
                READSR = .FALSE.
                ISREAD = 0
              ENDIF
            ELSE
              DO KSTR = ISTR-2,ISTR
                VALIDS = CHKSTR(PARSTR(KSTR),LENSTR(KSTR),ICHRNO)
                IF(.NOT.VALIDS)THEN
                  ERRSTR = 'Rate data format incorrect'
                  CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
              ENDDO
              READ(PARSTR(ISTR-2)(1:LENSTR(ISTR-2)),*)RPARAM(1,NSTEP)
              READ(PARSTR(ISTR-1)(1:LENSTR(ISTR-1)),*)RPARAM(2,NSTEP)
              READ(PARSTR(ISTR)(1:LENSTR(ISTR)),*)RPARAM(3,NSTEP)
            ENDIF

C           SET FLAG FOR LINDEMANN FORM
            IF(PARSTR(ISTR+1)(1:1).EQ.'L')READLN = .TRUE.
C           SET FLAG AND LINE COUNTER FOR TROE FORM
            IF(PARSTR(ISTR+1)(1:1).EQ.'T')THEN
              READTR = .TRUE.
              ITREAD = 1
            ENDIF
C           SET FLAG AND LINE COUNTER FOR SRI FORM
            IF(PARSTR(ISTR+1)(1:1).EQ.'S')THEN
              READSR = .TRUE.
              ISREAD = 1
            ENDIF
          ENDIF

C         READ ANOTHER LINE
          GOTO 1000

        ENDIF

C     END OF LOOP 1000

C     MECHANISM FILE ERROR HANDLER
      GOTO 1020
1010  ERRSTR = 'unexpected end of mechanism file'
      CALL ERHAND(ERRSTR,IEFATL)
1020  CONTINUE
C     END OF MECHANISM FILE ERROR HANDLER

C     =========================================================================

C     NSTEP CONTAINS THE NUMBER OF STEPS IN THE REACTION MECHANISM

C     NRTABL CONTAINS THE STEP-SPECIES TABLE FOR REACTANTS
C            ENTRY IS 0 IF SPECIES IS NOT INVOLVED IN THIS STEP
C                  OR STOIC COEFF IF SPECIES IS A REACTANT IN THIS STEP

C     NPTABL CONTAINS THE STEP-SPECIES TABLE FOR PRODUCTS
C            ENTRY IS 0 IF SPECIES IS NOT INVOLVED IN THIS STEP
C                  OR STOIC COEFF IF SPECIES IS A PRODUCT IN THIS STEP

C     MBLIST CONTAINS THE LIST OF STEPS INVOLVING THIRD BODIES
C            ENTRY IS 0 IF NO THIRD BODY IN THIS STEP
C                  OR NO. OF THIRD BODY IN THIS STEP

C     EFFY3B CONTAINS THE THIRD BODY EFFICIENCIES FOR EACH THIRD BODY

C     MGLIST CONTAINS THE LIST OF STEPS INVOLVING GIBBS FUNCTION EVALUATION
C            ENTRY IS 0 IF NO GIBBS FUNCTION IN THIS STEP
C            ENTRY IS STEP NO OF THE GIBBS STEP

C     MLLIST CONTAINS THE LIST OF STEPS INVOLVING LINDEMANN RATE PARAMETERS
C            ENTRY IS 0 IF NO LINDEMANN RATE PARAMETERS IN THIS STEP
C            ENTRY IS INDEX NO OF THE LINDEMANN STEP

C     RPARAM CONTAINS THE REACTION RATE PARAMETERS FOR EACH STEP

C     RCLIND CONTAINS THE LINDEMANN RATE PARAMETER DATA
C            FOR EACH LINDEMANN STEP

C     =========================================================================


      RETURN
      END
