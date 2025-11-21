      SUBROUTINE PPMBDY
 
C     *************************************************************************
C
C     PPMBDY
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     26-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     PREPROCESSES CHEMICAL MECHANISM THIRD BODY DATA FOR PPCHEM
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
      INTEGER NUMSTR,LENMAX
      PARAMETER(NUMSTR=12, LENMAX=50)


C     EXTERNAL FUNCTION
C     =================
      LOGICAL CHKSTR
      EXTERNAL CHKSTR


C     LOCAL DATA
C     ==========
      INTEGER LENSBC(NBDYMX)
      INTEGER LENSTR(NUMSTR)
      INTEGER ISPEC,JSPEC
      INTEGER IBODY,NSPBC
      INTEGER NPARSE
      INTEGER ISTR,ISTOTL
      INTEGER LENS3B
 
      CHARACTER*10 SBCSYM(NSPCMX)
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR
      CHARACTER*10 SYM3BL

      LOGICAL VALIDS


C     BEGIN
C     =====

C     =========================================================================

C     READ THE THIRD BODY LIST
C     ------------------------
C     READ AND IGNORE HEADER LINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE

C     READ AND PARSE EACH LINE OF THE THIRD BODY LIST
      NBODY = 0
      NSPBC = 0
1100  CONTINUE

C       READ LINE AND INITIALISE PARSER
        READ(NCMECH,'(A)',END=1010)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE

        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)

C       CONDITION FOR END OF THIRD BODY LIST
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
        IF(VALIDS)THEN

C         THIRD BODY LIST FINISHED - DROP THROUGH TO MECHANISM LIST
          CONTINUE

        ELSE

          IF(ISTOTL.EQ.0)THEN

            ERRSTR = 'blank line in third body list'
            CALL ERHAND(ERRSTR,IEWARN)

          ELSE

C           CONVERT STRINGS TO THIRD BODY DATA
            IF((PARSTR(1)(1:1).GE.'0')
     +        .AND.(PARSTR(1)(1:1).LE.'9'))THEN

C             1ST STRING BEGINS WITH A NUMBER => NEW THIRD BODY
C             INCREMENT THIRD BODY COUNTER AND CHECK
              NBODY = NBODY + 1
              IF(NBODY.GT.NBDYMX)THEN
                ERRSTR = 'max no. of third bodies exceeded'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
          
C             CHECK FOR VALID NUMBERS
              VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHINT)
              IF(.NOT.VALIDS)THEN
                ERRSTR = 'third body number incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
              READ(PARSTR(1)(1:LENSTR(1)),*)IBODY
              IF(IBODY.NE.NBODY)THEN
                ERRSTR = 'third body number out of order'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
              NUMBDY(NBODY) = IBODY

C             2ND STRING MUST BE THE THIRD BODY SYMBOL
              LENS3B = LENSTR(2)
              SYM3BL(1:LENS3B) = PARSTR(2)(1:LENS3B)
C             CHECK THIRD BODY SYMBOL IS NOT REPEATED
              IBODY = 1
1200          CONTINUE
              IF(SYM3BL(1:LENS3B).EQ.BDYSYM(IBODY)(1:LENBDY(IBODY)))THEN
                  ERRSTR = 'third body symbol duplicated'
                  CALL ERHAND(ERRSTR,IEFATL)
                ELSE
                  IF(IBODY.LT.NBODY)THEN
                    IBODY = IBODY + 1
                    GOTO 1200
                  ELSE
C                   DROP THROUGH TO SPECIES SYMBOL CHECK
                    CONTINUE
                  ENDIF
                ENDIF
C             END OF LOOP 1200
C             CHECK THIRD BODY SYMBOL IS NOT A SPECIES SYMBOL
              ISPEC = 1
1300          CONTINUE
              IF(SYM3BL(1:LENS3B).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN
                  ERRSTR = 'third body symbol same as a species symbol'
                  CALL ERHAND(ERRSTR,IEFATL)
                ELSE
                  IF(ISPEC.LT.NSPEC)THEN
                    ISPEC = ISPEC + 1
                    GOTO 1300
                  ELSE
                    BDYSYM(NBODY) = SYM3BL(1:LENS3B)
                    LENBDY(NBODY) = LENS3B
                  ENDIF
                ENDIF
C             END OF LOOP 1300

C             CLEAR THE SPECIES CHECKING STRINGS
              DO ISPEC = 1, NSPEC
                SBCSYM(ISPEC) = '          '
                LENSBC(ISPEC) = 0
              ENDDO
              NSPBC = 0

            ELSE

C             1ST STRING MUST BE A SPECIES SYMBOL
              LENS3B = LENSTR(1)
              SYM3BL(1:LENS3B) = PARSTR(1)(1:LENS3B)

C             CHECK SPECIES SYMBOL IS VALID AND ENUMERATE
              ISPEC = 1
1400          CONTINUE
              IF(SYM3BL(1:LENS3B).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN
                  CONTINUE
                ELSE
                  IF(ISPEC.LT.NSPEC)THEN
                    ISPEC = ISPEC + 1
                    GOTO 1400
                  ELSE
                    ERRSTR = 'species not found'
                    CALL ERHAND(ERRSTR,IEFATL)
                  ENDIF
                ENDIF
C             END OF LOOP 1400
              JSPEC = ISPEC

C             CHECK SPECIES SYMBOL IS NOT REPEATED
              ISPEC = 1
1500          CONTINUE
              IF(SYM3BL(1:LENS3B).EQ.SBCSYM(ISPEC)(1:LENSBC(ISPEC)))THEN
                  ERRSTR =
     +                  'species duplicated in 3rd body efficiency list'
                  CALL ERHAND(ERRSTR,IEFATL)
                ELSE
                  IF(ISPEC.LT.NSPBC)THEN
                    ISPEC = ISPEC + 1
                    GOTO 1500
                  ELSE
                    IF(NSPBC.LT.NSPEC)THEN
                      NSPBC = NSPBC + 1
                    ELSE
                      ERRSTR =
     +                  'too many species in 3rd body efficiency list'
                      CALL ERHAND(ERRSTR,IEFATL)
                    ENDIF
                    SBCSYM(NSPBC)(1:LENS3B) = SYM3BL(1:LENS3B)
                    LENSBC(NSPBC) = LENS3B
                  ENDIF
                ENDIF
C             END OF LOOP 1500

C             2ND STRING MUST BE A THIRD BODY EFFICIENCY - CONVERT
C             CHECK FOR VALID NUMBERS
              VALIDS = CHKSTR(PARSTR(2),LENSTR(2),ICHRNO)
              IF(.NOT.VALIDS)THEN
                ERRSTR = 'third body efficiency format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
              READ(PARSTR(2)(1:LENSTR(2)),*)EFFY3B(JSPEC,NBODY)

            ENDIF

          ENDIF

C         READ ANOTHER LINE
          GOTO 1100

        ENDIF

C     END OF LOOP 1100

C     NBODY NOW CONTAINS THE NUMBER OF THIRD BODIES

C     =========================================================================

C     MECHANISM END-OF-FILE ERROR HANDLER
      GOTO 1020
1010  ERRSTR = 'unexpected end of mechanism file'
      CALL ERHAND(ERRSTR,IEFATL)
1020  CONTINUE
C     END OF MECHANISM END-OF-FILE ERROR HANDLER

C     =========================================================================


      RETURN
      END
