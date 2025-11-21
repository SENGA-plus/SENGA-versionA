      SUBROUTINE PPSPEC
 
C     *************************************************************************
C
C     PPSPEC
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
C     PREPROCESSES SPECIES DATA FOR PPCHEM
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
C     ==================
      LOGICAL CHKSTR
      EXTERNAL CHKSTR


C     LOCAL DATA
C     ==========
      INTEGER LENSTR(NUMSTR)
      INTEGER ISPEC
      INTEGER NPARSE
      INTEGER ISTR,ISTOTL
      INTEGER IC,JC
      INTEGER ITINT
      INTEGER LENCHK
 
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR
      CHARACTER*10 SYMCHK

      LOGICAL SPCFLG(NSPCMX)
      LOGICAL VALIDS


C     BEGIN
C     =====

C     =========================================================================

C     READ THE SPECIES DATA
C     ---------------------
      OPEN(UNIT=NCSPEC,FILE=FNSPEC,STATUS='OLD',FORM='FORMATTED',
     +     ERR=1010)

C     FIVE LINES OF COMMENT ARE IGNORED
C     ---------------------------------
      READ(NCSPEC,'(A)',END=1020)ALINE
      READ(NCSPEC,'(A)',END=1020)ALINE
      READ(NCSPEC,'(A)',END=1020)ALINE
      READ(NCSPEC,'(A)',END=1020)ALINE
      READ(NCSPEC,'(A)',END=1020)ALINE

C     READ THE REFERENCE PRESSURE
      READ(NCSPEC,*)PREFGB

C     SPECIES-DETECTED FLAGS
      DO ISPEC = 1, NSPEC
        SPCFLG(ISPEC) = .FALSE.
      ENDDO

1000  CONTINUE

C       READ LINE AND INITIALISE PARSER
1100    CONTINUE
          READ(NCSPEC,'(A)',END=1020)ALINE
          NPARSE = LEN(ALINE)
          WRITE(NCREPT,'(A)')ALINE
          CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
          IF(ISTOTL.EQ.0)THEN
            ERRSTR = 'blank line in species data'
            CALL ERHAND(ERRSTR,IEWARN)
            GOTO 1100
          ENDIF
C       END OF LOOP 1100

C       CHECK FOR END OF FILE
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
        IF(VALIDS)THEN

C         END OF FILE DETECTED
C         CHECK THAT DATA HAS BEEN COLLECTED FOR ALL SPECIES
          ISPEC = 1
1200      CONTINUE
            IF(SPCFLG(ISPEC))THEN
              IF(ISPEC.LT.NSPEC)THEN
                ISPEC = ISPEC + 1
                GOTO 1200
              ENDIF
            ELSE
              ERRSTR = 'end of species file: species data incomplete'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
C         END OF LOOP 1200

        ELSE

C         CHECK SPECIES SYMBOL
          LENCHK = LENSTR(1)
          SYMCHK(1:LENCHK) = PARSTR(1)(1:LENCHK)

          ISPEC = 1
1300      CONTINUE
            IF(SYMCHK(1:LENCHK).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN

C             SPECIES SYMBOL RECOGNISED

C             SET SPECIES FLAG
              SPCFLG(ISPEC) = .TRUE.

C             PROCEED TO READ THE SPECIES DATA
              ERRSTR = 'invalid species data format'

C             LINE 1: WMOLAR
              READ(NCSPEC,'(A)',END=1020)ALINE
              NPARSE = LEN(ALINE)
              WRITE(NCREPT,'(A)')ALINE
              CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
              IF(ISTR.LT.1)THEN
                CALL ERHAND(ERRSTR,IEFATL)
              ELSE
                VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHRNO)
                IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
              READ(PARSTR(1),*)WMOLAR(ISPEC)

C             LINE 2: CLEWIS
              READ(NCSPEC,'(A)',END=1020)ALINE
              NPARSE = LEN(ALINE)
              WRITE(NCREPT,'(A)')ALINE
              CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
              IF(ISTR.LT.1)THEN
                CALL ERHAND(ERRSTR,IEFATL)
              ELSE
                VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHRNO)
                IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
              READ(PARSTR(1),*)CLEWIS(ISPEC)

C             LINE 3: NTINT
              READ(NCSPEC,'(A)',END=1020)ALINE
              NPARSE = LEN(ALINE)
              WRITE(NCREPT,'(A)')ALINE
              CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
              IF(ISTR.LT.1)THEN
                CALL ERHAND(ERRSTR,IEFATL)
              ELSE
                VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHINT)
                IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
              READ(PARSTR(1),*)NTINT(ISPEC)

C             LINE 4 TO END-OF-SPECIES
              DO IC = 1,NTINT(ISPEC)
                READ(NCSPEC,'(A)',END=1020)ALINE
                NPARSE = LEN(ALINE)
                WRITE(NCREPT,'(A)')ALINE
              CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
                IF(ISTR.LT.4)THEN
                  CALL ERHAND(ERRSTR,IEFATL)
                ELSE
                  VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHINT)
                  IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
                  VALIDS = CHKSTR(PARSTR(2),LENSTR(2),ICHRNO)
                  IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
                  VALIDS = CHKSTR(PARSTR(3),LENSTR(3),ICHRNO)
                  IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
                  VALIDS = CHKSTR(PARSTR(4),LENSTR(4),ICHINT)
                  IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
                ENDIF
                READ(PARSTR(1),*)ITINT
                READ(PARSTR(2),*)TINTLO(ITINT,ISPEC)
                READ(PARSTR(3),*)TINTHI(ITINT,ISPEC)
                READ(PARSTR(4),*)NCOFCP(ITINT,ISPEC)
                DO JC = 1,NCOFCP(ITINT,ISPEC)
                  READ(NCSPEC,'(A)',END=1020)ALINE
                  NPARSE = LEN(ALINE)
                  WRITE(NCREPT,'(A)')ALINE
              CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
                  IF(ISTR.LT.1)THEN
                    CALL ERHAND(ERRSTR,IEFATL)
                  ELSE
                    VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHRNO)
                    IF(.NOT.VALIDS)CALL ERHAND(ERRSTR,IEFATL)
                  ENDIF
                  READ(PARSTR(1),*)ACOFCP(JC,IC,ISPEC)
                ENDDO
              ENDDO

C             END-OF-SPECIES LINE
              READ(NCSPEC,'(A)',END=1020)ALINE
              NPARSE = LEN(ALINE)
              WRITE(NCREPT,'(A)')ALINE
              CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
              VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
              IF(.NOT.VALIDS)THEN
                ERRSTR = 'species file: missing END'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF

            ELSE

C             SPECIES SYMBOL NOT RECOGNISED

              IF(ISPEC.LT.NSPEC)THEN

C               CHECK THE NEXT SPECIES SYMBOL
                ISPEC = ISPEC + 1
                GOTO 1300

              ELSE

C               ALL SPECIES CHECKED - NO MATCH
C               READ AND IGNORE THE SPECIES DATA

C               READ ON UNTIL END OF SPECIES DATA BLOCK
1400            CONTINUE
                  READ(NCSPEC,'(A)',END=1020)ALINE
                  NPARSE = LEN(ALINE)
                  WRITE(NCREPT,'(A)')ALINE
              CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
                  VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
                  IF(.NOT.VALIDS)GOTO 1400
C               END OF LOOP 1400

              ENDIF

C           END OF LOOP 1300

          ENDIF
C         CHECK SPECIES SYMBOL

C         READ ANOTHER LINE
          GOTO 1000

        ENDIF
C       CHECK FOR END OF FILE

C     END OF LOOP 1000

C     =========================================================================

C     CLOSE THE SPECIES FILE
C     ----------------------
      CLOSE(NCSPEC)

C     =========================================================================

C     SPECIES FILE ERROR HANDLER
      GOTO 1030
1010  ERRSTR = 'Error opening species file:'
      CALL ERHAND(ERRSTR,IEFATL)
1020  ERRSTR = 'unexpected end of species data file'
      CALL ERHAND(ERRSTR,IEFATL)
1030  CONTINUE
C     END OF SPECIES FILE ERROR HANDLER

C     =========================================================================


      RETURN
      END
