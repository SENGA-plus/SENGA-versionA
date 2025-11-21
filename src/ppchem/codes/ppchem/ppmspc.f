      SUBROUTINE PPMSPC
 
C     *************************************************************************
C
C     PPMSPC

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
C     PREPROCESSES CHEMICAL MECHANISM SPECIES DATA FOR PPCHEM
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
      INTEGER ISTR,JSTR,ISTOTL
      INTEGER NAMLEN
      INTEGER LENCHK
 
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR
      CHARACTER*10 SYMCHK

      LOGICAL VALIDS


C     BEGIN
C     =====

C     =========================================================================

C     READ THE SPECIES LIST
C     ---------------------
C     READ AND IGNORE HEADER LINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE

C     READ AND PARSE EACH LINE OF THE SPECIES LIST
      NSPEC = 0
1000  CONTINUE

C       INCREMENT SPECIES COUNTER AND CHECK
        NSPEC = NSPEC + 1
        IF(NSPEC.GT.NSPCMX)THEN
          ERRSTR = 'max no. of species exceeded'
          CALL ERHAND(ERRSTR,IEFATL)
        ENDIF

C       READ LINE AND INITIALISE PARSER
        READ(NCMECH,'(A)',END=1010)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE

        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)

C       CONDITION FOR END OF SPECIES LIST
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
        IF(VALIDS)THEN

C         SPECIES LIST FINISHED - DROP THROUGH TO MECHANISM LIST
          NSPEC = NSPEC - 1

        ELSE

          IF(ISTOTL.EQ.0)THEN

            ERRSTR = 'blank line in species list'
            CALL ERHAND(ERRSTR,IEWARN)

          ELSE

C           CONVERT STRINGS TO SPECIES DATA

C           SPECIES NUMBER
C           CHECK FOR VALID NUMBERS
            DO JSTR = 1, LENSTR(1)
              IF((PARSTR(1)(JSTR:JSTR).LT.'0')
     +          .OR.(PARSTR(1)(JSTR:JSTR).GT.'9'))THEN
                ERRSTR = 'species number format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
            ENDDO
            READ(PARSTR(1)(1:LENSTR(1)),*)ISPEC
            IF(ISPEC.NE.NSPEC)THEN
              ERRSTR = 'species number out of order'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF
            NUMSPC(NSPEC) = ISPEC

C           SPECIES SYMBOL
            LENCHK = LENSTR(2)
            SYMCHK(1:LENCHK) = PARSTR(2)(1:LENCHK)

C           CHECK SPECIES SYMBOL DOES NOT END WITH "D" OR "E"
C           REQUIRED TO AVOID CONFLICT WITH RATE DATA FORMAT
            IF((SYMCHK(LENCHK:LENCHK).EQ.'D')
     +      .OR.(SYMCHK(LENCHK:LENCHK).EQ.'E'))THEN
              ERRSTR = 'species name cannot end with D or E'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF

C           CHECK SPECIES SYMBOL IS NOT REPEATED
            ISPEC = 1
1100        CONTINUE
              IF(SYMCHK(1:LENCHK).EQ.SPCSYM(ISPEC)(1:LENSYM(ISPEC)))THEN
                ERRSTR = 'species symbol duplicated'
                CALL ERHAND(ERRSTR,IEFATL)
              ELSE
                IF(ISPEC.LT.NSPEC)THEN
                  ISPEC = ISPEC + 1
                  GOTO 1100
                ELSE
                  SPCSYM(NSPEC) = SYMCHK(1:LENCHK)
                  LENSYM(NSPEC) = LENCHK
                ENDIF
              ENDIF
C           END OF LOOP 1100

C           SPECIES NAME
            SPCNAM(NSPEC) = PARSTR(3)(1:LENSTR(3))
            NAMLEN = LENSTR(3)
            DO JSTR = 4,ISTR
              SPCNAM(NSPEC) = SPCNAM(NSPEC)(1:NAMLEN)
     +          //' '//PARSTR(JSTR)(1:LENSTR(JSTR))
              NAMLEN = NAMLEN + 1 + LENSTR(ISTR)
            ENDDO
            LENNAM(NSPEC) = NAMLEN
 
          ENDIF

C         READ ANOTHER LINE
          GOTO 1000

        ENDIF

C     END OF LOOP 1000

C     NSPEC NOW CONTAINS THE NUMBER OF SPECIES

C     =========================================================================

C     END-OF-FILE ERROR HANDLER
      GOTO 1020
1010  ERRSTR = 'unexpected end of mechanism file'
      CALL ERHAND(ERRSTR,IEFATL)
1020  CONTINUE
C     END OF END-OF-FILE ERROR HANDLER

C     =========================================================================


      RETURN
      END
