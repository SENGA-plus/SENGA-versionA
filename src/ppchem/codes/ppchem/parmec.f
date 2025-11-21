      SUBROUTINE PARMEC(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)
 
C     *************************************************************************
C
C     PARMEC
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     21-SEP-2002:  CREATED
C      
C     DESCRIPTION
C     -----------
C     PARSES A LINE OF DATA IN THE REACTION MECHANISM FOR PPCHEM
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      INTEGER IEFATL
      PARAMETER(IEFATL=5)
      CHARACTER*50 BLANKS
      PARAMETER(
     +   BLANKS='                                                  ')


C     ARGUMENTS
C     =========
      INTEGER NPARSE,NUMSTR,ISTR,ISTOTL
      INTEGER LENSTR(NUMSTR)
      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR


C     LOCAL DATA
C     ==========
      INTEGER IPARSE,ILEN
      CHARACTER*1 PARCHR
      LOGICAL ONWORD,ONSIGN


C     BEGIN
C     =====

C     CLEAR ALL STRINGS
      DO ISTR = 1,NUMSTR
        PARSTR(ISTR) = BLANKS
      ENDDO

      ONWORD = .FALSE.
      ONSIGN = .FALSE.
      IPARSE = 1
      ISTOTL = 0
      ISTR = 0
      ILEN = 0
1000  CONTINUE

C       PICK THE NEXT CHARACTER
        PARCHR = ALINE(IPARSE:IPARSE)

        IF(ONSIGN)THEN
          IF(PARCHR.EQ.'>')THEN
            ISTR = ISTR + 1
            ILEN = 2
            PARSTR(ISTR) = BLANKS
            PARSTR(ISTR)(1:2) = '=>'
            LENSTR(ISTR) = ILEN
            ONSIGN = .FALSE.
          ELSE IF(PARCHR.EQ.'=')THEN
            ISTR = ISTR + 1
            ILEN = 2
            PARSTR(ISTR) = BLANKS
            PARSTR(ISTR)(1:2) = '=='
            LENSTR(ISTR) = ILEN
            PARCHR = '>'
            ONSIGN = .FALSE.
          ELSE
            ERRSTR = 'malformed sign in mechanism list'
            CALL ERHAND(ERRSTR,IEFATL)
          ENDIF
        ENDIF

        IF(ONWORD)THEN

C         CURRENTLY BUILDING A STRING
          IF(PARCHR.EQ.' ')THEN

C           STRING IS FINISHED
            ONWORD = .FALSE.
            LENSTR(ISTR) = ILEN
            ISTOTL = ISTOTL + LENSTR(ISTR)
 
          ELSE IF(PARCHR.EQ.'+')THEN

C           CHECK CONTEXT
            IF((PARSTR(ISTR)(ILEN:ILEN).EQ.'D')
     +      .OR.(PARSTR(ISTR)(ILEN:ILEN).EQ.'E'))THEN

C             RATE DATA FORMAT - STRING IS NOT FINISHED
              ILEN = ILEN + 1
              PARSTR(ISTR)(ILEN:ILEN) = PARCHR

            ELSE

C             STRING IS FINISHED
              ONWORD = .FALSE.
              LENSTR(ISTR) = ILEN
              ISTOTL = ISTOTL + LENSTR(ISTR)

            ENDIF
 
          ELSE IF(PARCHR.EQ.'=')THEN

C           STRING IS FINISHED, BEGUN ARROW
            ONWORD = .FALSE.
            ONSIGN = .TRUE.
            LENSTR(ISTR) = ILEN
            ISTOTL = ISTOTL + LENSTR(ISTR)
 
          ELSE

C           STRING IS NOT FINISHED
            ILEN = ILEN + 1
            PARSTR(ISTR)(ILEN:ILEN) = PARCHR

          ENDIF

        ELSE

C         NOT CURRENTLY BUILDING A STRING
          IF(PARCHR.EQ.' ')THEN

C           STRIP THE BLANK
            CONTINUE

          ELSE IF(PARCHR.EQ.'+')THEN

C           STRIP THE PLUS SIGN
            CONTINUE

          ELSE IF(PARCHR.EQ.'=')THEN

C           BEGUN ARROW
            ONSIGN = .TRUE.

          ELSE IF(PARCHR.EQ.'>')THEN

C           STRIP THE ARROWHEAD
            CONTINUE

          ELSE

C           BEGIN A NEW STRING
            ONWORD = .TRUE.
            ISTR = ISTR + 1
            ILEN = 1
            PARSTR(ISTR) = BLANKS
            PARSTR(ISTR)(1:1) = PARCHR

          ENDIF

        ENDIF

C       PREPARE FOR NEXT CHARACTER
        IPARSE = IPARSE + 1
        IF(IPARSE.LE.NPARSE)THEN
          GOTO 1000
        ELSE
C         END OF LINE
          CONTINUE
        ENDIF
             
C     END OF LOOP 1000


      RETURN
      END
