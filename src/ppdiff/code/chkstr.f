      FUNCTION CHKSTR(STRING,LENSTR,ICHECK)
 
C     *************************************************************************
C
C     CHKSTR
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
C     CHECKS FORMAT OF CHARACTER STRING IN PPTHRM
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      INTEGER ICHNUL,ICHARS,ICHEND,ICHINT,ICHSNT,ICHRNO,ICHFNO
      PARAMETER(ICHNUL=0,
     +          ICHARS=1, ICHEND=2, ICHINT=3, ICHSNT=4,
     +          ICHRNO=5, ICHFNO=6)


C     FUNCTION
C     ========
      LOGICAL CHKSTR


C     ARGUMENTS
C     =========
      INTEGER LENSTR,ICHECK
      CHARACTER*50 STRING


C     LOCAL DATA
C     ==========
      INTEGER ISTR
      LOGICAL FLGSTR


C     BEGIN
C     =====
      CHKSTR = .FALSE.

C     NULL CASE
      IF(ICHECK.EQ.ICHNUL)THEN
        CHKSTR = .TRUE.
      ENDIF

C     CHECK STRING IS STRICTLY A CHARACTER STRING
      IF(ICHECK.EQ.ICHARS)THEN
        ISTR = 0
1000    CONTINUE
          ISTR = ISTR + 1
          FLGSTR
     +     = ((STRING(ISTR:ISTR).GE.'A').AND.(STRING(ISTR:ISTR).LE.'Z'))
          FLGSTR = FLGSTR
     +   .OR.((STRING(ISTR:ISTR).GE.'a').AND.(STRING(ISTR:ISTR).LE.'z'))
          FLGSTR = FLGSTR
     +   .OR.((STRING(ISTR:ISTR).GE.'0').AND.(STRING(ISTR:ISTR).LE.'9'))
          IF(FLGSTR)THEN
            IF(ISTR.LT.LENSTR)THEN
              GOTO 1000
            ELSE
              CHKSTR = .TRUE.
            ENDIF
          ENDIF
C       END OF LOOP 1000
      ENDIF
      
C     CHECK FOR END
      IF(ICHECK.EQ.ICHEND)THEN
        IF((STRING(1:LENSTR).EQ.'END')      
     +    .OR.(STRING(1:LENSTR).EQ.'end')      
     +      .OR.(STRING(1:LENSTR).EQ.'End'))CHKSTR = .TRUE.
      ENDIF

C     CHECK STRING REPRESENTS AN UNSIGNED INTEGER
      IF(ICHECK.EQ.ICHINT)THEN
        ISTR = 0
3000    CONTINUE
          ISTR = ISTR + 1
          FLGSTR 
     +     = ((STRING(ISTR:ISTR).GE.'0').AND.(STRING(ISTR:ISTR).LE.'9'))
          IF(FLGSTR)THEN
            IF(ISTR.LT.LENSTR)THEN
              GOTO 3000
            ELSE
              CHKSTR = .TRUE.
            ENDIF
          ENDIF
C       END OF LOOP 3000
      ENDIF

C     CHECK STRING REPRESENTS A SIGNED INTEGER
      IF(ICHECK.EQ.ICHSNT)THEN
        ISTR = 0
4000    CONTINUE
          ISTR = ISTR + 1
          FLGSTR
     +     = ((STRING(ISTR:ISTR).GE.'0').AND.(STRING(ISTR:ISTR).LE.'9'))
          IF(ISTR.EQ.1)FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'-')
          IF(FLGSTR)THEN
            IF(ISTR.LT.LENSTR)THEN
              GOTO 4000
            ELSE
              CHKSTR = .TRUE.
            ENDIF
          ENDIF
C       END OF LOOP 4000
      ENDIF

C     CHECK STRING COULD REPRESENT A REAL NUMBER
C     IN F,D OR E FORMAT
      IF(ICHECK.EQ.ICHRNO)THEN
        ISTR = 0
5000    CONTINUE
          ISTR = ISTR + 1
          FLGSTR 
     +     = ((STRING(ISTR:ISTR).GE.'0').AND.(STRING(ISTR:ISTR).LE.'9'))
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'-')
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'+')
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'.')
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'D')
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'E')
          IF(FLGSTR)THEN
            IF(ISTR.LT.LENSTR)THEN
              GOTO 5000
            ELSE
              CHKSTR = .TRUE.
            ENDIF
          ENDIF
C       END OF LOOP 5000
      ENDIF

C     CHECK STRING COULD REPRESENT A REAL NUMBER
C     IN F FORMAT ONLY
      IF(ICHECK.EQ.ICHFNO)THEN
        ISTR = 0
6000    CONTINUE
          ISTR = ISTR + 1
          FLGSTR 
     +     = ((STRING(ISTR:ISTR).GE.'0').AND.(STRING(ISTR:ISTR).LE.'9'))
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'-')
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'+')
          FLGSTR = FLGSTR.OR.(STRING(ISTR:ISTR).EQ.'.')
          IF(FLGSTR)THEN
            IF(ISTR.LT.LENSTR)THEN
              GOTO 6000
            ELSE
              CHKSTR = .TRUE.
            ENDIF
          ENDIF
C       END OF LOOP 6000
      ENDIF


      RETURN
      END
