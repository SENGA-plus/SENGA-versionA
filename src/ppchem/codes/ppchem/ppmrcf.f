      SUBROUTINE PPMRCF
 
C     *************************************************************************
C
C     PPMRCF
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     01-MAR-2003:  CREATED
C      
C     DESCRIPTION
C     -----------
C     APPLIES CONVERSION FACTORS TO REACTION RATE DATA
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
      DOUBLE PRECISION RFACTR,RFACTS,RFACTT,RFACTM,RFACTV
      DOUBLE PRECISION RFACTA,RFACTE,RFACTL,RFPOWR

      INTEGER LENSTR(NUMSTR)
      INTEGER NPARSE
      INTEGER ISTR,KSTR,ISTOTL
      INTEGER ISTEP

      CHARACTER*80 ALINE
      CHARACTER*50 PARSTR(NUMSTR)
      CHARACTER*50 ERRSTR

      LOGICAL VALIDS


C     BEGIN
C     =====

C     =========================================================================

C     READ THE CONVERSION FACTOR LIST
C     -------------------------------
C     READ AND IGNORE HEADER LINE
      READ(NCMECH,'(A)',END=1010)ALINE
      WRITE(NCREPT,'(A)')ALINE

C     READ AND PARSE THE CONVERSION FACTOR LIST

1000  CONTINUE

C       READ LINE AND INITIALISE PARSER
        READ(NCMECH,'(A)',END=1010)ALINE
        NPARSE = LEN(ALINE)
        WRITE(NCREPT,'(A)')ALINE

        CALL PARSER(ALINE,NPARSE,NUMSTR,LENSTR,PARSTR,ISTR,ISTOTL)

C       CONDITION FOR END OF CONVERSION FACTOR LIST
        VALIDS = CHKSTR(PARSTR(1),LENSTR(1),ICHEND)
        IF(VALIDS)THEN

C         CONVERSION FACTOR LIST FINISHED

        ELSE

          IF(ISTOTL.EQ.0)THEN

            ERRSTR = 'blank line in conversion factor list'
            CALL ERHAND(ERRSTR,IEWARN)

          ELSE

C           PROCESS THE STRINGS

C           CHECK NUMBER OF STRINGS
            IF(ISTR.NE.NFCMAX)THEN
              ERRSTR = 'conversion factor data format incorrect'
              CALL ERHAND(ERRSTR,IEFATL)
            ENDIF

C           CONVERT ALL STRINGS TO CONVERSION FACTOR DATA
C           CHECK FOR VALID NUMBERS
            DO KSTR = 1,ISTR
              VALIDS = CHKSTR(PARSTR(KSTR),LENSTR(KSTR),ICHRNO)
              IF(.NOT.VALIDS)THEN
                ERRSTR = 'conversion factor data format incorrect'
                CALL ERHAND(ERRSTR,IEFATL)
              ENDIF
              READ(PARSTR(KSTR)(1:LENSTR(KSTR)),*)RFCONV(KSTR)
            ENDDO

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

C     APPLY THE CONVERSION FACTORS
C     ----------------------------
C     RFCONV(1): length to metres
C     RFCONV(2): amount of substance to kmol
C     RFCONV(3): time to seconds
C     RFCONV(4): temperature to Kelvin
C     RFCONV(5): Activation energy to Joules

      DO ISTEP = 1,NSTEP

C       CONVERSION FACTOR FOR UNITS OF TIME
        RFACTS = RFCONV(3)

C       CONVERSION FACTOR FOR UNITS OF TEMPERATURE
        RFACTT = EXP(RPARAM(2,ISTEP)*LOG(RFCONV(4)))

C       CONVERSION FACTOR FOR UNITS OF AMOUNT OF SUBSTANCE
        RFACTM = RFCONV(2)

C       CONVERSION FACTOR FOR UNITS OF VOLUME
        RFACTV = EXP(THREE*LOG(RFCONV(1)))
     
C       CONVERSION FACTOR FOR UNITS OF CONCENTRATION
        RFACTR = RFACTM/RFACTV

C       POWER AND CONVERSION FACTOR FOR UNITS OF PRE-EXPONENTIAL FACTOR A
        RFPOWR = REAL(NORDER(ISTEP)-1)
        RFACTL = EXP(-RFPOWR*LOG(RFACTR))/RFACTS
        RFACTA = RFACTL/RFACTT

C       CONVERSION FACTOR FOR UNITS OF ENERGY
        RFACTE = RFCONV(5)
     
C       CONVERSION FACTOR FOR UNITS OF ENERGY PER AMOUNT OF SUBSTANCE
        RFACTE = RFACTE/RFACTM
     
C       CONVERT STANDARD RATE PARAMETERS
        RPARAM(1,ISTEP) = RPARAM(1,ISTEP)*RFACTA
        RPARAM(3,ISTEP) = RPARAM(3,ISTEP)*RFACTE

C       CONVERT LINDEMANN RATE PARAMETERS
        IF(MLLIST(ISTEP).NE.0)THEN

          RPARAM(1,ISTEP) = RPARAM(1,ISTEP)*RFACTR

          RFACTT = EXP(RCLIND(2,ISTEP)*LOG(RFCONV(4)))
          RFACTL = RFACTL/RFACTT
          RCLIND(1,MLLIST(ISTEP)) = RCLIND(1,MLLIST(ISTEP))*RFACTL
          RCLIND(3,MLLIST(ISTEP)) = RCLIND(3,MLLIST(ISTEP))*RFACTE

        ENDIF

C       CONVERT TROE RATE PARAMETERS
        IF(MTLIST(ISTEP).NE.0)THEN

          RPARAM(1,ISTEP) = RPARAM(1,ISTEP)*RFACTR

          RFACTT = EXP(RCTROE(2,ISTEP)*LOG(RFCONV(4)))
          RFACTL = RFACTL/RFACTT
          RCTROE(1,MTLIST(ISTEP)) = RCTROE(1,MTLIST(ISTEP))*RFACTL
          RCTROE(3,MTLIST(ISTEP)) = RCTROE(3,MTLIST(ISTEP))*RFACTE

        ENDIF

C       CONVERT SRI RATE PARAMETERS
        IF(MSLIST(ISTEP).NE.0)THEN

          RPARAM(1,ISTEP) = RPARAM(1,ISTEP)*RFACTR

          RFACTT = EXP(RCSRIF(2,ISTEP)*LOG(RFCONV(4)))
          RFACTL = RFACTL/RFACTT
          RCSRIF(1,MSLIST(ISTEP)) = RCSRIF(1,MSLIST(ISTEP))*RFACTL
          RCSRIF(3,MSLIST(ISTEP)) = RCSRIF(3,MSLIST(ISTEP))*RFACTE

        ENDIF

      ENDDO

C     =========================================================================


      RETURN
      END
