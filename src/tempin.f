      SUBROUTINE TEMPIN
 
C     *************************************************************************
C
C     TEMPIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     15-MAY-2004:  CREATED
C     
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES TEMPERATURE AND PRESSURE
C     USES BISECTION METHOD FOR ROBUSTNESS
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
      DOUBLE PRECISION TOLTMP
      PARAMETER(TOLTMP = 1.0D-4)
      DOUBLE PRECISION TININC
      PARAMETER(TININC = 5.0D1)
      DOUBLE PRECISION TLIMLO,TLIMHI
      PARAMETER(TLIMLO = 2.0D2, TLIMHI = 3.0D3)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TCOEFF(0:NCTMAX)
      DOUBLE PRECISION UKUK
      DOUBLE PRECISION TEMPOR,TUPPER,TLOWER,TRESID,TUK2ME,CPFORY
      INTEGER IC,JC,KC,ISPEC,ITINT,ICP
      INTEGER IINDEX,IPOWER
      LOGICAL FNCONV 


C     BEGIN
C     =====

C     =========================================================================

C     TEMPERATURE AND PRESSURE
C     ------------------------

C     TEMPERATURE AND PRESSURE ARE PARALLEL
      DO KC = KSTALT,KSTOLT
        DO JC = JSTALT,JSTOLT
          DO IC = ISTALT,ISTOLT

C           ===================================================================

C           KINETIC ENERGY TERM
            UKUK = (URHS(IC,JC,KC)*URHS(IC,JC,KC)
     +           +  VRHS(IC,JC,KC)*VRHS(IC,JC,KC)
     +           +  WRHS(IC,JC,KC)*WRHS(IC,JC,KC))/DRHS(IC,JC,KC)

            TUK2ME = HALF*UKUK - ERHS(IC,JC,KC)

C           ===================================================================

C           SOLVE FOR TEMPERATURE
C           =====================
C           USING BISECTION
            TLOWER = TLIMLO
            TUPPER = TLIMHI

C           SET NON-CONVERGENCE FLAG
            FNCONV = .TRUE.

C           SET TEMPERATURE FROM INITIAL GUESS
            TEMPOR = TRIN

C           ===================================================================

C           FIRST BRACKET THE ROOT
C           ======================

C           -------------------------------------------------------------------

C           INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
            TCOEFF(0) = TUK2ME
            DO ICP = 1, NCTMAX
              TCOEFF(ICP) = ZERO
            ENDDO

C           FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
            DO ISPEC = 1,NSPEC

              ITINT = 1
1100          CONTINUE
                IF(TEMPOR.GT.TINTHI(ITINT,ISPEC))THEN
                  IF(ITINT.LT.NTINT(ISPEC))THEN
                    ITINT = ITINT + 1
                    GOTO 1100
                  ENDIF
                ENDIF
C             END OF LOOP 1100

C             CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
              TCOEFF(0) = TCOEFF(0)
     +                  + YRHS(IC,JC,KC,ISPEC)*
     +                    AMASCP(NCENTH(ITINT,ISPEC),ITINT,ISPEC)

              TCOEFF(1) = TCOEFF(1)
     +                  + YRHS(IC,JC,KC,ISPEC)*AMASCT(1,ITINT,ISPEC)
              DO ICP = 2, NCPOLY(ITINT,ISPEC)
                TCOEFF(ICP) = TCOEFF(ICP)
     +                  + YRHS(IC,JC,KC,ISPEC)*AMASCT(ICP,ITINT,ISPEC)
              ENDDO

            ENDDO
C           END OF RUN THROUGH ALL SPECIES

C           -------------------------------------------------------------------

C           EVALUATE TEMPERATURE RESIDUAL
            TRESID = TCOEFF(NCTMAX)
            DO ICP = NCTMM1,1,-1
              TRESID = TCOEFF(ICP) + TRESID*TEMPOR
            ENDDO
            TRESID = TCOEFF(0) + TRESID*TEMPOR

C           -------------------------------------------------------------------

C           CHECK INITIAL GUESS FOR CONVERGENCE
            IF(ABS(TRESID).LT.TOLTMP)THEN

C             -----------------------------------------------------------------

C             CONVERGED ON FIRST PASS
              FNCONV = .FALSE.

C             -----------------------------------------------------------------

            ELSEIF(TRESID.LT.ZERO)THEN

C             -----------------------------------------------------------------

C             INITIAL GUESS IS TOO LOW
2000          CONTINUE

C               SET INITIAL GUESS AS LOWER LIMIT AND TRY AGAIN
                TLOWER = TEMPOR
                TEMPOR = TEMPOR + TININC

C               ---------------------------------------------------------------

C               INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                TCOEFF(0) = TUK2ME
                DO ICP = 1, NCTMAX
                  TCOEFF(ICP) = ZERO
                ENDDO

C               FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
                DO ISPEC = 1,NSPEC

                  ITINT = 1
2100              CONTINUE
                    IF(TEMPOR.GT.TINTHI(ITINT,ISPEC))THEN
                      IF(ITINT.LT.NTINT(ISPEC))THEN
                        ITINT = ITINT + 1
                        GOTO 2100
                      ENDIF
                    ENDIF
C                 END OF LOOP 2100

C                 CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                  TCOEFF(0) = TCOEFF(0)
     +                      + YRHS(IC,JC,KC,ISPEC)*
     +                        AMASCP(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
                  TCOEFF(1) = TCOEFF(1)
     +                      + YRHS(IC,JC,KC,ISPEC)*AMASCT(1,ITINT,ISPEC)
                  DO ICP = 2, NCPOLY(ITINT,ISPEC)
                    TCOEFF(ICP) = TCOEFF(ICP)
     +                    + YRHS(IC,JC,KC,ISPEC)*AMASCT(ICP,ITINT,ISPEC)
                  ENDDO

                ENDDO
C               END OF RUN THROUGH ALL SPECIES

C               ---------------------------------------------------------------

C               EVALUATE TEMPERATURE RESIDUAL
                TRESID = TCOEFF(NCTMAX)
                DO ICP = NCTMM1,1,-1
                  TRESID = TCOEFF(ICP) + TRESID*TEMPOR
                ENDDO
                TRESID = TCOEFF(0) + TRESID*TEMPOR

C               ---------------------------------------------------------------

C               CHECK NEW GUESS FOR CONVERGENCE
                IF(ABS(TRESID).LT.TOLTMP)THEN

C                 -------------------------------------------------------------

C                 NEW GUESS HAS CONVERGED
                  FNCONV = .FALSE.

C                 -------------------------------------------------------------

                ELSEIF(TRESID.LT.ZERO)THEN

C                 -------------------------------------------------------------

C                 NEW GUESS IS STILL TOO LOW: GO ROUND AGAIN
                  IF(TEMPOR.LT.TLIMHI)THEN
                    GO TO 2000
                  ELSE
                    WRITE(6,*)
     +              'Fatal: TEMPIN: T upper bracket failed to converge'
                    WRITE(6,*)'processor:',IPROC
                    WRITE(6,*)'at point:',IC,JC,KC
                    WRITE(6,*)'with values:',TEMPOR,TRESID
                    WRITE(6,*)DRHS(IC,JC,KC)
                    WRITE(6,*)URHS(IC,JC,KC)
                    WRITE(6,*)VRHS(IC,JC,KC)
                    WRITE(6,*)WRHS(IC,JC,KC)
                    WRITE(6,*)ERHS(IC,JC,KC)
                    DO ISPEC = 1, NSPEC
                      WRITE(6,*)YRHS(IC,JC,KC,ISPEC)
                    ENDDO
                    STOP
                  ENDIF

C                 -------------------------------------------------------------

                ELSEIF(TRESID.GT.ZERO)THEN

C                 -------------------------------------------------------------

C                 ROOT IS BRACKETED
                  TUPPER = TEMPOR

C                 -------------------------------------------------------------

                ENDIF

C               ---------------------------------------------------------------

C             END OF LOOP 2000

C             -----------------------------------------------------------------

            ELSEIF(TRESID.GT.ZERO)THEN

C             -----------------------------------------------------------------

C             INITIAL GUESS IS TOO HIGH
3000          CONTINUE

C               SET INITIAL GUESS AS UPPER LIMIT AND TRY AGAIN
                TUPPER = TEMPOR
                TEMPOR = TEMPOR - TININC

C               ---------------------------------------------------------------

C               INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                TCOEFF(0) = TUK2ME
                DO ICP = 1, NCTMAX
                  TCOEFF(ICP) = ZERO
                ENDDO

C               FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
                DO ISPEC = 1,NSPEC

                  ITINT = 1
3100              CONTINUE
                    IF(TEMPOR.GT.TINTHI(ITINT,ISPEC))THEN
                      IF(ITINT.LT.NTINT(ISPEC))THEN
                        ITINT = ITINT + 1
                        GOTO 3100
                      ENDIF
                    ENDIF
C                 END OF LOOP 3100

C                 CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                  TCOEFF(0) = TCOEFF(0)
     +                      + YRHS(IC,JC,KC,ISPEC)*
     +                        AMASCP(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
                  TCOEFF(1) = TCOEFF(1)
     +                      + YRHS(IC,JC,KC,ISPEC)*AMASCT(1,ITINT,ISPEC)
                  DO ICP = 2, NCPOLY(ITINT,ISPEC)
                    TCOEFF(ICP) = TCOEFF(ICP)
     +                    + YRHS(IC,JC,KC,ISPEC)*AMASCT(ICP,ITINT,ISPEC)
                  ENDDO

                ENDDO
C               END OF RUN THROUGH ALL SPECIES

C               ---------------------------------------------------------------

C               EVALUATE TEMPERATURE RESIDUAL
                TRESID = TCOEFF(NCTMAX)
                DO ICP = NCTMM1,1,-1
                  TRESID = TCOEFF(ICP) + TRESID*TEMPOR
                ENDDO
                TRESID = TCOEFF(0) + TRESID*TEMPOR

C               ---------------------------------------------------------------

C               CHECK NEW GUESS FOR CONVERGENCE
                IF(ABS(TRESID).LT.TOLTMP)THEN

C                 -------------------------------------------------------------

C                 NEW GUESS HAS CONVERGED
                  FNCONV = .FALSE.

C                 -------------------------------------------------------------

                ELSEIF(TRESID.GT.ZERO)THEN

C                 -------------------------------------------------------------

C                 NEW GUESS IS STILL TOO HIGH: GO ROUND AGAIN
                  IF(TEMPOR.GT.TLIMLO)THEN
                    GO TO 3000
                  ELSE
                    WRITE(6,*)
     +              'Fatal: TEMPIN: T lower bracket failed to converge'
                    WRITE(6,*)'processor:',IPROC
                    WRITE(6,*)'at point:',IC,JC,KC
                    WRITE(6,*)'with values:',TEMPOR,TRESID
                    WRITE(6,*)DRHS(IC,JC,KC)
                    WRITE(6,*)URHS(IC,JC,KC)
                    WRITE(6,*)VRHS(IC,JC,KC)
                    WRITE(6,*)WRHS(IC,JC,KC)
                    WRITE(6,*)ERHS(IC,JC,KC)
                    DO ISPEC = 1, NSPEC
                      WRITE(6,*)YRHS(IC,JC,KC,ISPEC)
                    ENDDO
                    STOP
                  ENDIF

C                 -------------------------------------------------------------

                ELSEIF(TRESID.LT.ZERO)THEN

C                 -------------------------------------------------------------

C                 ROOT IS BRACKETED
                  TLOWER = TEMPOR

C                 -------------------------------------------------------------

                ENDIF

C               ---------------------------------------------------------------

C             END OF LOOP 3000

C             -----------------------------------------------------------------

            ENDIF
C           END OF CHECK INITIAL GUESS FOR CONVERGENCE

C           ===================================================================

C           ROOT IS BRACKETED
C           =================
C           NOW USE BISECTION TO REFINE THE ROOT

            IF(FNCONV)THEN

4000          CONTINUE

C               BISECT
                TEMPOR = HALF*(TLOWER+TUPPER)

C               ---------------------------------------------------------------

C               INITIALISE COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                TCOEFF(0) = TUK2ME
                DO ICP = 1, NCTMAX
                  TCOEFF(ICP) = ZERO
                ENDDO

C               FOR EACH SPECIES LOCATE TEMPERATURE IN AN INTERVAL
                DO ISPEC = 1,NSPEC

                  ITINT = 1
4100              CONTINUE
                    IF(TEMPOR.GT.TINTHI(ITINT,ISPEC))THEN
                      IF(ITINT.LT.NTINT(ISPEC))THEN
                        ITINT = ITINT + 1
                        GOTO 4100
                      ENDIF
                    ENDIF
C                 END OF LOOP 4100

C                 CONSTRUCT COEFFICIENTS OF TEMPERATURE POLYNOMIAL
                  TCOEFF(0) = TCOEFF(0)
     +                      + YRHS(IC,JC,KC,ISPEC)*
     +                        AMASCP(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
                  TCOEFF(1) = TCOEFF(1)
     +                      + YRHS(IC,JC,KC,ISPEC)*AMASCT(1,ITINT,ISPEC)
                  DO ICP = 2, NCPOLY(ITINT,ISPEC)
                    TCOEFF(ICP) = TCOEFF(ICP)
     +                    + YRHS(IC,JC,KC,ISPEC)*AMASCT(ICP,ITINT,ISPEC)
                  ENDDO

                ENDDO
C               END OF RUN THROUGH ALL SPECIES

C               ---------------------------------------------------------------

C               EVALUATE TEMPERATURE RESIDUAL
                TRESID = TCOEFF(NCTMAX)
                DO ICP = NCTMM1,1,-1
                  TRESID = TCOEFF(ICP) + TRESID*TEMPOR
                ENDDO
                TRESID = TCOEFF(0) + TRESID*TEMPOR

C               ---------------------------------------------------------------

                IF(ABS(TRESID).LT.TOLTMP)THEN

C                 CONVERGED
                  TRUN(IC,JC,KC) = TEMPOR

                ELSEIF(TRESID.LT.ZERO)THEN

                  TLOWER = TEMPOR
                  GOTO 4000

                ELSEIF(TRESID.GT.ZERO)THEN

                  TUPPER = TEMPOR
                  GOTO 4000

                ENDIF

C             -----------------------------------------------------------------
C             END OF LOOP 4000

            ENDIF
C           END OF BISECTION

C           ===================================================================

C           SET THE NEW TEMPERATURE
            TRUN(IC,JC,KC) = TEMPOR

C           ===================================================================

C           CONSTRUCT THE TEMPERATURE INTERVAL INDEX
C           EVALUATE PRESSURE
C           EVALUATE MIXTURE SPECIFIC HEAT CP
            DO IINDEX = 1,NINTMX
              ITNDEX(IC,JC,KC,IINDEX) = 0
            ENDDO
            STORE7(IC,JC,KC) = ZERO
            TRANSP(IC,JC,KC) = ZERO
            DO ISPEC = 1,NSPEC

              ITINT = 1
5100          CONTINUE
                IF(TRUN(IC,JC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                  IF(ITINT.LT.NTINT(ISPEC))THEN
                    ITINT = ITINT + 1
                    GOTO 5100
                  ENDIF
                ENDIF
C             END OF LOOP 5100

C             SET THE TEMPERATURE INDEX
              IINDEX = 1 + (ISPEC-1)/NSPIMX
              IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
              ITNDEX(IC,JC,KC,IINDEX) = ITNDEX(IC,JC,KC,IINDEX)
     +                                + (ITINT-1)*NTBASE**IPOWER

C             =================================================================

C             EVALUATE MIXTURE SPECIFIC HEAT CP
              CPFORY = AMASCP(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
              DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                CPFORY = CPFORY*TRUN(IC,JC,KC) + AMASCP(ICP,ITINT,ISPEC)
              ENDDO
              TRANSP(IC,JC,KC) = TRANSP(IC,JC,KC)
     +                         + YRHS(IC,JC,KC,ISPEC)*CPFORY

C             =================================================================

C             EVALUATE (DENSITY TIMES) MIXTURE GAS CONSTANT FOR PRESSURE
              STORE7(IC,JC,KC) = STORE7(IC,JC,KC)
     +                         + YRHS(IC,JC,KC,ISPEC)*RGSPEC(ISPEC)

C             =================================================================

            ENDDO
            TRANSP(IC,JC,KC) = TRANSP(IC,JC,KC)/DRHS(IC,JC,KC)

C           ===================================================================

C           EVALUATE PRESSURE
            PRUN(IC,JC,KC) = TRUN(IC,JC,KC)*STORE7(IC,JC,KC)

C           ===================================================================

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================


      RETURN
      END
