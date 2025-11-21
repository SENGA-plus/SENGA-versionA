      SUBROUTINE TURBIN 
 
C     *************************************************************************
C
C     TURBIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     28-MAR-1997:  CREATED
C     20-OCT-1999:  KWJ PARALLEL FFT IMPLEMENTATION
C     24-MAY-2003:  RSC UPDATED FOR SENGA2
C     08-SEP-2012:  RSC/RACG BUG FIX SPECIAL CASE OF ZERO 12-PLANE VECTOR
C 
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIAL TURBULENCE GENERATOR
C     COMPUTES A HOMOGENEOUS ISOTROPIC TURBULENT VELOCITY FIELD
C     SATISFYING CONTINUITY, ISOTROPY AND TOTAL ENERGY CONDITIONS
C
C     FIELD IS BUILT IN FOURIER SPACE ACCORDING TO ROGALLO (NASA TM 81315)
C     
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     FUNCTIONS
C     =========
      DOUBLE PRECISION ESPECT,RANUNI
      EXTERNAL ESPECT,RANUNI


C     PARAMETERS
C     ==========
      DOUBLE PRECISION VECTOL,TOLIMG
      PARAMETER(VECTOL=1.0D-5, TOLIMG=1.0E-6)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION VECTK1,VECTK2,VECTK3,VKSIZE,OVKSIZ
      DOUBLE PRECISION VELMAG,VFACTR,PLNMAG
      DOUBLE PRECISION OVPLMG,AZIANG,COSAZI,SINAZI
      DOUBLE PRECISION PHANG1,PHANG2,COSPH1,SINPH1,COSPH2,SINPH2
      DOUBLE PRECISION VATERM,VBTERM
      DOUBLE PRECISION TKLODD
      DOUBLE PRECISION TWOPI,OVTOPI
      DOUBLE PRECISION UBART,VBART,WBART,UVART,VVART,WVART,TKET
      DOUBLE PRECISION UBARTT,VBARTT,WBARTT,UVARTT,VVARTT,WVARTT,TKETOT
      DOUBLE PRECISION UBARTL,VBARTL,WBARTL,UVARTL,VVARTL,WVARTL,TKETL
      DOUBLE PRECISION UBARTG,VBARTG,WBARTG,UVARTG,VVARTG,WVARTG,TKETG
      DOUBLE PRECISION UDEV,VDEV,WDEV,FACLAV,FACGAV
      INTEGER IC,JC,KC,IX,JX,KX,ICPROC
      INTEGER IGOFST,JGOFST,KGOFST,IGOFM1,JGOFM1,KGOFM1
      INTEGER IGSTAL,JGSTAL,KGSTAL,IGSTOL,JGSTOL,KGSTOL
      INTEGER NODBLX,NODBLY,NODBLZ
      INTEGER NGODDX,NGODDY,NGODDZ
      INTEGER ISEED
      LOGICAL FLAGIM
      LOGICAL FLRANI,FLRANJ,FLRANK


C     BEGIN
C     =====

C     =========================================================================

C     INDEXING
C     --------

C     SET ODD-NUMBER INDICATORS
      NGODDX = MOD(NXGLBL,2)
      NGODDY = MOD(NYGLBL,2)
      NGODDZ = MOD(NZGLBL,2)

C     SET ODDBALL WAVENUMBER INDICES
      NODBLX = NXGLBL/2 - 1 + NGODDX
      NODBLY = NYGLBL/2 - 1 + NGODDY
      NODBLZ = NZGLBL/2 - 1 + NGODDZ

C     PHYSICAL-SPACE GLOBAL INDEX OFFSETS
      IGOFST = 0
      DO ICPROC = 0, IXPROC-1
        IGOFST = IGOFST + NPMAPX(ICPROC)
      ENDDO
      IGSTAL = IGOFST+1
      IGSTOL = IGOFST + NPMAPX(IXPROC)
      IGOFM1 = IGOFST-1

      JGOFST = 0
      DO ICPROC = 0, IYPROC-1
        JGOFST = JGOFST + NPMAPY(ICPROC)
      ENDDO
      JGSTAL = JGOFST+1
      JGSTOL = JGOFST + NPMAPY(IYPROC)
      JGOFM1 = JGOFST-1

      KGOFST = 0
      DO ICPROC = 0, IZPROC-1
        KGOFST = KGOFST + NPMAPZ(ICPROC)
      ENDDO
      KGSTAL = KGOFST+1
      KGSTOL = KGOFST + NPMAPZ(IZPROC)
      KGOFM1 = KGOFST-1

C     =========================================================================

C     SET OTHER FACTORS
      TWOPI = TWO*PI
      OVTOPI = ONE/TWOPI
      TKLODD = ZERO

C     =========================================================================

C     INITIALISE THE RANDOM NUMBER GENERATOR
C     --------------------------------------
      IF(INSEED.GE.0)THEN

C       LOCAL INITIALISATION
C       --------------------
C       RANDOM SEED MUST BE DIFFERENT FOR EVERY PROCESSOR
        ISEED = INSEED + IPROC
        CALL RANINI(ISEED)

C       SWEEP THROUGH LOCAL PHYSICAL SPACE
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

C             GET AND SAVE THREE RANDOM NUMBERS
              UTMP(IC,JC,KC) = RANUNI(ISEED)
              VTMP(IC,JC,KC) = RANUNI(ISEED)
              WTMP(IC,JC,KC) = RANUNI(ISEED)
              
            ENDDO
          ENDDO
        ENDDO

      ELSE

C       GLOBAL INITIALISATION
C       ---------------------
C       RANDOM SEED MUST BE IDENTICAL FOR EVERY PROCESSOR
        ISEED = INSEED
        CALL RANINI(ISEED)

C       SWEEP THROUGH GLOBAL PHYSICAL SPACE
        DO KC = 1,NZGLBL
          DO JC = 1,NYGLBL
            DO IC = 1,NXGLBL

C             GET THREE RANDOM NUMBERS
              AZIANG = RANUNI(ISEED)
              PHANG1 = RANUNI(ISEED)
              PHANG2 = RANUNI(ISEED)

C             CHECK PHYSICAL-SPACE INDEX RANGE FOR THIS PROCESSOR
              FLRANI = (IC.GE.IGSTAL).AND.(IC.LE.IGSTOL)
              FLRANJ = (JC.GE.JGSTAL).AND.(JC.LE.JGSTOL)
              FLRANK = (KC.GE.KGSTAL).AND.(KC.LE.KGSTOL)

              IF(FLRANI.AND.FLRANJ.AND.FLRANK)THEN

C               SET LOCAL INDEXING
                IX = IC - IGOFST
                JX = JC - JGOFST
                KX = KC - KGOFST

C               SAVE THE THREE RANDOM NUMBERS
                UTMP(IX,JX,KX) = AZIANG
                VTMP(IX,JX,KX) = PHANG1
                WTMP(IX,JX,KX) = PHANG2
              
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     INITIALISE THE ENERGY SPECTRUM
C     ------------------------------
      CALL ESPINI

C     =========================================================================

C     EVALUATE VELOCITY COMPONENTS IN FOURIER SPACE
C     ---------------------------------------------

C     SWEEP THROUGH LOCAL PHYSICAL SPACE
C     ----------------------------------
      DO KC = KSTAL,KSTOL

C       FOURIER-SPACE GLOBAL INDEXING
        KX = KGOFM1 + KC
        IF(KX.GT.NODBLZ)KX = KX-NZGLBL

        DO JC = JSTAL,JSTOL

C         FOURIER-SPACE GLOBAL INDEXING
          JX = JGOFM1 + JC
          IF(JX.GT.NODBLY)JX = JX-NYGLBL

          DO IC = ISTAL,ISTOL

C           FOURIER-SPACE GLOBAL INDEXING
            IX = IGOFM1 + IC
            IF(IX.GT.NODBLX)IX = IX-NXGLBL

C           CHECK LOCATION IN FOURIER SPACE
            FLRANK = (KX.GT.0)
            FLRANJ = (KX.EQ.0).AND.(JX.GT.0)
            FLRANI = (KX.EQ.0).AND.(JX.EQ.0).AND.(IX.GE.0)

            IF(FLRANI.OR.FLRANJ.OR.FLRANK)THEN
 
C             IN UPPER-CENTRAL HALF OF FOURIER SPACE

C             EVALUATE THE WAVENUMBER VECTOR COMPONENTS
              VECTK1 = REAL(IX)
              VECTK2 = REAL(JX)
              VECTK3 = REAL(KX)
              VKSIZE = SQRT(VECTK1*VECTK1+VECTK2*VECTK2+VECTK3*VECTK3)
C             SPECIAL CASE OF ZERO VECTOR
              OVKSIZ = ONE
              IF(VKSIZE.GT.VECTOL)OVKSIZ = ONE/VKSIZE

C             WAVENUMBER VECTOR MAGNITUDE IN THE 12-PLANE
              PLNMAG = SQRT(VECTK1*VECTK1 + VECTK2*VECTK2)
C             SPECIAL CASE OF ZERO 12-PLANE VECTOR
C             RSC/RACG 08-SEP-2012 BUG FIX
C              OVPLMG = ONE
C              IF(PLNMAG.GT.VECTOL)OVPLMG = ONE/PLNMAG
              IF(PLNMAG.LT.VECTOL)THEN
C               AZIMUTHAL ANGLE IS ARBITRARY/STOCHASTIC: CHOOSE K1=0 K2=1
                OVPLMG = ONE
                VECTK1 = ZERO
                VECTK2 = ONE
              ELSE
                OVPLMG = ONE/PLNMAG
              ENDIF
                  
C             EVALUATE THE ENERGY SPECTRUM FUNCTION
C             ENERGY SPECTRUM SPECIFIED EXTERNALLY IN FUNCTION ESPECT
              VELMAG = ESPECT(VKSIZE)

C             EVALUATE THE FOURIER-VELOCITY MAGNITUDE
              VELMAG = SQRT(VELMAG*OVTOPI)*OVKSIZ

C             GENERATE THE (RANDOM) AZIMUTH ANGLE
C             RANGE IS ZERO TO 2*PI
              AZIANG = UTMP(IC,JC,KC)
              AZIANG = TWOPI*AZIANG
              COSAZI = COS(AZIANG)
              SINAZI = SIN(AZIANG)

C             GENERATE THE (RANDOM) PHASE ANGLES
C             RANGE IS -PI TO PI
              PHANG1 = VTMP(IC,JC,KC)
              PHANG1 = PI*(TWO*PHANG1-ONE)
              COSPH1 = COS(PHANG1)
              SINPH1 = SIN(PHANG1)
              PHANG2 = WTMP(IC,JC,KC)
              PHANG2 = PI*(TWO*PHANG2-ONE)
              COSPH2 = COS(PHANG2)
              SINPH2 = SIN(PHANG2)

C             EVALUATE COMMON FACTOR FOR ALL COMPONENTS
              VFACTR = VELMAG*OVKSIZ*OVPLMG

C             U-COMPONENT
              VATERM = VFACTR*VKSIZE*VECTK2*COSAZI 
              VBTERM = VFACTR*VECTK1*VECTK3*SINAZI 
C             REAL PART
              URUN(IC,JC,KC) = VATERM*COSPH1 + VBTERM*COSPH2
C             IMAGINARY PART
              UTMP(IC,JC,KC) = VATERM*SINPH1 + VBTERM*SINPH2

C             V-COMPONENT
              VATERM = VFACTR*VKSIZE*VECTK1*COSAZI 
              VBTERM = VFACTR*VECTK2*VECTK3*SINAZI 
C             REAL PART
              VRUN(IC,JC,KC) = VBTERM*COSPH2 - VATERM*COSPH1
C             IMAGINARY PART
              VTMP(IC,JC,KC) = VBTERM*SINPH2 - VATERM*SINPH1

C             W-COMPONENT
C             VATERM IS ZERO
              VBTERM = -VFACTR*PLNMAG*PLNMAG*SINAZI 
C             REAL PART
              WRUN(IC,JC,KC) = VBTERM*COSPH2
C             IMAGINARY PART
              WTMP(IC,JC,KC) = VBTERM*SINPH2

            ELSE

C             IN LOWER-CENTRAL HALF OF FOURIER SPACE

C             SET VELOCITY COMPONENT VALUES TO ZERO FOR NOW
              URUN(IC,JC,KC) = ZERO
              UTMP(IC,JC,KC) = ZERO
              VRUN(IC,JC,KC) = ZERO
              VTMP(IC,JC,KC) = ZERO
              WRUN(IC,JC,KC) = ZERO
              WTMP(IC,JC,KC) = ZERO

            ENDIF

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================

C     CHECK ENERGY CONTENT
C     --------------------
      TKET = ZERO
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

             TKET = TKET
     +            + URUN(IC,JC,KC)*URUN(IC,JC,KC)
     +            + UTMP(IC,JC,KC)*UTMP(IC,JC,KC)
     +            + VRUN(IC,JC,KC)*VRUN(IC,JC,KC)
     +            + VTMP(IC,JC,KC)*VTMP(IC,JC,KC)
     +            + WRUN(IC,JC,KC)*WRUN(IC,JC,KC)
     +            + WTMP(IC,JC,KC)*WTMP(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     SUM OVER ALL PROCESSORS
      CALL P_SUMM(TKET,TKETOT)

C     REPORT
      IF(IPROC.EQ.0)THEN
        WRITE(NCREPT,*)'Fourier-space turbulence kinetic energy'
        WRITE(NCREPT,'(1PE12.4)')TKETOT
      ENDIF

C     =========================================================================

C     CARRY OUT AN INVERSE FOURIER TRANSFORM
C     ======================================
C     WITH CONJUGATE ANTISYMMETRY
      CALL TURBFT

C     =========================================================================

C     CHECK THAT TRANSFORMED DATA IS REAL
C     -----------------------------------

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            FLAGIM = .FALSE.
            IF(UTMP(IC,JC,KC).GT.TOLIMG)FLAGIM = .TRUE.
            IF(VTMP(IC,JC,KC).GT.TOLIMG)FLAGIM = .TRUE.
            IF(WTMP(IC,JC,KC).GT.TOLIMG)FLAGIM = .TRUE.
            IF(FLAGIM)THEN
              WRITE(6,*)'Warning: Imaginary part too large',IPROC
              WRITE(6,'(3I5,3(1PE12.4))')IC,JC,KC,
     +          UTMP(IC,JC,KC),VTMP(IC,JC,KC),WTMP(IC,JC,KC)
            ENDIF

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================

C     CHECK TURBULENCE STATISTICS
C     ---------------------------

C     AVERAGING FACTORS
      FACLAV = ONE/REAL(NXNODE)/REAL(NYNODE)/REAL(NZNODE)
      FACGAV = ONE/REAL(NXGLBL)/REAL(NYGLBL)/REAL(NZGLBL)


C     VELOCITY MEANS
C     --------------
C     SUM OVER LOCAL VALUES
      UBART = ZERO
      VBART = ZERO
      WBART = ZERO
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            UBART = UBART + URUN(IC,JC,KC)
            VBART = VBART + VRUN(IC,JC,KC)
            WBART = WBART + WRUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     LOCAL AVERAGES
      UBARTL = UBART*FACLAV
      VBARTL = VBART*FACLAV
      WBARTL = WBART*FACLAV

C     SUM OVER ALL PROCESSORS
      CALL P_SUMM(UBART,UBARTT)
      CALL P_SUMM(VBART,VBARTT)
      CALL P_SUMM(WBART,WBARTT)

C     GLOBAL AVERAGES
      UBARTG = UBARTT*FACGAV
      VBARTG = VBARTT*FACGAV
      WBARTG = WBARTT*FACGAV


C     VELOCITY VARIANCES
C     ------------------
C     SUM OVER LOCAL VALUES
      UVART = ZERO
      VVART = ZERO
      WVART = ZERO
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            UDEV = URUN(IC,JC,KC)-UBARTG
            UVART = UVART + UDEV*UDEV
            VDEV = VRUN(IC,JC,KC)-VBARTG
            VVART = VVART + VDEV*VDEV
            WDEV = WRUN(IC,JC,KC)-WBARTG
            WVART = WVART + WDEV*WDEV

          ENDDO
        ENDDO
      ENDDO

C     TURBULENCE TOTAL ENERGY
      TKET   = HALF*(UVART+VVART+WVART)

C     LOCAL AVERAGES
      UVARTL = UVART*FACLAV
      VVARTL = VVART*FACLAV
      WVARTL = WVART*FACLAV
      TKETL  = HALF*(UVARTL+VVARTL+WVARTL)

C     SUM OVER ALL PROCESSORS
      CALL P_SUMM(UVART,UVARTT)
      CALL P_SUMM(VVART,VVARTT)
      CALL P_SUMM(WVART,WVARTT)
      CALL P_SUMM(TKET,TKETOT)

C     GLOBAL AVERAGES
      UVARTG = UVARTT*FACGAV
      VVARTG = VVARTT*FACGAV
      WVARTG = WVARTT*FACGAV
      TKETG  = TKETOT*FACGAV

C     REPORT
      IF(IPROC.EQ.0)THEN
        WRITE(NCREPT,*)'Physical-space turbulence kinetic energy'
        WRITE(NCREPT,'(1PE12.4)')TKETG
        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Velocity means:'
        WRITE(NCREPT,'(3(1PE12.4))')UBARTG,VBARTG,WBARTG
        WRITE(NCREPT,*)'Velocity variances:'
        WRITE(NCREPT,'(3(1PE12.4))')UVARTG,VVARTG,WVARTG
        WRITE(NCREPT,*)
      ENDIF
       
C     =========================================================================


      RETURN
      END
