      SUBROUTINE DIFFIN 
 
C     *************************************************************************
C
C     DIFFIN
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     06-JAN-2013:  CREATED
C     01-APR-2022:  RSC REFITTED TRANSPORT COEFFICIENTS
C     26-MAY-2023:  RSC BUG FIX THERMAL DIFFUSION RATIO
C     22-MAR-2026:  RSC BUG FIX THERMAL DIFFUSION RATIO
C      
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES DATA FOR MULTI-SPECIES MIXTURE AVERAGED MOLECULAR TRANSPORT
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     PARAMETER
C     =========
      DOUBLE PRECISION WMLTOL
      PARAMETER(WMLTOL = 1.0D-3)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION PCOUNT,FORNOW
      INTEGER ISPEC,JSPEC,KSPEC,ICOEFF,NDSPEC
      INTEGER NCOUNT
      INTEGER IROOT
      CHARACTER*10 SPCDIF(NSPCMX)
 

C     BEGIN
C     =====

C     =========================================================================

C     CONTROL FLAGS
C     DEFAULTS: MIXTURE AVERAGED TRANSPORT NOT ACTIVATED
      FLMAVT = .FALSE.
      FLMIXW = .FALSE.
      FLMIXP = .FALSE.
      FLMIXT = .FALSE.
      DO ISPEC = 1, NSPEC
        FLMSOR(ISPEC) = .FALSE.
        FLMDUF(ISPEC) = .FALSE.
        FLMTDR(ISPEC) = .FALSE.
      ENDDO

C     =========================================================================

C     CHECK SWITCH FOR MIXTURE AVERAGED TRANSPORT
      IF(NFMAVT.EQ.1)THEN

C       =======================================================================

C       DATA IS READ BY THE LOWEST-RANKED PROCESSOR AND BROADCAST TO THE OTHERS
C       -----------------------------------------------------------------------

C       SET ID OF ROOT (BROADCASTING) PROCESSOR
        IROOT = 0

C       INITIALISE THE BROADCAST COUNTER
        PCOUNT = ZERO

C       INITIALISE THE BROADCAST ARRAY
        DO NCOUNT = 1, NPARAY
          PARRAY(NCOUNT) = ZERO
        ENDDO

C       =======================================================================

C       LOWEST-RANKED PROCESSOR
        IF(IPROC.EQ.0)THEN

C       =======================================================================

C         READ THE MOLECULAR TRANSPORT DATA FILE
C         --------------------------------------
C         ---------------------------------------------------------------------

          OPEN(UNIT=NCDIFF,FILE=FNDIFF,STATUS='OLD',FORM='FORMATTED')

C         ---------------------------------------------------------------------

C         READ AND IGNORE THE HEADER
          READ(NCDIFF,*)
          READ(NCDIFF,*)
          READ(NCDIFF,*)
          READ(NCDIFF,*)
          READ(NCDIFF,*)

C         ---------------------------------------------------------------------

C         READ SPECIES LIST AND CHECK AGAINST CHEMICAL DATA
          READ(NCDIFF,*)
          READ(NCDIFF,'(I5)')NDSPEC
          IF(NDSPEC.NE.NSPEC)THEN
            WRITE(NCREPT,*)
            WRITE(NCREPT,*)'Warning: DIFFIN: mismatch in no. of species'
            WRITE(NCREPT,*)'CHEMIN: ',NSPEC,'; DIFFIN: ',NDSPEC
            WRITE(NCREPT,*)
          ENDIF

          DO JSPEC = 1, NSPEC
            READ(NCDIFF,'(I5,3X,A)')ISPEC,SPCDIF(ISPEC)
            IF(SPCSYM(ISPEC).NE.SPCDIF(ISPEC))THEN
              WRITE(NCREPT,*)
            WRITE(NCREPT,*)'Warning: DIFFIN: mismatch in species symbol'
              WRITE(NCREPT,*)'CHEMIN: ',SPCSYM(ISPEC),
     +                     '; DIFFIN: ',SPCDIF(ISPEC)
              WRITE(NCREPT,*)
            ENDIF
          ENDDO

C         MOLECULAR TRANSPORT DATA
          READ(NCDIFF,*)
          READ(NCDIFF,*)PDIFGB,TDIFGB

          READ(NCDIFF,*)
          DO ISPEC = 1,NSPEC
            READ(NCDIFF,*)JSPEC,NCOVIS
            IF(JSPEC.NE.ISPEC)THEN
              WRITE(NCREPT,*)
              WRITE(NCREPT,*)'Warning: DIFFIN: species ID mismatch'
              WRITE(NCREPT,*)'Viscosity data'
              WRITE(NCREPT,*)'Expected: ',ISPEC,'; Input: ',JSPEC
            ENDIF
            READ(NCDIFF,*)
     +          (VISCCO(ICOEFF,ISPEC),ICOEFF=1,NCOVIS)
          ENDDO

          READ(NCDIFF,*)
          DO ISPEC = 1,NSPEC
            READ(NCDIFF,*)JSPEC,NCOCON
            IF(JSPEC.NE.ISPEC)THEN
              WRITE(NCREPT,*)
              WRITE(NCREPT,*)'Warning: DIFFIN: species ID mismatch'
              WRITE(NCREPT,*)'Thermal conductivity data'
              WRITE(NCREPT,*)'Expected: ',ISPEC,'; Input: ',JSPEC
            ENDIF
            READ(NCDIFF,*)
     +          (CONDCO(ICOEFF,ISPEC),ICOEFF=1,NCOCON)
          ENDDO

          READ(NCDIFF,*)
          DO ISPEC = 1,NSPEC
            READ(NCDIFF,*)JSPEC,NCODIF
            IF(JSPEC.NE.ISPEC)THEN
              WRITE(NCREPT,*)
              WRITE(NCREPT,*)'Warning: DIFFIN: species ID mismatch'
              WRITE(NCREPT,*)'Binary diffusion coefficient data ISPEC'
              WRITE(NCREPT,*)'Expected: ',ISPEC,'; Input: ',JSPEC
            ENDIF
            DO JSPEC = 1,ISPEC
              READ(NCDIFF,*)KSPEC,
     +            (DIFFCO(ICOEFF,JSPEC,ISPEC),ICOEFF=1,NCODIF)
              IF(KSPEC.NE.JSPEC)THEN
                WRITE(NCREPT,*)
                WRITE(NCREPT,*)'Warning: DIFFIN: species ID mismatch'
                WRITE(NCREPT,*)'Binary diffusion coefficient data JSPEC'
                WRITE(NCREPT,*)'Expected: ',JSPEC,'; Input: ',KSPEC
              ENDIF
            ENDDO
          ENDDO

          READ(NCDIFF,*)
          DO ISPEC = 1,NSPEC
            READ(NCDIFF,*)JSPEC,NCOTDR
            IF(JSPEC.NE.ISPEC)THEN
              WRITE(NCREPT,*)
              WRITE(NCREPT,*)'Warning: DIFFIN: species ID mismatch'
              WRITE(NCREPT,*)'Thermal diffusion coefficient data ISPEC'
              WRITE(NCREPT,*)'Expected: ',ISPEC,'; Input: ',JSPEC
            ENDIF
            DO JSPEC = 1,ISPEC-1
              READ(NCDIFF,*)KSPEC,
     +            (TDRCCO(ICOEFF,JSPEC,ISPEC),ICOEFF=1,NCOTDR)
              IF(KSPEC.NE.JSPEC)THEN
                WRITE(NCREPT,*)
                WRITE(NCREPT,*)'Warning: DIFFIN: species ID mismatch'
              WRITE(NCREPT,*)'Thermal diffusion coefficient data JSPEC'
                WRITE(NCREPT,*)'Expected: ',ISPEC,'; Input: ',JSPEC
              ENDIF
            ENDDO
          ENDDO

C         ---------------------------------------------------------------------

C         READ END-OF-FILE LINE
          READ(NCDIFF,*)
          CLOSE(NCDIFF)

C         =====================================================================

C         COMPRESS THE DATA
C         -----------------

C         MOLECULAR TRANSPORT DATA
          NCOUNT = 1
          PARRAY(NCOUNT) = PDIFGB
          NCOUNT = NCOUNT + 1
          PARRAY(NCOUNT) = TDIFGB
          DO ISPEC = 1,NSPEC
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NCOVIS)
            DO ICOEFF = 1, NCOVIS
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = VISCCO(ICOEFF,ISPEC)
            ENDDO
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NCOCON)
            DO ICOEFF = 1, NCOCON
              NCOUNT = NCOUNT + 1
              PARRAY(NCOUNT) = CONDCO(ICOEFF,ISPEC)
            ENDDO
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NCODIF)
            DO JSPEC = 1, ISPEC
              DO ICOEFF = 1, NCODIF
                NCOUNT = NCOUNT + 1
                PARRAY(NCOUNT) = DIFFCO(ICOEFF,JSPEC,ISPEC)
              ENDDO
            ENDDO
            NCOUNT = NCOUNT + 1
            PARRAY(NCOUNT) = REAL(NCOTDR)
            DO JSPEC = 1, ISPEC-1
              DO ICOEFF = 1, NCOTDR
                NCOUNT = NCOUNT + 1
                PARRAY(NCOUNT) = TDRCCO(ICOEFF,JSPEC,ISPEC)
              ENDDO
            ENDDO
          ENDDO

C         ---------------------------------------------------------------------

C         BROADCAST COUNTER
          PCOUNT = REAL(NCOUNT)

C         CHECK BROADCAST COUNTER AGAINST BROADCAST ARRAY SIZE
          IF(NCOUNT.GT.NPARAY)THEN
          WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Warning: DIFFIN: broadcast array size too small'
          WRITE(NCREPT,*)'Actual: ',NCOUNT,'; Available: ',NPARAY
          WRITE(NCREPT,*)
          ENDIF

C         =====================================================================

        ENDIF

C       =======================================================================

C       BROADCAST (OR RECEIVE) THE COUNTER
C       ----------------------------------
        CALL P_BCST(PCOUNT,1,1,IROOT)
        NCOUNT = NINT(PCOUNT)

C       BROADCAST (OR RECEIVE) THE DATA
C       -------------------------------
        CALL P_BCST(PARRAY,NPARAY,NCOUNT,IROOT)

C       =======================================================================

C       NOT THE LOWEST-RANKED PROCESSOR
        IF(IPROC.NE.0)THEN

C         =====================================================================

C         UNCOMPRESS THE DATA
C         -------------------

C         MOLECULAR TRANSPORT DATA
          NCOUNT = 1
          PDIFGB = PARRAY(NCOUNT)
          NCOUNT = NCOUNT + 1
          TDIFGB = PARRAY(NCOUNT)
          DO ISPEC = 1,NSPEC
            NCOUNT = NCOUNT + 1
            NCOVIS = NINT(PARRAY(NCOUNT))
            DO ICOEFF = 1, NCOVIS
              NCOUNT = NCOUNT + 1
              VISCCO(ICOEFF,ISPEC) = PARRAY(NCOUNT)
            ENDDO
            NCOUNT = NCOUNT + 1
            NCOCON = NINT(PARRAY(NCOUNT))
            DO ICOEFF = 1, NCOCON
              NCOUNT = NCOUNT + 1
              CONDCO(ICOEFF,ISPEC) = PARRAY(NCOUNT)
            ENDDO
            NCOUNT = NCOUNT + 1
            NCODIF = NINT(PARRAY(NCOUNT))
            DO JSPEC = 1, ISPEC
              DO ICOEFF = 1, NCODIF
                NCOUNT = NCOUNT + 1
                DIFFCO(ICOEFF,JSPEC,ISPEC) = PARRAY(NCOUNT)
              ENDDO
            ENDDO
            NCOUNT = NCOUNT + 1
            NCOTDR = NINT(PARRAY(NCOUNT))
            DO JSPEC = 1, ISPEC-1
              DO ICOEFF = 1, NCOTDR
                NCOUNT = NCOUNT + 1
                TDRCCO(ICOEFF,JSPEC,ISPEC) = PARRAY(NCOUNT)
              ENDDO
            ENDDO
          ENDDO

C         =====================================================================

        ENDIF

C       =======================================================================

C       EVALUATE DERIVED QUANTITIES
C       ===========================

C       COEFFICIENT COUNTERS
        NCOVM1 = NCOVIS - 1
        NCOCM1 = NCOCON - 1
        NCODM1 = NCODIF - 1
        NCOTM1 = NCOTDR - 1

C       REFERENCE TEMPERATURE AND PRESSURE
        TDIFGB = ONE/TDIFGB
        PDIFGB = ONE/PDIFGB

C       FILL OUT THE DIFFUSION COEFFICIENT MATRIX
C       MATRIX IS SYMMETRIC
        DO ISPEC = 1, NSPEC
          DO JSPEC = ISPEC+1, NSPEC
            DO ICOEFF = 1, NCODIF
              FORNOW = DIFFCO(ICOEFF,ISPEC,JSPEC)
              DIFFCO(ICOEFF,JSPEC,ISPEC) = FORNOW
            ENDDO
          ENDDO
        ENDDO

C       FILL OUT THE THERMAL DIFFUSION COEFFICIENT MATRIX
C       RSC 26-MAY-2023
        DO ISPEC = 1, NSPEC
C         LEADING DIAGONAL IS ZERO
          DO ICOEFF = 1, NCOTDR
            TDRCCO(ICOEFF,ISPEC,ISPEC) = ZERO
          ENDDO
C         RSC 22-MAR-2026
C         MATRIX IS ANTISYMMETRIC
C         ELEMENTS OF THE MATRIX TRANSPOSE ARE STORED
C          DO JSPEC = 1, ISPEC-1
C            DO ICOEFF = 1, NCOTDR
C              FORNOW = TDRCCO(ICOEFF,JSPEC,ISPEC)
C              TDRCCO(ICOEFF,JSPEC,ISPEC) = -FORNOW
C            ENDDO
C          ENDDO
          DO JSPEC = ISPEC+1, NSPEC
            DO ICOEFF = 1, NCOTDR
              FORNOW = TDRCCO(ICOEFF,ISPEC,JSPEC)
C              TDRCCO(ICOEFF,JSPEC,ISPEC) = FORNOW
              TDRCCO(ICOEFF,JSPEC,ISPEC) = -FORNOW
            ENDDO
          ENDDO
        ENDDO

C       PRECOMPUTE ELEMENTS OF WILKES COMBINATION RULE FOR VISCOSITY
C       NOTE THAT THE WILKES MATRIX IS NOT SYMMETRIC
C       ELEMENTS OF THE MATRIX TRANSPOSE ARE STORED
        DO ISPEC = 1, NSPEC
          DO JSPEC = 1, NSPEC
            FORNOW = WMOLAR(ISPEC)/WMOLAR(JSPEC)
            WILKO2(JSPEC,ISPEC) = ONE/SQRT(SQRT(FORNOW))
            FORNOW = EIGHT*(ONE + FORNOW)
            WILKO1(JSPEC,ISPEC) = ONE/SQRT(FORNOW)
          ENDDO
        ENDDO

C       =======================================================================

C       CONTROL FLAGS
        FLMAVT = .TRUE.
        IF(NFMIXW.EQ.1)FLMIXW = .TRUE.
        IF(NFMIXP.EQ.1)FLMIXP = .TRUE.

        IF(NFMSOR.EQ.1)THEN
          FLMIXT = .TRUE.
          DO ISPEC = 1, NSPEC
            IF(WMOLAR(ISPEC).LE.(WMLTDR+WMLTOL))THEN
              FLMSOR(ISPEC) = .TRUE.
              FLMTDR(ISPEC) = .TRUE.
            ENDIF
          ENDDO
        ENDIF

        IF(NFMDUF.EQ.1)THEN
          FLMIXT = .TRUE.
          DO ISPEC = 1, NSPEC
            IF(WMOLAR(ISPEC).LE.(WMLTDR+WMLTOL))THEN
              FLMDUF(ISPEC) = .TRUE.
              FLMTDR(ISPEC) = .TRUE.
            ENDIF
          ENDDO
        ENDIF

C       =======================================================================

      ENDIF

C     =========================================================================


      RETURN
      END
