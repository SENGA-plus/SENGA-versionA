      SUBROUTINE INDATA
 
C     *************************************************************************
C
C     INDATA
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT
C
C     CHANGE RECORD
C     -------------
C     01-AUG-1996:  CREATED
C     12-MAR-2003:  RSC UPDATED FOR SENGA2
C     11-JUL-2009:  RSC ADD A DUMP FORMAT SWITCH
C     29-AUG-2009:  RSC UPDATE NUMBER OF PROCESSORS
C     06-JAN-2013:  RSC MIXTURE-AVERAGED TRANSPORT
C     14-JUL-2013:  RSC RADIATION TREATMENT
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     INITIALISES ALL DATA
C
C     FILES: CHEMICAL DATA FILE  - UNIT NCCHEM - FORMATTED INPUT
C            TRANSPORT DATA FILE - UNIT NCDIFF - FORMATTED INPUT
C            RADIATION DATA FILE - UNIT NCRADN - FORMATTED INPUT
C            REPORT FILE         - UNIT NCREPT - FORMATTED OUTPUT
C            STATISTICS FILE     - UNIT NCSTAT - FORMATTED OUTPUT
C            DUMP OUTPUT FILES   - UNIT NCDMPO - UNFORMATTED      
C            DUMP INPUT FILES    - UNIT NCDMPI - UNFORMATTED  (RESTART ONLY)     
C        
C     *************************************************************************

!START LR: Removal
!#ifdef HDF5
!      USE HDF5IO
!#endif
!END

C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
      INCLUDE 'com_espect.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
C     FILENAME COMPONENTS
!START: LR: Replacement LR:Reverting
      CHARACTER*10 PNCONT,PNCHEM,PNDIFF,PNRADN
      CHARACTER*11 PNREPT,PNSTAT,PNDMPI,PNDMPO
!      CHARACTER*4 PNCONT,PNCHEM,PNDIFF,PNRADN
!      CHARACTER*4 PNREPT,PNSTAT,PNDMPI,PNDMPO
!END

!START: LR: Replacement
!      CHARACTER*4 PNXDAT,PNXRES
      CHARACTER*4 PNXDAT,PNXRES,PNTIXL,PNTCXL
      CHARACTER*4 PNXHDF,PNXBIN,PNEXT
      CHARACTER*6 PNTIME
!END

!START LR: Replacement LR:Reverting
      PARAMETER(PNCONT = 'input/cont',
     +          PNCHEM = 'input/chem',
     +          PNDIFF = 'input/diff', PNRADN = 'input/radn',
     +          PNREPT = 'output/rept', PNSTAT = 'output/stat',
     +          PNDMPI = 'output/dmpi', PNDMPO = 'output/dmpo',
!      PARAMETER(PNCONT = 'cont',
!     +          PNCHEM = 'chem', PNDIFF = 'diff', PNRADN = 'radn',
!     +          PNREPT = 'rept', PNSTAT = 'stat',
!     +          PNDMPI = 'dmpi', PNDMPO = 'dmpo',
!END
!START LR: Removal
!     +          PNXDAT = '.dat', PNXRES = '.res')
!END
!START LR: Addition
     +          PNXDAT = '.dat', PNXRES = '.res',
     +          PNXHDF = '.h5',  PNXBIN = '.bin',
     +          PNTIXL = 'tixl', PNTCXL = 'tcxl')
!END

C     REPORT FILE LINE LENGTH
      INTEGER LENLIN
      PARAMETER(LENLIN = 79)

C     THERMO DATA INTERVAL-MATCHING TOLERANCE
      DOUBLE PRECISION TTTOL
      PARAMETER(TTTOL = 1.0D-3)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION DYRIN(NSPCMX)
      DOUBLE PRECISION CTRANS(NSPCMX,4)
      DOUBLE PRECISION DTRANS(NSPCMX)
      DOUBLE PRECISION RTIN,DURIN,DVRIN,DWRIN,DERIN
      DOUBLE PRECISION TTEMP(5),TTOLD(5)
!START: LR: Replacement LR:Reverting
      DOUBLE PRECISION FORNOW
!      DOUBLE PRECISION FORNOW,TMPLOG
!END

      DOUBLE PRECISION COMBO1,COMBO2,COMBO3
      INTEGER IINDEX,IPOWER,ICOEF1,ICOEF2
      INTEGER IC,JC,KC,ISPEC,ISTEP,ITINT,ICP
!START: LR: Replacement LR: Reverting
      INTEGER JSPEC,NCOUNT
!      INTEGER JSPEC,NCOUNT,ICOEFF
!END
      INTEGER NXDMAX,NYDMAX,NZDMAX,NDSPEC
      CHARACTER*132 STROUT
      CHARACTER*79 CLSTAR,CLDASH
      CHARACTER*10 STREAC
      CHARACTER*5 STRSTP
      CHARACTER*4 STRARO
      CHARACTER*1 STRCOF
C     RSC UPDATE NUMBER OF PROCESSORS
      CHARACTER*6 PNPROC
!START: LR: Removal LR: Reverting
      CHARACTER*1 PNFLAG
!END
      CHARACTER*20 FNAME
      LOGICAL FXDUMP
      LOGICAL FLGREQ
      

C     BEGIN
C     =====

C     =========================================================================
C     SET UP HDF5 I/O
C     ===============
!START LR: Removal
!#ifdef HDF5
!      CALL HDF5_INIT
!#endif
!END

C     SET UP FILE I/O
C     ===============

C     UNIT NUMBERS
      NCCONT = 1
      NCCHEM = 2
      NCDIFF = 3
      NCRADN = 4
      NCDMPO = 7
      NCREPT = 8
      NCSTAT = 9
!START LR: Addition LR:Reverting
!      NCTIXL = 10
!END

C     BUILD THE FILENAMES
      WRITE(PNPROC,'(I6.6)')IPROC
      FNCONT = PNCONT//PNXDAT
      FNCHEM = PNCHEM//PNXDAT
      FNDIFF = PNDIFF//PNXDAT
      FNRADN = PNRADN//PNXDAT
      FNREPT = PNREPT//PNXRES
      FNSTAT = PNSTAT//PNXRES
      IDFLAG = 0
!START LR: Removal
!#ifdef HDF5
!      H5_FILENAME(1) = PNDMPI//PNFLAG//".h5"
!#endif
!END
!START LR: Removal LR: Reverting     
      WRITE(PNFLAG,'(I1)')IDFLAG
      FNDMPO(1) = PNDMPI//PNPROC//PNFLAG//PNXDAT
      IDFLAG = 1
      WRITE(PNFLAG,'(I1)')IDFLAG
      FNDMPO(2) = PNDMPI//PNPROC//PNFLAG//PNXDAT
!END
!START LR: Removal
!#ifdef HDF5      
!      H5_FILENAME(2) = PNDMPI//PNFLAG//".h5"
!#endif
!END

C     =========================================================================

C     GET THE RUN CONTROL DATA
C     ========================
      CALL CONTIN

      IF(NXLPRM(1) .EQ. 4) THEN
C     VM: SYNTHETIC DIGITAL FILTERING METHOD
        DO JC=JSTAL,JSTOL
          DO KC=KSTAL,KSTOL
            UINF1(JC,KC) = ZERO
            UINF2(JC,KC) = ZERO
            VINF1(JC,KC) = ZERO
            VINF2(JC,KC) = ZERO
            WINF1(JC,KC) = ZERO
            WINF2(JC,KC) = ZERO
            DO ISPEC=1,NSPEC
              YINF1(JC,KC,ISPEC)=ZEROS
              YINF2(JC,KC,ISPEC)=ZEROS
            ENDDO
          ENDDO
        ENDDO
C     ATTENTION WHEN RESTART: THE SAME FILE EXTENSION MUST BE ENTERED HERE AS IN INFLOW!
        IF (NCDMPI.EQ.1 .AND. IPROC.EQ.0)THEN
C     ACTUALLY, INTRAN0.DAT OR INTRAN1.DAT SHOULD BE HERE!
          OPEN(UNIT=16,FILE='output/intran.dat',FORM='FORMATTED') 
          READ(16,*)INTRAN
          CLOSE(16)    
          PASSER=DBLE(INTRAN)
          CALL P_BCST(PASSER,1,1,0)
        ELSE
          INTRAN=-1
        ENDIF

        CALL P_SYNC
        INTRAN=NINT(PASSER)

        IF(NCDMPI.EQ.1)THEN
          IF (IXPROC.EQ.0)THEN
            WRITE(FNAME,'(A,I4.4,A)')'output/inflow'//PNPROC//'.dat'
            OPEN(UNIT=16,FILE=FNAME,STATUS='unknown',FORM='FORMATTED')
            DO KC = 1,NZSIZE
              DO JC = 1,NYSIZE
                READ(16,*)UINF2(JC,KC),VINF2(JC,KC),WINF2(JC,KC),
     +                    (YINF2(JC,KC,ISPEC),ISPEC=1,NSPEC)
              ENDDO
            ENDDO
            CLOSE(16)

            DO JC=JSTAL,JSTOL
              DO KC=KSTAL,KSTOL
                UINF1(JC,KC)=UINF2(JC,KC)
                VINF1(JC,KC)=VINF2(JC,KC)
                WINF1(JC,KC)=WINF2(JC,KC)
                DO ISPEC=1,NSPEC
                  YINF1(JC,KC,ISPEC)=YINF2(JC,KC,ISPEC)
                ENDDO
              ENDDO
            ENDDO
          END IF
        END IF
      END IF

C     =========================================================================

C     WRITE THE INITIAL REPORT
C     ========================
     
      IF(IPROC.EQ.0)THEN

C       SET UP SEPARATOR LINES
        DO IC = 1,LENLIN
          CLSTAR(IC:IC) = '*'
          CLDASH(IC:IC) = '-'
        ENDDO

        OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='NEW',FORM='FORMATTED')

C       WRITE A HEADER
        WRITE(NCREPT,9000)CLSTAR
        WRITE(NCREPT,9010)
        WRITE(NCREPT,9020)
        WRITE(NCREPT,9010)
        WRITE(NCREPT,9030)
        WRITE(NCREPT,9010)
        WRITE(NCREPT,9040)
        WRITE(NCREPT,9010)
        WRITE(NCREPT,9000)CLSTAR
        WRITE(NCREPT,*)

        WRITE(NCREPT,*)'==========='
        WRITE(NCREPT,*)'REPORT FILE'
        WRITE(NCREPT,*)'==========='
        WRITE(NCREPT,*)

        WRITE(NCREPT,9000)CLDASH

C       DOMAIN DECOMPOSITION DATA AND CHECKING
        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Domain decomposition data:'
        WRITE(NCREPT,*)'-------------------------'
        WRITE(NCREPT,*)

C       GLOBAL GRID SIZES
        WRITE(NCREPT,*)'NXGLBL,NYGLBL,NZGLBL'
        WRITE(NCREPT,9200)NXGLBL,NYGLBL,NZGLBL
        WRITE(NCREPT,*)

C       CHECK GLOBAL GRID SIZE
        FLGREQ = (NXGREQ.EQ.NXGLBL)
     +      .AND.(NYGREQ.EQ.NYGLBL)
     +      .AND.(NZGREQ.EQ.NZGLBL)
        IF(.NOT.FLGREQ)THEN
          WRITE(NCREPT,*)'Warning: INDATA: global domain size mismatch'
          WRITE(NCREPT,*)'Required values in control file:'
          WRITE(NCREPT,9200)NXGREQ,NYGREQ,NZGREQ
          WRITE(NCREPT,*)'Global values take precedence'
          WRITE(NCREPT,*)
        ENDIF

C       NUMBERS OF PROCESSORS
        WRITE(NCREPT,*)'NXPROC,NYPROC,NZPROC,NPROC'
        WRITE(NCREPT,9200)NXPROC,NYPROC,NZPROC,NPROC
        WRITE(NCREPT,*)
        IF(NPROC.NE.NXPROC*NYPROC*NZPROC)THEN
          WRITE(NCREPT,*)
     +    'Warning: INDATA: mismatched total no. of processors'
          WRITE(NCREPT,*)
        ENDIF

C       CHECK NUMBER OF PROCESSORS
        FLGREQ = (NXPREQ.EQ.NXPROC)
     +      .AND.(NYPREQ.EQ.NYPROC)
     +      .AND.(NZPREQ.EQ.NZPROC)
        IF(.NOT.FLGREQ)THEN
          WRITE(NCREPT,*)'Warning: INDATA: no. of processors mismatch'
          WRITE(NCREPT,*)'Required values in control file:'
          WRITE(NCREPT,9200)NXPREQ,NYPREQ,NZPREQ
          WRITE(NCREPT,*)'Global values take precedence'
          WRITE(NCREPT,*)
        ENDIF

C       LOCAL ARRAY SIZES
        WRITE(NCREPT,*)'NXSIZE,NYSIZE,NZSIZE'
        WRITE(NCREPT,9200)NXSIZE,NYSIZE,NZSIZE
        WRITE(NCREPT,*)

C       LOCAL GRID SIZES
        WRITE(NCREPT,*)'NPMAPX'
        WRITE(NCREPT,9200)(NPMAPX(IC),IC=0,NXPRM1)
        WRITE(NCREPT,*)'NPMAPY'
        WRITE(NCREPT,9200)(NPMAPY(IC),IC=0,NYPRM1)
        WRITE(NCREPT,*)'NPMAPZ'
        WRITE(NCREPT,9200)(NPMAPZ(IC),IC=0,NZPRM1)

C       DOMAIN DECOMPOSITION CHECKING
C       X-DIRECTION
        IINDEX = 0
        DO IC = 0,NXPRM1 
          IINDEX = IINDEX + NPMAPX(IC)
          IF(NPMAPX(IC).GT.NXSIZE)THEN
            WRITE(NCREPT,*)
     +      'Warning: local array bounds exceeded in x direction'
            WRITE(NCREPT,9200)IC,NPMAPX(IC),NXSIZE
          ENDIF
        ENDDO
        IF(IINDEX.NE.NXGLBL)THEN
          WRITE(NCREPT,*)
     +    'Warning: global grid size mismatch in x direction'
          WRITE(NCREPT,9200)IC,IINDEX,NXGLBL
        ENDIF

C       Y-DIRECTION
        IINDEX = 0
        DO IC = 0,NYPRM1 
          IINDEX = IINDEX + NPMAPY(IC)
          IF(NPMAPY(IC).GT.NYSIZE)THEN
            WRITE(NCREPT,*)
     +      'Warning: local array bounds exceeded in y direction'
            WRITE(NCREPT,9200)IC,NPMAPY(IC),NYSIZE
          ENDIF
        ENDDO
        IF(IINDEX.NE.NYGLBL)THEN
          WRITE(NCREPT,*)
     +    'Warning: global grid size mismatch in y direction'
          WRITE(NCREPT,9200)IC,IINDEX,NYGLBL
        ENDIF

C       Z-DIRECTION
        IINDEX = 0
        DO IC = 0,NZPRM1 
          IINDEX = IINDEX + NPMAPZ(IC)
          IF(NPMAPZ(IC).GT.NZSIZE)THEN
            WRITE(NCREPT,*)
     +      'Warning: local array bounds exceeded in z direction'
            WRITE(NCREPT,9200)IC,NPMAPZ(IC),NZSIZE
          ENDIF
        ENDDO
        IF(IINDEX.NE.NZGLBL)THEN
          WRITE(NCREPT,*)
     +    'Warning: global grid size mismatch in z direction'
          WRITE(NCREPT,9200)IC,IINDEX,NZGLBL
        ENDIF

C       PARALLEL TRANSFER ARRAY SIZE CHECKING
        NCOUNT = NHALOX*NYSIZE*NZSIZE
        IF(NCOUNT.GT.NPARAY)THEN
          WRITE(NCREPT,*)
     + 'Warning: parallel transfer array bounds exceeded in x direction'
          WRITE(NCREPT,9200)NCOUNT,NPARAY
        ENDIF
        NCOUNT = NHALOY*NXNBIG*NZSIZE
        IF(NCOUNT.GT.NPARAY)THEN
          WRITE(NCREPT,*)
     + 'Warning: parallel transfer array bounds exceeded in y direction'
          WRITE(NCREPT,9200)NCOUNT,NPARAY
        ENDIF
        NCOUNT = NHALOZ*NXNBIG*NYNBIG
        IF(NCOUNT.GT.NPARAY)THEN
          WRITE(NCREPT,*)
     + 'Warning: parallel transfer array bounds exceeded in z direction'
          WRITE(NCREPT,9200)NCOUNT,NPARAY
        ENDIF

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'End of domain decomposition data'
        WRITE(NCREPT,*)'--------------------------------'
        WRITE(NCREPT,*)
        WRITE(NCREPT,9000)CLDASH

C       =======================================================================

C       INITIAL DATA
        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Initial data:'
        WRITE(NCREPT,*)'------------'
        WRITE(NCREPT,*)

C       GLOBAL DOMAIN SIZE (x,y,z)
        WRITE(NCREPT,*)'XGDLEN,YGDLEN,ZGDLEN'
        WRITE(NCREPT,9400)XGDLEN,YGDLEN,ZGDLEN
        WRITE(NCREPT,*)

C       NUMERICAL CONTROL
        WRITE(NCREPT,*)'TSTEP,NTIME1,NTIME,NSTPSW,NTDUMP,NTREPT,NTSTAT'
        WRITE(NCREPT,9300)TSTEP,NTIME1,NTIME,NSTPSW,NTDUMP,NTREPT,NTSTAT
        WRITE(NCREPT,*)

C       COLD START SWITCH
        WRITE(NCREPT,*)'NCDMPI'
        WRITE(NCREPT,9200)NCDMPI

C       INITIAL TURBULENCE GENERATOR 
        WRITE(NCREPT,*)'INTURB, INSEED, SPARAM'
        WRITE(NCREPT,9270)INTURB,INSEED,(SPARAM(IC),IC=1,NSPARM)

C       FLAME START OPTION
        WRITE(NCREPT,*)'INFLAM'
        WRITE(NCREPT,9200)INFLAM

C       DEFAULT INITIAL DATA
        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'PRIN,TRIN,URIN,VRIN,WRIN'
        WRITE(NCREPT,9400)PRIN,TRIN,URIN,VRIN,WRIN
C       CHECK NUMBER OF SPECIES
        IF(NSPREQ.NE.NSPCMX)THEN
          WRITE(NCREPT,*)'Warning: INDATA: no. of species mismatch'
          WRITE(NCREPT,*)'Global value:'
          WRITE(NCREPT,9200)NSPCMX
          WRITE(NCREPT,*)'Required value in control file:'
          WRITE(NCREPT,9200)NSPREQ
          WRITE(NCREPT,*)'Global value takes precedence'
          WRITE(NCREPT,*)
        ENDIF
C       INITIAL SPECIES MASS FRACTIONS
        WRITE(NCREPT,*)'ISPEC,YRIN'
!START LR: Replacement LR:Reverting
        DO ISPEC = 1,NSPEC
!        DO ISPEC = 1,NSPREQ
!END
          WRITE(NCREPT,9450)ISPEC,YRIN(ISPEC)
        ENDDO
C       CHECK INITIAL SPECIES MASS FRACTIONS
        FORNOW = ZERO
!START LR: Replacement LR: Reverting
        DO ISPEC = 1,NSPEC
!        DO ISPEC = 1,NSPREQ
!END
          FORNOW = FORNOW + YRIN(ISPEC)
        ENDDO
        IF(ABS(FORNOW-ONE).GT.YTOLER)THEN
          WRITE(NCREPT,*)
     +    'Warning: INDATA: initial mass fractions do not sum to unity'
          WRITE(NCREPT,*)
        ENDIF

C       GLOBAL BOUNDARY CONDITION TYPES
        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'NGBCXL,NGBCXR,NGBCYL,NGBCYR,NGBCZL,NGBCZR'
        WRITE(NCREPT,9250)NGBCXL,NGBCXR,NGBCYL,NGBCYR,NGBCZL,NGBCZR

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'End of initial data'
        WRITE(NCREPT,*)'-------------------'
        WRITE(NCREPT,*)
        WRITE(NCREPT,9000)CLDASH

      ENDIF

C     =========================================================================
 
C     GET THE CHEMICAL DATA
C     =====================
      CALL CHEMIN

C     =========================================================================

      IF(IPROC.EQ.0)THEN

C       REPORT THE CHEMICAL DATA
C       ========================

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Chemical data:'
        WRITE(NCREPT,*)'-------------'
        WRITE(NCREPT,*)

C       SPECIES LIST
        WRITE(NCREPT,*)'Number of species:'
        WRITE(NCREPT,'(I5)')NSPEC
C       CHECK NUMBER OF SPECIES
        IF(NSPEC.NE.NSPCMX)THEN
          WRITE(NCREPT,*)'Warning: INDATA: no. of species mismatch'
          WRITE(NCREPT,*)'Global value:'
          WRITE(NCREPT,9200)NSPCMX
          WRITE(NCREPT,*)'Required value in chemical data file:'
          WRITE(NCREPT,9200)NSPEC
          WRITE(NCREPT,*)'Chemical data file value takes precedence'
          WRITE(NCREPT,*)
        ENDIF
        WRITE(NCREPT,*)'List of species:'
        WRITE(NCREPT,*)'  No.  Symbol         Mol. Mass'
        DO ISPEC = 1,NSPEC
          WRITE(NCREPT,'(I5,3X,A10,3X,1PE12.4)')
     +          ISPEC,SPCSYM(ISPEC),WMOLAR(ISPEC)
        ENDDO
        WRITE(NCREPT,*)

C       THERMODYNAMIC DATA
        WRITE(NCREPT,*)'Species thermodynamic data:'
        WRITE(NCREPT,'("  Reference pressure: ",1PE12.4)')PREFGB
        WRITE(NCREPT,*)' Spec.  No of T intervals   ',
     +                'Interval   T low       T high       No of coeffs'
        DO ISPEC = 1,NSPEC
          IC = 1
          WRITE(NCREPT,'(I5,6X,I5,9X,I5,8X,2(1PE12.4),I8)')
     +          ISPEC,NTINT(ISPEC),
     +          IC,TINTLO(IC,ISPEC),TINTHI(IC,ISPEC),NCOFCP(IC,ISPEC)
          DO IC = 2,NTINT(ISPEC)
            WRITE(NCREPT,'(25X,I5,8X,2(1PE12.4),I8)')
     +          IC,TINTLO(IC,ISPEC),TINTHI(IC,ISPEC),NCOFCP(IC,ISPEC)

          ENDDO
        ENDDO
        WRITE(NCREPT,*)'Cp coeffs by mass'
        WRITE(NCREPT,*)' Spec.  T int.  Coeff no.  Coeff.'
        DO ISPEC = 1,NSPEC
          DO IC = 1,NTINT(ISPEC)
            JC = 1
            WRITE(NCREPT,'(I5,2X,I5,5X,I5,4X,1PE15.7)')
     +            ISPEC,IC,JC,AMASCP(JC,IC,ISPEC)
            DO JC = 2,NCOFCP(IC,ISPEC)
              WRITE(NCREPT,'(17X,I5,4X,1PE15.7)')JC,AMASCP(JC,IC,ISPEC)
            ENDDO
          ENDDO
        ENDDO

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Mass-specific Cp, Enthalpy, Entropy;',
     +                 '  Molar Gibbs fn.:'
        WRITE(NCREPT,*)' Spec.  T int.  Temp.       Cp',
     +                 '          Enthalpy    Entropy     Molar Gibbs'
        DO ISPEC = 1,NSPEC
          DO IC = 1,NTINT(ISPEC)
            DO ICP = 1,2

C             TEMPERATURE
              IF(ICP.EQ.1)TTEMP(1) = TINTLO(IC,ISPEC)
              IF(ICP.EQ.2)TTEMP(1) = TINTHI(IC,ISPEC)

C             MASS-SPECIFIC CP
              FORNOW = AMASCP(NCPOLY(IC,ISPEC),IC,ISPEC)
              DO JC = NCPOM1(IC,ISPEC),1,-1
                FORNOW = FORNOW*TTEMP(1) + AMASCP(JC,IC,ISPEC)
              ENDDO
              TTEMP(2) = FORNOW

C             MASS-SPECIFIC ENTHALPY
              FORNOW = AMASCH(NCPOLY(IC,ISPEC),IC,ISPEC)
              DO JC = NCPOM1(IC,ISPEC),1,-1
                FORNOW = FORNOW*TTEMP(1) + AMASCH(JC,IC,ISPEC)
              ENDDO
              FORNOW = AMASCH(NCENTH(IC,ISPEC),IC,ISPEC)
     +               + FORNOW*TTEMP(1)
              TTEMP(3) = FORNOW

C             MASS-SPECIFIC ENTROPY
              FORNOW = AMASCS(NCPOLY(IC,ISPEC),IC,ISPEC)
              DO JC = NCPOM1(IC,ISPEC),2,-1
                FORNOW = FORNOW*TTEMP(1) + AMASCS(JC,IC,ISPEC)
              ENDDO
              FORNOW = AMASCS(NCENPY(IC,ISPEC),IC,ISPEC)
     +               + FORNOW*TTEMP(1)
     +               + AMASCS(1,IC,ISPEC)*LOG(TTEMP(1))
              TTEMP(4) = FORNOW

C             MOLAR GIBBS FUNCTION
C             ACTUALLY GIBBS/(R^0 T) WITH PRESSURE TERM
              FORNOW = AMOLGB(NCPOLY(IC,ISPEC),IC,ISPEC)
              DO JC = NCPOM1(IC,ISPEC),1,-1
                FORNOW = AMOLGB(JC,IC,ISPEC)
     +                 + FORNOW*TTEMP(1)
              ENDDO
              FORNOW = AMOLGB(NCENTH(IC,ISPEC),IC,ISPEC)
     +                /TTEMP(1)
     +               - AMOLGB(NCENPY(IC,ISPEC),IC,ISPEC)
     +                *LOG(TTEMP(1))
     +               - FORNOW
              TTEMP(5) = FORNOW

              IF(ICP.EQ.1)THEN
                IF(IC.EQ.1)THEN
                  WRITE(NCREPT,'(I5,2X,I5,X,"l",X,5(1PE12.4))')
     +                  ISPEC,IC,(TTEMP(JC),JC=1,5)
                ELSE
                  WRITE(NCREPT,'(7X,I5,X,"l",X,5(1PE12.4))')
     +                        IC,(TTEMP(JC),JC=1,5)
                  DO JC = 1,5
      IF(ABS(TTOLD(JC)-TTEMP(JC)).GT.ABS(TTTOL*TTEMP(JC)))THEN
        WRITE(NCREPT,*)'Warning: INDATA: Mismatched thermo data'
        WRITE(NCREPT,'(I7,1PE12.4))')JC,TTEMP(JC)
      ENDIF
                  ENDDO
                ENDIF
              ELSE
                WRITE(NCREPT,'(7X,I5,X,"h",X,5(1PE12.4))')
     +                        IC,(TTEMP(JC),JC=1,5)
              ENDIF

              IF(ICP.EQ.2)THEN
                DO JC = 1,5
                  TTOLD(JC) = TTEMP(JC)
                ENDDO
              ENDIF

            ENDDO
          ENDDO
        ENDDO
        WRITE(NCREPT,*)

C       REACTION MECHANISM
C       STEP LIST
        WRITE(NCREPT,*)'Reaction mechanism:'
        WRITE(NCREPT,*)'  Number of steps:'
        WRITE(NCREPT,'(I5)')NSTEP
!START LR: Replacement LR:Reverting
        IF(NSTEP.NE.NSTPMX)THEN
!        IF(MAX(1,NSTEP).NE.NSTPMX)THEN
!END
          WRITE(NCREPT,*)'Warning: mismatch in number of steps:' 
          WRITE(NCREPT,9600)NSTEP,NSTPMX
        ENDIF
        DO ISTEP = 1,NSTEP

          WRITE(STRSTP,'(I5)')ISTEP
          ICOEF1 = 1
          ICOEF2 = LEN(STRSTP) + 3
          STROUT(ICOEF1:ICOEF2) = STRSTP//'   '
          DO ISPEC = 1,NRSLEN(ISTEP)
            STREAC = SPCSYM(NRSPEC(ISPEC,ISTEP))
            ICOEF1 = ICOEF2 + 1
            ICOEF2 = ICOEF1 + LEN(STREAC)
            STROUT(ICOEF1:ICOEF2) = STREAC
            ICOEF1 = ICOEF2 + 1
            ICOEF2 = ICOEF1
            STROUT(ICOEF1:ICOEF2) = '+'
          ENDDO
          IF(MBLIST(ISTEP).GT.0)THEN
            STREAC = BDYSYM(MBLIST(ISTEP))
            ICOEF1 = ICOEF2 + 1
            ICOEF2 = ICOEF1 + LEN(STREAC)
            STROUT(ICOEF1:ICOEF2) = STREAC
          ELSE
            ICOEF2 = ICOEF2 - 1
          ENDIF
          STRARO = ' => '
          IF(MGLIST(ISTEP).GT.0)STRARO = ' == '
          ICOEF1 = ICOEF2 + 1
          ICOEF2 = ICOEF1 + LEN(STRARO)
          STROUT(ICOEF1:ICOEF2) = STRARO
          DO ISPEC = 1,NSSLEN(ISTEP)
            ICP = 0
            DO JSPEC = 1,NRSLEN(ISTEP)
              IF(NRSPEC(JSPEC,ISTEP).EQ.NSSPEC(ISPEC,ISTEP))THEN
                ICP = ICP + 1
              ENDIF
            ENDDO
            ICP = NINT(DIFFMU(ISPEC,ISTEP)) + ICP
            IF(ICP.GT.0)THEN
              IF(ICP.GT.1)THEN
                WRITE(STRCOF,'(I1)')ICP
                ICOEF1 = ICOEF2 + 1
                ICOEF2 = ICOEF1 + LEN(STRCOF)
                STROUT(ICOEF1:ICOEF2) = STRCOF
              ENDIF
              STREAC = SPCSYM(NSSPEC(ISPEC,ISTEP))
              ICOEF1 = ICOEF2 + 1
              ICOEF2 = ICOEF1 + LEN(STREAC)
              STROUT(ICOEF1:ICOEF2) = STREAC
              ICOEF1 = ICOEF2 + 1
              ICOEF2 = ICOEF1
              STROUT(ICOEF1:ICOEF2) = '+'
            ENDIF
          ENDDO
          IF(MBLIST(ISTEP).GT.0)THEN
            STREAC = BDYSYM(MBLIST(ISTEP))
            ICOEF1 = ICOEF2 + 1
            ICOEF2 = ICOEF1 + LEN(STREAC)
            STROUT(ICOEF1:ICOEF2) = STREAC
          ELSE
            ICOEF2 = ICOEF2 - 1
          ENDIF
          
          WRITE(NCREPT,'(A)')STROUT(1:ICOEF2)

        ENDDO

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Reaction parameters A,n,E:'
        DO ISTEP = 1,NSTEP
          WRITE(NCREPT,'(I5,3(2X,1PE12.4))')
     +          ISTEP,(RPARAM(ICP,ISTEP),ICP=1,3)
        ENDDO

        IF(NLIND.GT.0)THEN
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Lindemann parameters:'
          DO ISTEP = 1,NSTEP
            IF(MLLIST(ISTEP).NE.0)THEN
              WRITE(NCREPT,'(I5,3(2X,1PE12.4))')
     +              ISTEP,(RCLIND(ICP,MLLIST(ISTEP)),ICP=1,3)
            ENDIF
          ENDDO
        ENDIF

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Third body efficiencies:'
        DO JSPEC = 1,NBODY
          WRITE(NCREPT,'(I5,3X,A10)')JSPEC,BDYSYM(JSPEC)
          DO ISPEC = 1,NSPEC
            WRITE(NCREPT,'(I5,3X,A10,1PE12.4)')
     +            ISPEC,SPCSYM(ISPEC),EFFY3B(ISPEC,JSPEC)
          ENDDO
        ENDDO

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'End of chemical data'
        WRITE(NCREPT,*)'--------------------'
        WRITE(NCREPT,*)
        WRITE(NCREPT,9000)CLDASH

C       =======================================================================

      ENDIF

C     =========================================================================

C     GET THE TRANSPORT DATA
C     ======================
C     RSC 06-JAN-2013 MIXTURE AVERAGED TRANSPORT
      CALL DIFFIN

      IF(FLMAVT)THEN

C       =======================================================================

        IF(IPROC.EQ.0)THEN

C         REPORT THE TRANSPORT DATA
C         =========================

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Transport data: mixture averaged transport'
          WRITE(NCREPT,*)'--------------'
          WRITE(NCREPT,*)

          WRITE(NCREPT,'(2(1PE12.4))')PDIFGB,TDIFGB
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Viscosity'
          DO ISPEC = 1,NSPEC
            WRITE(NCREPT,'(2I5)')ISPEC,NCOVIS
            WRITE(NCREPT,'(5(1PE15.7))')
     +           (VISCCO(ICOEFF,ISPEC),ICOEFF=1,NCOVIS)
          ENDDO
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Thermal conductivity'
          DO ISPEC = 1,NSPEC
            WRITE(NCREPT,'(2I5)')ISPEC,NCOCON
            WRITE(NCREPT,'(5(1PE15.7))')
     +           (CONDCO(ICOEFF,ISPEC),ICOEFF=1,NCOCON)
          ENDDO
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Binary diffusion coefficient'
          DO ISPEC = 1,NSPEC
            WRITE(NCREPT,'(2I5)')ISPEC,NCODIF
            DO JSPEC = 1,ISPEC
              WRITE(NCREPT,'(I5,5(1PE15.7))')JSPEC,
     +             (DIFFCO(ICOEFF,JSPEC,ISPEC),ICOEFF=1,NCODIF)
            ENDDO
          ENDDO
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Thermal diffusion ratio'
          DO ISPEC = 1,NSPEC
            WRITE(NCREPT,'(2I5)')ISPEC,NCOTDR
            DO JSPEC = 1,ISPEC-1
              WRITE(NCREPT,'(I5,5(1PE15.7))')JSPEC,
     +             (TDRCCO(ICOEFF,JSPEC,ISPEC),ICOEFF=1,NCOTDR)
            ENDDO
          ENDDO

C         EVALUATE TRANSPORT COEFFICIENTS FOR THE DEFAULT INITIAL STATE
          TMPLOG = LOG(TRIN/TDIFGB)

C         INITIAL DENSITY
          RTIN = ZERO
          DO ISPEC = 1,NSPEC
            RTIN = RTIN + RGSPEC(ISPEC)*YRIN(ISPEC)
          ENDDO
          RTIN = RTIN*TRIN
          DRIN = PRIN/RTIN

C         VISCOSITY FOR EACH SPECIES
          DO ISPEC = 1, NSPEC
            FORNOW = VISCCO(NCOVIS,ISPEC)
            DO ICP = NCOVM1,1,-1
              FORNOW = FORNOW*TMPLOG + VISCCO(ICP,ISPEC)
            ENDDO
            CTRANS(ISPEC,1) = EXP(FORNOW)
          ENDDO

C         COMBINATION RULE FOR VISCOSITY
          COMBO1 = ZERO
          DO ISPEC = 1, NSPEC
            COMBO2 = ZERO
            DO JSPEC = 1, NSPEC
              FORNOW = SQRT(CTRANS(ISPEC,1)/CTRANS(JSPEC,1))
              FORNOW = ONE + FORNOW*WILKO2(JSPEC,ISPEC)
              FORNOW = WILKO1(JSPEC,ISPEC)*FORNOW*FORNOW
              COMBO2 = COMBO2
     +               + YRIN(JSPEC)*OVWMOL(JSPEC)*FORNOW
            ENDDO
            FORNOW = CTRANS(ISPEC,1)/COMBO2
            COMBO1 = COMBO1
     +             + YRIN(ISPEC)*OVWMOL(ISPEC)*FORNOW

          ENDDO
          TTEMP(1) = COMBO1

C         THERMAL CONDUCTIVITY FOR EACH SPECIES
          DO ISPEC = 1, NSPEC
            FORNOW = CONDCO(NCOCON,ISPEC)
            DO ICP = NCOCM1,1,-1
              FORNOW = FORNOW*TMPLOG + CONDCO(ICP,ISPEC)
            ENDDO
            CTRANS(ISPEC,2) = EXP(FORNOW)
          ENDDO

C         COMBINATION RULE FOR CONDUCTIVITY
          COMBO1 = ZERO
          COMBO2 = ZERO
          COMBO3 = ZERO
          DO ISPEC = 1, NSPEC
            FORNOW = YRIN(ISPEC)*OVWMOL(ISPEC)
            COMBO1 = COMBO1 + FORNOW*CTRANS(ISPEC,2)
            COMBO2 = COMBO2 + FORNOW/CTRANS(ISPEC,2)
            COMBO3 = COMBO3 + FORNOW
          ENDDO
          COMBO3 = ONE/COMBO3
          COMBO1 = COMBO1*COMBO3
          COMBO2 = COMBO2*COMBO3
          TTEMP(2) = HALF*(COMBO1 + ONE/COMBO2)

C         MASS-SPECIFIC CP
          COMBO2 = ZERO
          IC = 1
          DO ISPEC = 1, NSPEC
            FORNOW = AMASCP(NCPOLY(IC,ISPEC),IC,ISPEC)
            DO JC = NCPOM1(IC,ISPEC),1,-1
              FORNOW = FORNOW*TRIN + AMASCP(JC,IC,ISPEC)
            ENDDO
            COMBO1 = FORNOW
            COMBO2 = COMBO2 + YRIN(ISPEC)*COMBO1
          ENDDO
          TTEMP(3) = COMBO2

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'         Viscosity         Conductivity'
          DO ISPEC = 1, NSPEC
            WRITE(NCREPT,'(I5,2(1PE18.7))')ISPEC,CTRANS(ISPEC,1),
     +                                           CTRANS(ISPEC,2)
          ENDDO
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)
     + '         Mix viscosity     Conductivity      Prandtl no'
          WRITE(NCREPT,'(5X,3(1PE18.7))')TTEMP(1),TTEMP(2),
     +                                   TTEMP(1)*COMBO2/TTEMP(2)


C         MASS DIFFUSIVITY FOR EACH SPECIES
          DO ISPEC = 1, NSPEC

            DO JSPEC = 1, NSPEC
              FORNOW = DIFFCO(NCOCON,JSPEC,ISPEC)
              DO ICP = NCODM1,1,-1
              FORNOW = FORNOW*TMPLOG + DIFFCO(ICP,JSPEC,ISPEC)
              ENDDO
              DTRANS(JSPEC) = EXP(FORNOW)*PDIFGB/PRIN
            ENDDO

C           COMBINATION RULE FOR MASS DIFFUSIVITY
            COMBO1 = ZERO
            COMBO2 = ZERO
            COMBO3 = ZERO
            DO JSPEC = 1, NSPEC
              FORNOW = YRIN(JSPEC) + DFCTOL
              COMBO1 = COMBO1 + FORNOW
              COMBO2 = COMBO2 + FORNOW*OVWMOL(JSPEC)/DTRANS(JSPEC)
              COMBO3 = COMBO3 + YRIN(JSPEC)*OVWMOL(JSPEC)
            ENDDO
            FORNOW = YRIN(ISPEC) + DFCTOL
            COMBO1 = COMBO1 - FORNOW
            COMBO2 = COMBO2 - FORNOW*OVWMOL(ISPEC)/DTRANS(ISPEC)
            COMBO2 = COMBO2/COMBO3
            CTRANS(ISPEC,3) = DRIN*COMBO1/COMBO2

          ENDDO

C         THERMAL DIFFUSION RATIO FOR EACH SPECIES
          TMPLOG = TRIN/TDIFGB
          DO ISPEC = 1, NSPEC

            DO JSPEC = 1, NSPEC
              FORNOW = TDRCCO(NCOCON,JSPEC,ISPEC)
              DO ICP = NCOTM1,1,-1
                FORNOW = FORNOW*TMPLOG + TDRCCO(ICP,JSPEC,ISPEC)
              ENDDO
              DTRANS(JSPEC) = FORNOW
            ENDDO

C           COMBINATION RULE FOR THERMAL DIFFUSION RATIO
            COMBO1 = ZERO
            COMBO2 = ZERO
            DO JSPEC = 1, NSPEC
              FORNOW = YRIN(JSPEC)*OVWMOL(JSPEC)
              COMBO1 = COMBO1 + FORNOW*DTRANS(JSPEC)
              COMBO2 = COMBO2 + FORNOW
            ENDDO
            CTRANS(ISPEC,4) = COMBO1/COMBO2

          ENDDO

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'         Diffusivity       Schmidt No',
     +                    '        Lewis No          Thermal diff'
          DO ISPEC = 1, NSPEC
            WRITE(NCREPT,'(I5,4(1PE18.7))')ISPEC,CTRANS(ISPEC,3),
     +                                  TTEMP(1)/CTRANS(ISPEC,3),
     +                        TTEMP(2)/(CTRANS(ISPEC,3)*TTEMP(3)),
     +                                           CTRANS(ISPEC,4)
          ENDDO

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'End of transport data'
          WRITE(NCREPT,*)'---------------------'
          WRITE(NCREPT,*)
          WRITE(NCREPT,9000)CLDASH

        ENDIF

C       =======================================================================

      ELSE

C       =======================================================================

        IF(IPROC.EQ.0)THEN

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Transport data: constant Lewis numbers'
          WRITE(NCREPT,*)'--------------'
          WRITE(NCREPT,*)

          WRITE(NCREPT,*)' constants A,r; ref T; Prandtl no.:'
          WRITE(NCREPT,'(4(1PE12.4))')ALAMDC,RLAMDA,TLAMDA,PRANTL
          WRITE(NCREPT,*)

          WRITE(NCREPT,*)'Spec.No.  Symbol         Lewis no.'
          DO ISPEC = 1,NSPEC
            WRITE(NCREPT,'(I8,3X,A10,3X,1PE12.4)')
     +            ISPEC,SPCSYM(ISPEC),CLEWIS(ISPEC)
          ENDDO

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'End of transport data'
          WRITE(NCREPT,*)'---------------------'
          WRITE(NCREPT,*)
          WRITE(NCREPT,9000)CLDASH

        ENDIF

C       =======================================================================

      ENDIF
C     MIXTURE AVERAGED TRANSPORT

C     =========================================================================

C     RADIATION TREATMENT
      CALL RADCIN

      IF(FLRADN)THEN

C       =======================================================================

        IF(IPROC.EQ.0)THEN

C         REPORT THE RADIATION DATA
C         =========================
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Radiation data'
          WRITE(NCREPT,*)'--------------'
          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Planck mean absorption coefficient data'
          DO ISPEC = 1, NSPRAD
            WRITE(NCREPT,'(2I6)')NSPRID(ISPEC),NKPRAD(ISPEC)
            WRITE(NCREPT,'(5(1PE15.7))')(AKPRAD(ICP,ISPEC),ICP=1,5)
            WRITE(NCREPT,'(5(1PE15.7))')(AKPRAD(ICP,ISPEC),
     +                                   ICP=6,NKPRAD(ISPEC))
          ENDDO

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Reference ambient temperature:'
          WRITE(NCREPT,'(1PE15.7)')TREFRN

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'Planck mean absorption coefficients:'
          DO ISPEC = 1, NSPRAD
            FORNOW = AKPRAD(NKPRAD(ISPEC),JSPEC)
            DO ICP = NKPRM1(ISPEC),1,-1
              FORNOW = FORNOW*TREFRN + AKPRAD(ICP,ISPEC)
            ENDDO
            JSPEC = NSPRID(ISPEC)
          WRITE(NCREPT,'(I6,5X,A10,1PE15.7)')JSPEC,SPCSYM(JSPEC),FORNOW
          ENDDO

          WRITE(NCREPT,*)
          WRITE(NCREPT,*)'End of radiation data'
          WRITE(NCREPT,*)'---------------------'
          WRITE(NCREPT,*)
          WRITE(NCREPT,9000)CLDASH

        ENDIF

C       =======================================================================

      ENDIF

C     =========================================================================

C     END OF REPORTING OF INITIAL DATA FOR NOW
C     NOTE THAT THE REPORT FILE IS LEFT OPEN
C     FOR FURTHER REPORTING FROM OTHER INITIALISATION ROUTINES
C     AND IS CLOSED BEFORE EXIT FROM THIS SUBROUTINE

C     =========================================================================

!START LR: Removal LR: Reverting
C     INITIALISE ONLINE STATISTICS
C     ============================

C     STATISTICS ON ONE PROCESSOR ONLY
      IF(IPROC.EQ.0)THEN

C       INITIALISE THE STATISTICS FILE
        OPEN(UNIT=NCSTAT,FILE=FNSTAT,STATUS='NEW',FORM='FORMATTED')
        ENDFILE(NCSTAT)
        CLOSE(NCSTAT)

      ENDIF

C     ZERO THE STATISTICS STORAGE COUNTER
      ITSTAT = 0

C     =========================================================================
!END

C     CHECK AND INITIALISE DUMP FILES
C     -------------------------------
#ifndef HDF5
!START LR: Removal LR: Reverting
      INQUIRE(FILE=FNDMPO(1),EXIST=FXDUMP)
      IF(.NOT.FXDUMP)THEN
        IF(NDOFMT.EQ.0)THEN
        OPEN(UNIT=NCDMPO,FILE=FNDMPO(1),STATUS='NEW',FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=NCDMPO,FILE=FNDMPO(1),STATUS='NEW',FORM='FORMATTED')
        ENDIF
        CLOSE(NCDMPO)
      ENDIF
      INQUIRE(FILE=FNDMPO(2),EXIST=FXDUMP)
      IF(.NOT.FXDUMP)THEN
        IF(NDOFMT.EQ.0)THEN
        OPEN(UNIT=NCDMPO,FILE=FNDMPO(2),STATUS='NEW',FORM='UNFORMATTED')
        ELSE
          OPEN(UNIT=NCDMPO,FILE=FNDMPO(2),STATUS='NEW',FORM='FORMATTED')
        ENDIF
        CLOSE(NCDMPO)
      ENDIF
!END
#else
!START LR: Replacement
!      CALL CREATE_H5DUMP_FILES

C     1. BUILD RESTART FILENAMES
C     EACH FILE CONSISTS OF FOLLOWING ELEMENTS:

C     TIME/STEP ELEMENT
C     -----------------
C     SET FNDMPO(1)
      IF(NCDMPI.EQ.0)THEN
         WRITE(PNTIME,'(I6.6)')
     $        FLOOR(REAL(ITIME)/(REAL(50*NTSTAT)))
      ELSE
         WRITE(PNTIME,'(I6.6)')
     $        CEILING(REAL(ITIME)/(REAL(50*NTSTAT)))
      ENDIF

C     EXTENSION
C     ---------
      IF(ABS(NDOFMT).EQ.1)THEN
         PNEXT = PNXBIN
      ELSE IF (ABS(NDOFMT).EQ.2)THEN
         PNEXT = PNXDAT
      ELSE
         PNEXT = PNXHDF
      ENDIF

C     PROC/COLLECTIVE
C     ---------------
      IF(NDOFMT.EQ.0)THEN
C     FNDMPO(1) = PNDMPI//PNTIME//PNEXT
!         WRITE(FNDMPO(1),'((1A4),(I6.6),(1A3))') PNDMPI,IDFLAG,PNXHDF   !LR: Commented out
         IDFLAG = 0                                                     !LR: New
         WRITE(FNDMPO(1),'((1A11),(I1),(1A3))') PNDMPI,IDFLAG,PNXHDF   !LR: New
         IDFLAG = 1                                                     !LR: New
         WRITE(FNDMPO(2),'((1A11),(I1),(1A3))') PNDMPI,IDFLAG,PNXHDF   !LR: New
      ELSE
         WRITE(PNPROC,'(I6.6)') IPROC
         FNDMPO(1) = PNDMPI//PNTIME//'.'//PNPROC//PNEXT
      ENDIF
C     USE FORMATTED ONLY FOR CASE 2 (STEP) AND -2 (TIME)
C     IN ALL OTHER CASES INITIALISE UNFORMATTED FILES

      IF(ABS(NDOFMT).GE.1)THEN
         INQUIRE(FILE=FNDMPO(1),EXIST=FXDUMP)
         IF(.NOT.FXDUMP)THEN
            IF(ABS(NDOFMT).EQ.1)THEN
               OPEN(UNIT=NCDMPO,FILE=FNDMPO(1),
     $              STATUS='NEW',FORM='UNFORMATTED')
            ELSE IF(ABS(NDOFMT).EQ.2)THEN
               OPEN(UNIT=NCDMPO,FILE=FNDMPO(1),
     $              STATUS='NEW',FORM='FORMATTED')
            ENDIF
            CLOSE(NCDMPO)
         ENDIF
      ENDIF

C     =========================================================================

C     BUILD THE FILENAMES FOR INLET/INITIAL TURBULENT VELOCITY FIELD      
      WRITE(PNPROC,'(I6.6)')IPROC

C     INLET RESTART FILE NAME
      IF(INTURB.EQ.2)THEN
         IF(NDIFMT.EQ.0)THEN

            WRITE(PNPROC,'(I6.6)') (FLOOR(REAL(IXPROC/NYPROC)))
            FNTCXL = PNTCXL//PNPROC//PNXHDF
C     WRITE(FNTCXL,'((1A4),(I6.6),(1A3))') PNTCXL,IDFLAG,PNXHDF

         ELSE

            WRITE(PNPROC,'(I6.6)') (MOD(IXPROC,NYPROC) +
     $           IYPROC*NXPROC+IZPROC*NXPROC*NYPROC)
            FNTCXL = PNTCXL//PNPROC//PNXDAT

         ENDIF
      ENDIF

C     INLET RESTART FILE NAME
      IF(INTURB.EQ.3)THEN
         IF(NDIFMT.EQ.0)THEN

C     READ FROM A SINGLE FILE
            WRITE(PNPROC,'(I6.6)') (FLOOR(REAL(IXPROC/NXPROC)))
            FNTCXL = PNTCXL//PNPROC//PNXHDF

         ELSE

            WRITE(PNPROC,'(I6.6)') (MOD(IXPROC,NXPROC) +
     $           IYPROC*NXPROC+IZPROC*NXPROC*NYPROC)
            FNTCXL = PNTCXL//PNPROC//PNXDAT

         ENDIF
      ENDIF

      IF(NGBCXL.EQ.NSBCI3.AND.NXLPRM(1).EQ.3)THEN

         IF(NDIFMT.EQ.0) THEN
            WRITE(FNTIXL,'((1A4),(I6.6),(1A3))') PNTIXL,IDFLAG,PNXHDF            
         ELSE
            FNTIXL = PNTIXL//PNPROC//PNXDAT
         ENDIF

      ENDIF
!END
#endif

C     ==========================================================================

C     INITIALISE ADDITIONAL PARAMETERS
C     ================================

C     SET VALUE OF PI
C     ---------------
      PI = FOUR*ATAN(ONE)

C     SET VALUE OF LN(10)
C     -------------------
      CLNTEN = LOG(1.0D1)

C     TIME STEPPING
C     -------------
      ETIME = ZERO
      NTIME2 = NTIME1+NTIME-1
!START LR: Replacement LR:Reverting
      ITIME = 0
!      ITIME = NTIME1 - 1
!END
!START LR: Addition LR:Reverting
!C     COUNTERS FOR TIME ACCURATE OUTPUT
!C     ---------------------------------
!      OTIME = 0
!      STIME = 0.0
!END
 
C     =========================================================================

C     INITIALISE BOUNDARIES
C     =====================
      CALL BCINIT

C     =========================================================================

C     INITIALISE THE VELOCITY FIELD
C     =============================
C     NOTE: ALL VARIABLES ARE INITIALISED ON STANDARD SIZE DOMAIN ONLY

C     PRE-INITIALISE THE VELOCITY FIELD TO ZERO
C     ---------------------------------
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL 

            URUN(IC,JC,KC) = ZERO
            VRUN(IC,JC,KC) = ZERO
            WRUN(IC,JC,KC) = ZERO

          ENDDO
        ENDDO
      ENDDO

C     -------------------------------------------------------------------------

C     INITIAL TURBULENCE GENERATOR
C     ----------------------------

!START LR: Replacement LR: Reverting
C     GENERATE FRESH INITIAL TURBULENCE
      IF(INTURB.EQ.1)CALL TURBIN

C     COPY TURBULENT INFLOW VELOCITY FIELD INTO THE DOMAIN 
      IF(INTURB.EQ.2)THEN

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL

              URUN(IC,JC,KC) = STORE4(IC,JC,KC)
              VRUN(IC,JC,KC) = STORE5(IC,JC,KC)
              WRUN(IC,JC,KC) = STORE6(IC,JC,KC)

            ENDDO
          ENDDO
        ENDDO

      ENDIF

C     -------------------------------------------------------------------------

!C     GENERATE FRESH INITIAL TURBULENCE
!     
!      IF(INTURB.EQ.1) THEN 
!     
!         CALL TURBIN
!     
!C     REPORT OUTPUT
!         IF(IPROC.EQ.0)THEN
!     
!            OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='OLD',
!     $           FORM='FORMATTED',POSITION='APPEND')
!            WRITE(NCREPT,9620)
!            WRITE(NCREPT,*)     
!     
!         ENDIF        
!     
!      ELSE IF (INTURB.GE.2)THEN
!     
!         CALL HC_TCXL(FNTCXL)
!
!         DO KC = KSTAL,KSTOL
!            DO JC = JSTAL,JSTOL
!               DO IC = ISTAL,ISTOL 
!     
!                  URUN(IC,JC,KC) = STORE4(IC,JC,KC)
!                  VRUN(IC,JC,KC) = STORE5(IC,JC,KC)
!                  WRUN(IC,JC,KC) = STORE6(IC,JC,KC)
!     
!               ENDDO
!            ENDDO
!         ENDDO        
!     
!      ELSE IF (INTURB.EQ.4) THEN 
!
!         CALL TGVRTX
!
!      ENDIF
!          
!C
!-------------------------------------------------------------------------
!END

C     ADD ON THE DEFAULT MEAN VELOCITY
C     --------------------------------
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL

            URUN(IC,JC,KC) = URIN + URUN(IC,JC,KC)
            VRUN(IC,JC,KC) = VRIN + VRUN(IC,JC,KC)
            WRUN(IC,JC,KC) = WRIN + WRUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

C     =========================================================================

C     INITIALISE THE SCALAR FIELD
C     ===========================
C     NOTE: ALL VARIABLES ARE INITIALISED ON STANDARD SIZE DOMAIN ONLY

C     PRE-INITIALISE SPECIES MASS FRACTIONS TO DEFAULT VALUES
C     -------------------------------------
      DO ISPEC = 1,NSPEC

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL 

              YRUN(IC,JC,KC,ISPEC) = YRIN(ISPEC)

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     PRE-INITIALISE PRESSURE AND TEMPERATURE TO DEFAULT VALUES
C     ---------------------------------------
C     BIGGER SIZE ARRAYS
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB

!START LR: Replacement LR: Reverting
            PRUN(IC,JC,KC) = PRIN
            TRUN(IC,JC,KC) = TRIN
!C     For a convecting vortex, change the pressure and density
!C     distribution according to analytical solution in vlamin
!            IF(INTURB.GE.0) PRUN(IC,JC,KC) = PRIN 
!            TRUN(IC,JC,KC) = TRIN 
!END

          ENDDO
        ENDDO
      ENDDO

C     INITIAL DENSITY
C     ---------------
      RTIN = ZERO
      DO ISPEC = 1,NSPEC
        RTIN = RTIN + RGSPEC(ISPEC)*YRIN(ISPEC)
      ENDDO
      RTIN = RTIN*TRIN
      DRIN = PRIN/RTIN

C     INITIAL INTERNAL ENERGY
C     -----------------------
      ERIN = ZERO
      DO ISPEC = 1,NSPEC

C       TEMPERATURE INTERVAL INDEX
        ITINT = 1
1010    CONTINUE
          IF(TRIN.GT.TINTHI(ITINT,ISPEC))THEN
            IF(ITINT.LT.NTINT(ISPEC))THEN
              ITINT = ITINT + 1
              GOTO 1010
            ENDIF
          ENDIF
C       END OF LOOP 1010

        FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
        DO ICP = NCPOM1(ITINT,ISPEC),1,-1
          FORNOW = FORNOW*TRIN + AMASCH(ICP,ITINT,ISPEC)
        ENDDO
        FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC) + FORNOW*TRIN
        ERIN = ERIN + FORNOW*YRIN(ISPEC)

      ENDDO
      FORNOW = URIN*URIN + VRIN*VRIN + WRIN*WRIN
      ERIN = ERIN - RTIN + HALF*FORNOW

C     -------------------------------------------------------------------------

C     INITIAL FLAME GENERATOR
C     -----------------------
      IF(INFLAM.EQ.1)THEN
      
C       GENERATE AN INITIAL THERMOCHEMICAL FIELD
C       TEMPERATURE, PRESSURE AND MASS FRACTIONS
C       MODIFY VELOCITIES AS REQUIRED
        CALL FLAMIN

      ENDIF

C     -------------------------------------------------------------------------

C     DENSITY FIELD
C     -------------
C     USE STORE1 AS AN ACCUMULATOR FOR MIXTURE GAS CONSTANT
C     INITIALISE TO ZERO
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL 

            STORE1(IC,JC,KC) = ZERO

          ENDDO
        ENDDO
      ENDDO

C     MIXTURE GAS CONSTANT 
      DO ISPEC = 1,NSPEC

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL 

               STORE1(IC,JC,KC) = STORE1(IC,JC,KC)
     +                          + RGSPEC(ISPEC)*YRUN(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     DENSITY
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL 

            STORE1(IC,JC,KC) = STORE1(IC,JC,KC)*TRUN(IC,JC,KC)
            DRUN(IC,JC,KC) = PRUN(IC,JC,KC)/STORE1(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

!START LR: Addition LR: Reverting
!      IF (INFLAM.GE.2) THEN 
!
!C     Generate laminar vortex core        
!     
!         CALL VLAMIN
!     
!      ENDIF
!END

C                                                             STORE1 = T*MIX RG
C     -------------------------------------------------------------------------

C     INTERNAL ENERGY FIELD
C     ---------------------
C     INITIALISE TEMPERATURE INTERVAL INDEX
C     FOR ALL SPECIES LOCATE TEMPERATURE IN AN INTERVAL
C     BIGGER SIZE ARRAY
      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB 

            DO IINDEX = 1,NINTMX
              ITNDEX(IC,JC,KC,IINDEX) = 0
            ENDDO

            DO ISPEC = 1,NSPEC

              ITINT = 1
1000          CONTINUE
                IF(TRUN(IC,JC,KC).GT.TINTHI(ITINT,ISPEC))THEN
                  IF(ITINT.LT.NTINT(ISPEC))THEN
                    ITINT = ITINT + 1
                    GOTO 1000
                  ENDIF
                ENDIF
C             END OF LOOP 1000

C             SET THE TEMPERATURE INDEX
              IINDEX = 1 + (ISPEC-1)/NSPIMX
              IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
              ITNDEX(IC,JC,KC,IINDEX) = ITNDEX(IC,JC,KC,IINDEX)
     +                          +(ITINT-1)*NTBASE**IPOWER

            ENDDO

          ENDDO
        ENDDO
      ENDDO

C     PRE-INITIALISE INTERNAL ENERGY TO ZERO
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL 

            ERUN(IC,JC,KC) = ZERO

          ENDDO
        ENDDO
      ENDDO

C     INITIALISE INTERNAL ENERGY
      DO ISPEC = 1,NSPEC

C       TEMPERATURE INTERVAL INDEXING 
        IINDEX = 1 + (ISPEC-1)/NSPIMX
        IPOWER = ISPEC - (IINDEX-1)*NSPIMX - 1
        ICOEF2 = NTBASE**IPOWER
        ICOEF1 = ICOEF2*NTBASE

C       USE THERMOCHEMICAL DATA
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL 

              ITINT = 1 + MOD(ITNDEX(IC,JC,KC,IINDEX),ICOEF1)/ICOEF2
              FORNOW = AMASCH(NCPOLY(ITINT,ISPEC),ITINT,ISPEC)
              DO ICP = NCPOM1(ITINT,ISPEC),1,-1
                FORNOW = FORNOW*TRUN(IC,JC,KC) + AMASCH(ICP,ITINT,ISPEC)
              ENDDO
              FORNOW = AMASCH(NCENTH(ITINT,ISPEC),ITINT,ISPEC)
     +               + FORNOW*TRUN(IC,JC,KC)

              ERUN(IC,JC,KC) = ERUN(IC,JC,KC)
     +                       + FORNOW*YRUN(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO

      ENDDO

C     FINALISE INTERNAL ENERGY
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL 

            FORNOW = URUN(IC,JC,KC)*URUN(IC,JC,KC)
     +             + VRUN(IC,JC,KC)*VRUN(IC,JC,KC)
     +             + WRUN(IC,JC,KC)*WRUN(IC,JC,KC)

            ERUN(IC,JC,KC) = ERUN(IC,JC,KC)
     +                     - STORE1(IC,JC,KC)
     +                     + HALF*FORNOW

          ENDDO
        ENDDO
      ENDDO
C                                                              ALL STORES CLEAR
C     =========================================================================

C     CONVERT VARIABLES TO CONSERVATIVE FORM
C     ======================================
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL 

            URUN(IC,JC,KC) = DRUN(IC,JC,KC)*URUN(IC,JC,KC)
            VRUN(IC,JC,KC) = DRUN(IC,JC,KC)*VRUN(IC,JC,KC)
            WRUN(IC,JC,KC) = DRUN(IC,JC,KC)*WRUN(IC,JC,KC)

            ERUN(IC,JC,KC) = DRUN(IC,JC,KC)*ERUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

      DO ISPEC = 1,NSPEC
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL 

              YRUN(IC,JC,KC,ISPEC) = DRUN(IC,JC,KC)*YRUN(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

C     =========================================================================

C     EXECUTE REQUIRED STARTUP
C     ========================

C     COLD START SWITCH
C     -----------------
      IF(NCDMPI.EQ.1)THEN

C       =======================================================================
C       WARM START
C       =======================================================================

C       -----------------------------------------------------------------------
C       THIS BLOCK MAY BE MODIFIED AS REQUIRED
C       TO BLEND INITIAL VELOCITY AND SCALAR FIELDS
C       WITH PREVIOUSLY DUMPED DATA
C       -----------------------------------------------------------------------

!START LR: Removal LR: Reverting
        !Vishnu: Prevent overwriting TSTEP
        TSOLD=TSTEP
!END

C       RESTART FROM FULL DUMP FILES
C       ----------------------------
C       READ THE DATA FROM DUMP INPUT FILE 1
C       NOTE THAT URUN,VRUN,WRUN,ERUN AND YRUN ARE ALL IN CONSERVATIVE FORM
C       RSC 11-JUL-2009 ADD A DUMP FORMAT SWITCH
#ifndef HDF5
        IF(NDIFMT.EQ.0)THEN

C         UNFORMATTED DUMP INPUT
          OPEN(UNIT=NCDMPI,FILE=FNDMPO(1),STATUS='OLD',
     +         FORM='UNFORMATTED')
          READ(NCDMPI)NXDMAX,NYDMAX,NZDMAX,NDSPEC,
     +                DRUN,URUN,VRUN,WRUN,ERUN,YRUN,
     +                ETIME,TSTEP,ERROLD,ERRLDR

C         SIZE ERROR CHECK
          IF(NXDMAX.NE.NXNODE)WRITE(6,*)'Dump input size error: x'
          IF(NYDMAX.NE.NYNODE)WRITE(6,*)'Dump input size error: y'
          IF(NZDMAX.NE.NZNODE)WRITE(6,*)'Dump input size error: z'
          IF(NDSPEC.NE.NSPEC)WRITE(6,*)'Dump input size error: species'

          CLOSE(NCDMPI)

        ELSE

C         FORMATTED DUMP INPUT
          OPEN(UNIT=NCDMPI,FILE=FNDMPO(1),STATUS='OLD',FORM='FORMATTED')
          READ(NCDMPI,*)NXDMAX,NYDMAX,NZDMAX,NDSPEC

C         SIZE ERROR CHECK
          IF(NXDMAX.NE.NXNODE)WRITE(6,*)'Dump input size error: x'
          IF(NYDMAX.NE.NYNODE)WRITE(6,*)'Dump input size error: y'
          IF(NZDMAX.NE.NZNODE)WRITE(6,*)'Dump input size error: z'
          IF(NDSPEC.NE.NSPEC)WRITE(6,*)'Dump input size error: species'

          DO KC = 1, NZNODE
            DO JC = 1, NYNODE
              DO IC = 1, NXNODE
                READ(NCDMPI,*)DRUN(IC,JC,KC),
     +                     URUN(IC,JC,KC),VRUN(IC,JC,KC),WRUN(IC,JC,KC),
     +                        ERUN(IC,JC,KC),
     +                       (YRUN(IC,JC,KC,ISPEC),ISPEC=1,NSPEC)
              ENDDO
            ENDDO
          ENDDO

          READ(NCDMPI,*)ETIME,TSTEP,ERROLD,ERRLDR

          CLOSE(NCDMPI)

        ENDIF
#else
!START LR: Replacement
!        CALL READ_H5DUMP_FILES

C       COLLECTIVE HDF DUMP INPUT            
        IF(NDIFMT.EQ.0)THEN

            CALL HC_RESI(FNDMPO(1))

C       UNFORMATTED DUMP INPUT            
        ELSE IF(NDIFMT.EQ.1)THEN

          OPEN(UNIT=NCDMPO,FILE=FNDMPO(1),STATUS='OLD',
     +         FORM='UNFORMATTED')
          READ(NCDMPO)NXDMAX,NYDMAX,NZDMAX,NDSPEC,
     +         DRUN,URUN,VRUN,WRUN,ERUN,YRUN,
     +         ETIME,TSTEP,ERROLD,ERRLDR

C         SIZE ERROR CHECK
          IF(NXDMAX.NE.NXNODE)WRITE(6,*)'Dump input size error: x'
          IF(NYDMAX.NE.NYNODE)WRITE(6,*)'Dump input size error: y'
          IF(NZDMAX.NE.NZNODE)WRITE(6,*)'Dump input size error: z'
          IF(NDSPEC.NE.NSPEC)WRITE(6,*)'Dump input size error: species'

          CLOSE(NCDMPO)

C       FORMATTED DUMP INPUT
        ELSE IF(NDIFMT.EQ.2)THEN

          OPEN(UNIT=NCDMPO,FILE=FNDMPO(1),STATUS='OLD',FORM='FORMATTED')
          READ(NCDMPO,*)NXDMAX,NYDMAX,NZDMAX,NDSPEC

C         SIZE ERROR CHECK
          IF(NXDMAX.NE.NXNODE)WRITE(6,*)'Dump input size error: x'
          IF(NYDMAX.NE.NYNODE)WRITE(6,*)'Dump input size error: y'
          IF(NZDMAX.NE.NZNODE)WRITE(6,*)'Dump input size error: z'
          IF(NDSPEC.NE.NSPEC)WRITE(6,*)'Dump input size error: species'

          DO KC = 1, NZNODE
             DO JC = 1, NYNODE
                DO IC = 1, NXNODE
                   READ(NCDMPI,*)DRUN(IC,JC,KC),
     +                  URUN(IC,JC,KC),VRUN(IC,JC,KC),WRUN(IC,JC,KC),
     +                  ERUN(IC,JC,KC),
     +                  (YRUN(IC,JC,KC,ISPEC),ISPEC=1,NSPEC)
                ENDDO
             ENDDO
          ENDDO

          READ(NCDMPO,*)ETIME,TSTEP,ERROLD,ERRLDR

          CLOSE(NCDMPO)

C       PARALLEL INDEPENDENT HDF INPUT
        ELSE IF (NDIFMT.EQ.3)THEN

          CALL HP_RESI(FNDMPO(1))

        ENDIF
!END
#endif

        !Vishnu: TSTEP not overwritten
        if(NSTPSW==0) TSTEP=TSOLD

!START LR: Addition LR: Reverting
!C     ASSIGN ELAPSED RESTART TIME TO INITIAL TIME
!C     -------------------------------------------
!       
!       STIME = ETIME
!       OTIME = FLOOR(ETIME/(MICRO*NTDUMP))
!       OTIME = (OTIME+1)*NTDUMP
!END

C       =======================================================================
C       WARM START COMPLETE
C       =======================================================================

!START LR: Addition LR: Reverting
!C     REPORT OUTPUT
!         IF(IPROC.EQ.0)THEN
!
!            OPEN(UNIT=NCREPT,FILE=FNREPT,STATUS='OLD',
!     $           FORM='FORMATTED',POSITION='APPEND')
!            WRITE(NCREPT,*) 'Restart completed successfully. 
!     $           Elapsed time: ',ETIME
!            WRITE(NCREPT,*) ''
!
!         ENDIF
!END


      ENDIF

C     ==========================================================================

C     INITIALISE TIME-STEPPING RHS TERMS
C     ==================================
C     BIGGER SIZE ARRAYS

C     PRE-INITIALISE TO DEFAULT VALUES
      DURIN = DRIN*URIN
      DVRIN = DRIN*VRIN
      DWRIN = DRIN*WRIN
      DERIN = DRIN*ERIN
      DO ISPEC = 1,NSPEC
        DYRIN(ISPEC) = DRIN*YRIN(ISPEC)
      ENDDO

      DO KC = KSTAB,KSTOB
        DO JC = JSTAB,JSTOB
          DO IC = ISTAB,ISTOB 

            DRHS(IC,JC,KC) = DRIN
            URHS(IC,JC,KC) = DURIN
            VRHS(IC,JC,KC) = DVRIN
            WRHS(IC,JC,KC) = DWRIN
            ERHS(IC,JC,KC) = DERIN

          ENDDO
        ENDDO
      ENDDO

      DO ISPEC = 1,NSPEC
        DO KC = KSTAB,KSTOB
          DO JC = JSTAB,JSTOB
            DO IC = ISTAB,ISTOB 

              YRHS(IC,JC,KC,ISPEC) = DYRIN(ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDDO


C     NOTE THAT THESE ARE BIGGER SIZE ARRAYS
C     BUT THAT NO HALO DATA IS YET AVAILABLE
C     THEREFORE ONLY THE STANDARD DOMAIN SIZE IS INITIALISED
      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          DO IC = ISTAL,ISTOL 

            DRHS(IC,JC,KC) = DRUN(IC,JC,KC)
            URHS(IC,JC,KC) = URUN(IC,JC,KC)
            VRHS(IC,JC,KC) = VRUN(IC,JC,KC)
            WRHS(IC,JC,KC) = WRUN(IC,JC,KC)
            ERHS(IC,JC,KC) = ERUN(IC,JC,KC)

          ENDDO
        ENDDO
      ENDDO

      DO ISPEC = 1,NSPEC
        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            DO IC = ISTAL,ISTOL 

              YRHS(IC,JC,KC,ISPEC) = YRUN(IC,JC,KC,ISPEC)

            ENDDO
          ENDDO
        ENDDO
      ENDDO

C     =========================================================================

C     INITIAL PARALLEL DATA TRANSFER
C     ==============================
C     INITIALISES TIME-STEPPING RHS TERMS TO BIGGER DOMAIN SIZE
      CALL PARFER

C     =========================================================================

C     INITIAL TEMPERATURE AND PRESSURE
C     ================================
      CALL TEMPIN

C     =========================================================================

C     INITIALISE TIME STEPPING
C     ========================
      CALL DTINIT

C     ==========================================================================

C     INITIALISE SPATIAL DIFFERENTIATORS
C     ==================================
      CALL DFINIT

C     ==========================================================================

C     INITIALISE SPATIAL FILTERING
C     ============================
      CALL FLINIT

C     ==========================================================================

C     CLOSE THE REPORT FILE
C     =====================
      IF(IPROC.EQ.0)THEN

        WRITE(NCREPT,9000)CLDASH
        WRITE(NCREPT,*)

        CLOSE(NCREPT)

      ENDIF

C     ==========================================================================


      RETURN

C     FORMATS FOR OUTPUT OF INITIAL DATA
9000  FORMAT(A79)
9010  FORMAT('**',75X,'**')
9020  FORMAT('**',18X,'DIRECT NUMERICAL SIMULATION CODE SENGA2',
     +            18X,'**')
9030  FORMAT('**',31X,'Stewart Cant',32X,'**')
9040  FORMAT('**',16X,'Cambridge University Engineering Department',
     +            16X,'**')
9200  FORMAT(12I6)
9250  FORMAT(6I7)
9270  FORMAT(2I5,5(1PE12.4))
9300  FORMAT(1PE12.4,6I8)
9400  FORMAT(6(1PE12.4))
9450  FORMAT(I7,1PE15.7)
9600  FORMAT(' Actual (NSPEC) = ',I5,2X,'Max (NSPCMX) = ',I5)
!START LR: Addition LR:Reverting
!9610  FORMAT('Turbulent field from restart.')
!9620  FORMAT('Turbulent field generated.')
!9630  FORMAT('Turbulent field read from: ',A) 
!END

      END
