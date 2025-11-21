      SUBROUTINE ESPINI
 
C     *************************************************************************
C
C     ESPINI
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  DEPARTMENT OF MECHANICAL ENGINEERING, UMIST
C
C     CHANGE RECORD
C     -------------
C     11-APR-2004:  CREATED
C     09-JUN-2008:  RSC BUG FIX TURBULENCE INTEGRAL SCALING
C     09-JUN-2008:  RSC ENHANCEMENTS TO LENGTH SCALE AND DISSIPATION DATA
C
C     DESCRIPTION
C     -----------
C     INITIALISES CONSTANTS FOR THE TURBULENT ENERGY SPECTRUM FUNCTION ESPECT
C
C     REFERENCES
C     ----------
C     1) BATCHELOR, G.K., TOWNSEND, A.A.: JFM 88(4), 685-709, 1948.
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
      INCLUDE 'com_espect.h'
C     -------------------------------------------------------------------------


C     PARAMETERS
C     ==========
C     RSC 09-JUN-2008 ENHANCEMENTS
      DOUBLE PRECISION TEN
      PARAMETER(TEN = 1.0D1)


C     FUNCTIONS
C     =========
      DOUBLE PRECISION ESPECT,ESPOVK,ESPKSQ
      EXTERNAL ESPECT,ESPOVK,ESPKSQ


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TENRGY,TINSCL,TDISIP
      DOUBLE PRECISION SPECLO,SPECHI
      DOUBLE PRECISION ARGMNT
C     RSC 09-JUN-2008 ENHANCEMENTS
      DOUBLE PRECISION VISCOS,VISKIN,TKOLMO,TAYLOR,RTIN,DVRIN
      INTEGER IMODES,NMODES
C     RSC 09-JUN-2008 ENHANCEMENTS
      INTEGER ISPEC


C     BEGIN
C     ===== 
C     =========================================================================

C     SET THE SPECTRUM PARAMETERS
C     ---------------------------
      CONST0 = SPARAM(1)
      CK0 = SPARAM(2)
      OVK0 = ONE/CK0
      COVK0 = CONST0/CK0

C     =========================================================================

C     CHECK THE ENERGY SPECTRUM
C     -------------------------
      IF(IPROC.EQ.0)THEN

        WRITE(NCREPT,*)
        WRITE(NCREPT,*)'Energy spectrum function data'
        WRITE(NCREPT,*)'-----------------------------'
        WRITE(NCREPT,*)

C       SET MAX NO OF MODES FOR SPECTRUM EVALUATION
        NMODES = NGZMAX/2

C       EVALUATE THE SPECTRUM
        WRITE(NCREPT,*)'spectrum function:'
        WRITE(NCREPT,*)'mode       k          E(k)'
        DO IMODES = 0, NMODES
C         RSC 09-JUN-2008 BUG FIX
C          ARGMNT = TWO*PI*REAL(IMODES)
          ARGMNT = REAL(IMODES)
          WRITE(NCREPT,'(I5,2(1PE12.4))')IMODES,ARGMNT,ESPECT(ARGMNT)
        ENDDO
        WRITE(NCREPT,*)

C       INTEGRATE THE SPECTRUM
C       LIMITS
        SPECLO = ZERO
C         RSC 09-JUN-2008 BUG FIX
C        SPECHI = TWO*PI*REAL(NMODES)
        SPECHI = REAL(NMODES)

C       INTEGRATE FOR THE ENERGY
        CALL INTEGF(ESPECT,SPECLO,SPECHI,TENRGY)
C       INTEGRATE FOR THE LENGTH SCALE
        CALL INTEGF(ESPOVK,SPECLO,SPECHI,TINSCL)
C       INTEGRATE FOR THE DISSIPATION
        CALL INTEGF(ESPKSQ,SPECLO,SPECHI,TDISIP)

C       WRITE OUT THE RESULTS
C       WRITE OUT THE RAW INTEGRALS
C       APPLIES TO SPHERICAL UNIT VOLUME IN K-SPACE
        WRITE(NCREPT,*)'raw integrals:'
        WRITE(NCREPT,*)'energy, length scale, dissipation'
        WRITE(NCREPT,'(3(1PE12.4))')TENRGY,TINSCL,TDISIP
        WRITE(NCREPT,*)

C       CONVERT TO LENGTH SCALE AND DISSIPATION
C       RSC 09-JUN-2008 BUG FIX
C        TINSCL = TINSCL*THREE*PI/TENRGY/TWO/TWO
C        TDISIP = TWO*TDISIP
        TINSCL = TINSCL*THREE/TENRGY/TWO/TWO/TWO
        TDISIP = TWO*TWO*TWO*PI*PI*TDISIP
        WRITE(NCREPT,*)'scaled values:'
        WRITE(NCREPT,*)'energy, length scale, dissipation/viscosity:'
        WRITE(NCREPT,'(3(1PE12.4))')TENRGY,TINSCL,TDISIP
        WRITE(NCREPT,*)

C       RSC 09-JUN-2008 ENHANCEMENTS
C       CONVERT TO DIMENSIONAL LENGTH SCALE AND DISSIPATION/VISCOSITY
        TINSCL = TINSCL*XGDLEN
        TDISIP = TDISIP/(XGDLEN*XGDLEN)
        WRITE(NCREPT,*)'dimensional values:'
        WRITE(NCREPT,*)'energy, length scale, dissipation/viscosity:'
        WRITE(NCREPT,'(3(1PE12.4))')TENRGY,TINSCL,TDISIP
        WRITE(NCREPT,*)

C       RSC 09-JUN-2008 ENHANCEMENTS
C       COMPUTE VISCOSITIES
        VISCOS = ALAMDA*EXP(RLAMDA*LOG(TRIN))
        VISCOS = VISCOS*PRANTL
C       REFERENCE DENSITY
        RTIN = ZERO
        DO ISPEC = 1,NSPEC
          RTIN = RTIN + RGSPEC(ISPEC)*YRIN(ISPEC)
        ENDDO
        RTIN = RTIN*TRIN
        DVRIN = PRIN/RTIN
        VISKIN = VISCOS/DVRIN

        WRITE(NCREPT,*)'viscosities: dynamic, kinematic:'
        WRITE(NCREPT,'(2(1PE12.4))')VISCOS,VISKIN
        WRITE(NCREPT,*)
C       RSC 09-JUN-2008 ENHANCEMENTS
C       CONVERT TO DIMENSIONAL DISSIPATION
        TDISIP = TDISIP*VISKIN
        WRITE(NCREPT,*)'energy, length scale, dissipation:'
        WRITE(NCREPT,'(3(1PE12.4))')TENRGY,TINSCL,TDISIP
        WRITE(NCREPT,*)

C       RSC 09-JUN-2008 ENHANCEMENTS
C       COMPUTE TAYLOR LENGTH SCALE
        TAYLOR = TEN*VISKIN*TENRGY/TDISIP
        TAYLOR = SQRT(TAYLOR)
C       COMPUTE KOLMOGOROV LENGTH SCALE
        TKOLMO = VISKIN*VISKIN*VISKIN/TDISIP
        TKOLMO = SQRT(SQRT(TKOLMO))
        WRITE(NCREPT,*)'length scales: int, Taylor, Kolmo:'
        WRITE(NCREPT,'(3(1PE12.4))')TINSCL,TAYLOR,TKOLMO
        WRITE(NCREPT,*)

        WRITE(NCREPT,*)'domain/length scale ratios: int, Taylor, Kolmo:'
        WRITE(NCREPT,'(3(1PE12.4))')XGDLEN/TINSCL,XGDLEN/TAYLOR,
     +                              XGDLEN/TKOLMO
        WRITE(NCREPT,*)

C       RSC 09-JUN-2008 ENHANCEMENTS
        WRITE(NCREPT,*)'------------------------------------'
        WRITE(NCREPT,*)

      ENDIF

C     =========================================================================


      RETURN
      END
