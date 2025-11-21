      SUBROUTINE BCUTXL
 
C     *************************************************************************
C
C     BCUTXL
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     30-DEC-2003:  CREATED
C     04-JAN-2007:  RSC REVISE PARALLEL RECEIVES
C
C     DESCRIPTION
C     -----------
C     DNS CODE SENGA2
C     EVALUATES TIME-DEPENDENT BOUNDARY CONDITIONS FOR VELOCITY COMPONENTS
C     AND THEIR TIME DERIVATIVES
C
C     X-DIRECTION LEFT-HAND END
C
C     *************************************************************************


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
CKA   FIX INFLOW BUG, BTIME IS DEFINED IN COM_SENGA2.H
CKA      DOUBLE PRECISION BTIME
      DOUBLE PRECISION FORNOW,ARGMNT,ARGVAL,REALKX 
      DOUBLE PRECISION COSVAL,SINVAL,COSTHT,SINTHT
      DOUBLE PRECISION PCOUNT
      INTEGER IC,JC,KC
      INTEGER IIC,IIM,KX,KXBASE
      INTEGER ICPROC,NCOUNT,IRPROC,IRTAG
C     VM: SYNTHETIC DIGITAL FILTERING METHOD
      DOUBLE PRECISION VFAC,DVFDT,COFLOW

C     FY - FOR NON-REFLECTING INLOW (INFLOW OPTION 4)
      DOUBLE PRECISION DELTAGY,YCOORD
      INTEGER IGOFSTY,IY
      DOUBLE PRECISION LAMBDA
      PARAMETER(LAMBDA = 3.74D-3) ! SOUND WAVELENGTH
      DOUBLE PRECISION PULRAT
      PARAMETER(PULRAT = 0.1D0) ! RATIO OF PULSE WIDTH TO SOUND WAVELENGTH
      DOUBLE PRECISION PTLY
      PARAMETER(PTLY = PULRAT*LAMBDA) ! BASE WIDTH OF PULSE
      DOUBLE PRECISION WIDTHP
      PARAMETER(WIDTHP = 0.5D0) ! ADDITIONAL WIDTH PARAMETER - DEFAULT 0.5
      DOUBLE PRECISION SLOPE
      PARAMETER(SLOPE = 2.0D4) ! LARGER VALUE RESULTS IN SHARPER SLOPE
      DOUBLE PRECISION HYPTAN ! FOR HYPERBOLIC TANGENT PROFILE

C     BEGIN
C     =====

C     =========================================================================

CKA   THIS WAS MOVED TO BOUNDT & BOUNTT TO FIX INFLOW SCANNING LOCATION
C     RK TIME INCREMENT IS HELD IN RKTIM(IRKSTP)
CKA      BTIME = ETIME + RKTIM(IRKSTP)

C     =========================================================================

C     CONSTANT U-VELOCITY
C     PARAMETER I1=1, R1=U-VELOCITY
      IF(NXLPRM(1).EQ.1)THEN

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL

            STRUXL(JC,KC) = RXLPRM(1)
C            FORNOW = REAL(JC)*DELTAY/(HALF*YGDLEN)
C            STRUXL(JC,KC) = RXLPRM(1)*TANH(FORNOW)
            STRVXL(JC,KC) = ZERO
            STRWXL(JC,KC) = ZERO

            DUDTXL(JC,KC) = ZERO
            DVDTXL(JC,KC) = ZERO
            DWDTXL(JC,KC) = ZERO

C           WRITE ETIME, PRESSURE
C            IF (IRKSTP.EQ.NRKSTP) THEN
C              OPEN(UNIT=16,FILE="vars.dat",ACCESS='APPEND')
C              WRITE(16,'(3E20.9E3)')ETIME,(STRPXL(JC,1)-PRIN),STRVEL
C              CLOSE(16)
C            ENDIF

          ENDDO
        ENDDO

      ENDIF

C     =========================================================================

C     SINUSOIDAL U-VELOCITY
C     PARAMETER I1=2, R1=AMPLITUDE, R2=PERIOD

C     FY - EDITED FOR USE WITH NON-REFLECTING INFLOW

      IF(NXLPRM(1).EQ.2)THEN
C
C        FORNOW = TWO*PI/RXLPRM(2)
C        ARGMNT = FORNOW*BTIME
C
C        DO KC = KSTAL,KSTOL
C          DO JC = JSTAL,JSTOL
C
C            STRUXL(JC,KC) = RXLPRM(1)*SIN(ARGMNT)
C            STRVXL(JC,KC) = ZERO
C            STRWXL(JC,KC) = ZERO
C
C            DUDTXL(JC,KC) = FORNOW*RXLPRM(1)*COS(ARGMNT)
C            DVDTXL(JC,KC) = ZERO
C            DWDTXL(JC,KC) = ZERO
C
C          ENDDO
C        ENDDO
        DELTAGY = YGDLEN/REAL(NYGLBL-1)

        IGOFSTY = 0
        DO ICPROC = 0, IYPROC-1
          IGOFSTY = IGOFSTY + NPMAPY(ICPROC)
        ENDDO
      
          FORNOW = TWO*PI/RXLPRC(8)
          ARGMNT = FORNOW*(ETIME+RXLPRC(9)*RXLPRC(8)) 

          DO KC = KSTAL,KSTOL
            DO JC = JSTAL,JSTOL
C             GLOBAL Y COORD
              IY = IGOFSTY + JC 

C             NXLPRM(4)=0 - OPTION FOR STANDARD FLAT SINUSOIDAL INFLOW VELOCITY IN TIME ON WHOLE XL FACE 
              IF (NXLPRM(4).EQ.0) THEN

                STRUXL(JC,KC) = RXLPRM(1)+RXLPRC(7)*SIN(ARGMNT)
                STRVXL(JC,KC) = ZERO
                STRWXL(JC,KC) = ZERO

                DUDTXL(JC,KC) = FORNOW*RXLPRC(7)*COS(ARGMNT)
                DVDTXL(JC,KC) = ZERO
                DWDTXL(JC,KC) = ZERO

C               WRITE ETIME, PRESSURE
C                IF (IRKSTP.EQ.NRKSTP) THEN
C                  OPEN(UNIT=16,FILE="vars.dat",ACCESS='APPEND')
C                  WRITE(16,'(3E20.9E3)')ETIME,(STRPXL(JC,1)-PRIN),STRVEL
C                  CLOSE(16)
C                ENDIF

C             NXLPRM(4)=4 - OPTION FOR SINUSOIDAL VELOCITY IN TIME FOR PART OF XL FACE - ACTS AS A POINT SOURCE 
              ELSEIF (NXLPRM(4).EQ.4) THEN

C               DEFINE HYPERBOLIC TANGENT PROFILE
                HYPTAN=(0.5D0*TANH(SLOPE*(DBLE(IY-1)*DELTAGY-
     +          (YGDLEN/2.0D0-WIDTHP*PTLY)))+0.5D0)*
     +                (-0.5D0*TANH(SLOPE*(DBLE(IY-1)*DELTAGY-
     +          (YGDLEN/2.0D0+WIDTHP*PTLY)))+0.5D0)

C               SET VELOCITY ON FACE BASED ON PROFILE
                STRUXL(JC,KC) = RXLPRM(1)+HYPTAN*RXLPRC(7)*SIN(ARGMNT)
                STRVXL(JC,KC) = ZERO
                STRWXL(JC,KC) = ZERO
     
                DUDTXL(JC,KC) = HYPTAN*FORNOW*RXLPRC(7)*COS(ARGMNT)
                DVDTXL(JC,KC) = ZERO
                DWDTXL(JC,KC) = ZERO

C               WRITE ETIME, PRESSURE
C                IF (IRKSTP.EQ.NRKSTP) THEN
C                  OPEN(UNIT=16,FILE="vars.dat",ACCESS='APPEND')
C                  WRITE(16,'(3E20.9E3)')ETIME,(STRPXL(JC,1)-PRIN),STRVEL
C                  CLOSE(16)
C                ENDIF


              END IF

            ENDDO
          ENDDO

      ENDIF

C     =========================================================================

CVMC     TURBULENT VELOCITY FIELD
CVMC     PARAMETER I1=3
CVM      IF(NXLPRM(1).EQ.3)THEN
CVM
CVMC       INTERPOLATE STORED TURBULENT VELOCITY FIELD ONTO INLET PLANE
CVMC       DO THE INTERPOLATION BY DFT: LOCAL-PROCESSOR CONTRIBUTION
CVM
CVMC       -----------------------------------------------------------------------
CVM
CVMC       UPDATE THE SCANNING PLANE LOCATION
CVM        SLOCXL = ELOCXL - SVELXL*BTIME
CVM        IF(SLOCXL.LT.ZERO)SLOCXL = XGDLEN + SLOCXL
CVMCKA     FIX INFLOW
CVMCKA        IF(IRKSTP.EQ.NRKSTP)ELOCXL = SLOCXL
CVM        IF(FUPELC)ELOCXL = SLOCXL
CVM
CVMC       INITIALISE THE PHASE ANGLE TERMS
CVM        ARGMNT = TPOVXG*SLOCXL
CVM        COSTHT = COS(ARGMNT)
CVM        SINTHT = SIN(ARGMNT)
CVM        KXBASE = KMINXL
CVM
CVMC       ZERO THE LOCAL-PROCESSOR CONTRIBUTION TO THE DFT
CVM        DO KC = KSTAL,KSTOL
CVM          DO JC = JSTAL,JSTOL
CVM
CVM            STRUXL(JC,KC) = ZERO
CVM            STRVXL(JC,KC) = ZERO
CVM            STRWXL(JC,KC) = ZERO
CVM
CVM            DUDTXL(JC,KC) = ZERO
CVM            DVDTXL(JC,KC) = ZERO
CVM            DWDTXL(JC,KC) = ZERO
CVM
CVM          ENDDO
CVM        ENDDO
CVM
CVMC       -----------------------------------------------------------------------
CVM
CVMC       SPECIAL CASE OF LEADING IMAGINARY TERM
CVM        IF(FLLIXL)THEN
CVM
CVM          KX = KXBASE
CVM          REALKX = REAL(KX)
CVM          ARGVAL = ARGMNT*REALKX
CVM          COSVAL = COS(ARGVAL)
CVM          SINVAL = SIN(ARGVAL)
CVM          IIC = 1
CVM
CVM          DO KC = KSTAL,KSTOL
CVM            DO JC = JSTAL,JSTOL
CVM
CVM              STRUXL(JC,KC) = STRUXL(JC,KC)
CVM     +                      + UFXL(IIC,JC,KC)*SINVAL
CVM              STRVXL(JC,KC) = STRVXL(JC,KC)
CVM     +                      + VFXL(IIC,JC,KC)*SINVAL
CVM              STRWXL(JC,KC) = STRWXL(JC,KC)
CVM     +                      + WFXL(IIC,JC,KC)*SINVAL
CVM
CVM              DUDTXL(JC,KC) = DUDTXL(JC,KC)
CVM     +                      - REALKX*UFXL(IIC,JC,KC)*COSVAL
CVM              DVDTXL(JC,KC) = DVDTXL(JC,KC)
CVM     +                      - REALKX*VFXL(IIC,JC,KC)*COSVAL
CVM              DWDTXL(JC,KC) = DWDTXL(JC,KC)
CVM     +                      - REALKX*WFXL(IIC,JC,KC)*COSVAL
CVM
CVM            ENDDO
CVM          ENDDO
CVM
CVM          KXBASE = KXBASE + 1
CVM
CVM        ENDIF
CVM
CVMC       -----------------------------------------------------------------------
CVM
CVMC       STANDARD LOCAL CONTRIBUTION
CVM
CVMC       ZEROTH WAVENUMBER
CVM        IF(KXBASE.EQ.0)THEN
CVM
CVM          DO KC = KSTAL,KSTOL
CVM            DO JC = JSTAL,JSTOL
CVM
CVM              KX = KXBASE
CVM              REALKX = REAL(KX)
CVM              ARGVAL = ARGMNT*REALKX
CVM              COSVAL = COS(ARGVAL)
CVM              SINVAL = SIN(ARGVAL)
CVM              IIM = 1
CVM
CVM              STRUXL(JC,KC) = STRUXL(JC,KC)
CVM     +                      + HALF*UFXL(IIM,JC,KC)*COSVAL
CVM              STRVXL(JC,KC) = STRVXL(JC,KC)
CVM     +                      + HALF*VFXL(IIM,JC,KC)*COSVAL
CVM              STRWXL(JC,KC) = STRWXL(JC,KC)
CVM     +                      + HALF*WFXL(IIM,JC,KC)*COSVAL
CVM
CVM            ENDDO
CVM          ENDDO
CVM
CVM          KXBASE = KXBASE + 1
CVM
CVM        ENDIF
CVM
CVMC       ALL OTHER WAVENUMBERS
CVM        DO KC = KSTAL,KSTOL
CVM          DO JC = JSTAL,JSTOL
CVM
CVM            KX = KXBASE
CVM            REALKX = REAL(KX)
CVM            ARGVAL = ARGMNT*REALKX
CVM            COSVAL = COS(ARGVAL)
CVM            SINVAL = SIN(ARGVAL)
CVM
CVM            DO IC = ISTAXL,ISTOXL,2
CVM
CVM              IIM = IC
CVM              IIC = IC+1
CVM
CVM              STRUXL(JC,KC) = STRUXL(JC,KC)
CVM     +                      + UFXL(IIM,JC,KC)*COSVAL
CVM     +                      + UFXL(IIC,JC,KC)*SINVAL
CVM              STRVXL(JC,KC) = STRVXL(JC,KC)
CVM     +                      + VFXL(IIM,JC,KC)*COSVAL
CVM     +                      + VFXL(IIC,JC,KC)*SINVAL
CVM              STRWXL(JC,KC) = STRWXL(JC,KC)
CVM     +                      + WFXL(IIM,JC,KC)*COSVAL
CVM     +                      + WFXL(IIC,JC,KC)*SINVAL
CVM
CVM              DUDTXL(JC,KC) = DUDTXL(JC,KC)
CVM     +                      + REALKX*(UFXL(IIM,JC,KC)*SINVAL
CVM     +                              - UFXL(IIC,JC,KC)*COSVAL)
CVM              DVDTXL(JC,KC) = DVDTXL(JC,KC)
CVM     +                      + REALKX*(VFXL(IIM,JC,KC)*SINVAL
CVM     +                              - VFXL(IIC,JC,KC)*COSVAL)
CVM              DWDTXL(JC,KC) = DWDTXL(JC,KC)
CVM     +                      + REALKX*(WFXL(IIM,JC,KC)*SINVAL
CVM     +                              - WFXL(IIC,JC,KC)*COSVAL)
CVM
CVM              KX = KX + 1
CVM              REALKX = REAL(KX)
CVM              FORNOW = COSVAL
CVM              COSVAL = COSTHT*COSVAL - SINTHT*SINVAL
CVM              SINVAL = SINTHT*FORNOW + COSTHT*SINVAL
CVM
CVM            ENDDO
CVM
CVM          ENDDO
CVM        ENDDO
CVM
CVMC       -----------------------------------------------------------------------
CVM
CVMC       SPECIAL CASE OF TRAILING REAL TERM
CVM        IF(FLTRXL)THEN
CVM
CVM          KX = KXBASE + ISTOXL/2
CVM          REALKX = REAL(KX)
CVM          ARGVAL = ARGMNT*REALKX
CVM          COSVAL = COS(ARGVAL)
CVM          SINVAL = SIN(ARGVAL)
CVM          IIM = ISTOXL + 1
CVM
CVM          DO KC = KSTAL,KSTOL
CVM            DO JC = JSTAL,JSTOL
CVM
CVM              STRUXL(JC,KC) = STRUXL(JC,KC)
CVM     +                      + UFXL(IIM,JC,KC)*COSVAL
CVM              STRVXL(JC,KC) = STRVXL(JC,KC)
CVM     +                      + VFXL(IIM,JC,KC)*COSVAL
CVM              STRWXL(JC,KC) = STRWXL(JC,KC)
CVM     +                      + WFXL(IIM,JC,KC)*COSVAL
CVM
CVM              DUDTXL(JC,KC) = DUDTXL(JC,KC)
CVM     +                      + REALKX*UFXL(IIM,JC,KC)*SINVAL
CVM              DVDTXL(JC,KC) = DVDTXL(JC,KC)
CVM     +                      + REALKX*VFXL(IIM,JC,KC)*SINVAL
CVM              DWDTXL(JC,KC) = DWDTXL(JC,KC)
CVM     +                      + REALKX*WFXL(IIM,JC,KC)*SINVAL
CVM
CVM            ENDDO
CVM          ENDDO
CVM
CVM        ENDIF
CVM
CVMC       -----------------------------------------------------------------------
CVM
CVMC       PARALLEL TRANSFER
CVMC       RSC 04-JAN-2007 REVISE PARALLEL RECEIVES
CVM        IF(IXPROC.EQ.0)THEN
CVM
CVMC         LEFTMOST PROCESSOR IN X
CVMC         RECEIVE FROM ALL OTHER PROCESSORS IN X
CVM          DO ICPROC = 1,NXPRM1
CVM
CVM            IRPROC = NPROCX(ICPROC)
CVM            IRTAG = IRPROC*NPROC+IPROC
CVM            CALL P_RECV(PCOUNT,1,IRPROC,IRTAG)
CVM            CALL P_RECV(PARRAY,NPARAY,IRPROC,IRTAG)
CVM
CVM            NCOUNT = 0
CVM            DO KC = KSTAL,KSTOL
CVM              DO JC = JSTAL,JSTOL
CVM
CVM                NCOUNT = NCOUNT + 1
CVM                STRUXL(JC,KC) = STRUXL(JC,KC) + PARRAY(NCOUNT)
CVM                NCOUNT = NCOUNT + 1
CVM                STRVXL(JC,KC) = STRVXL(JC,KC) + PARRAY(NCOUNT)
CVM                NCOUNT = NCOUNT + 1
CVM                STRWXL(JC,KC) = STRWXL(JC,KC) + PARRAY(NCOUNT)
CVM                NCOUNT = NCOUNT + 1
CVM                DUDTXL(JC,KC) = DUDTXL(JC,KC) + PARRAY(NCOUNT)
CVM                NCOUNT = NCOUNT + 1
CVM                DVDTXL(JC,KC) = DVDTXL(JC,KC) + PARRAY(NCOUNT)
CVM                NCOUNT = NCOUNT + 1
CVM                DWDTXL(JC,KC) = DWDTXL(JC,KC) + PARRAY(NCOUNT)
CVM
CVM              ENDDO
CVM            ENDDO
CVM          
CVM          ENDDO
CVM
CVMC         SCALING OF DFT
CVM          DO KC = KSTAL,KSTOL
CVM            DO JC = JSTAL,JSTOL
CVM
CVMC             VELOCITIES
CVM              STRUXL(JC,KC) = STRUXL(JC,KC)*SCAUXL
CVM              STRVXL(JC,KC) = STRVXL(JC,KC)*SCAUXL
CVM              STRWXL(JC,KC) = STRWXL(JC,KC)*SCAUXL
CVM
CVMC             DERIVATIVES
CVM              DUDTXL(JC,KC) = DUDTXL(JC,KC)*SCDUXL
CVM              DVDTXL(JC,KC) = DVDTXL(JC,KC)*SCDUXL
CVM              DWDTXL(JC,KC) = DWDTXL(JC,KC)*SCDUXL
CVM
CVMC             ADD MEAN VELOCITY
CVM              STRUXL(JC,KC) = STRUXL(JC,KC) + BVELXL
CVM
CVMC             CONVERT SPATIAL TO TEMPORAL DERIVATIVES
CVM              DUDTXL(JC,KC) = DUDTXL(JC,KC)*SVELXL
CVM              DVDTXL(JC,KC) = DVDTXL(JC,KC)*SVELXL
CVM              DWDTXL(JC,KC) = DWDTXL(JC,KC)*SVELXL
CVM
CVM            ENDDO
CVM          ENDDO
CVM
CVM        ELSE 
CVM
CVMC         NOT THE LEFTMOST PROCESSOR IN X
CVMC         SEND TO LEFTMOST PROCESSOR IN X
CVM          NCOUNT = 0
CVM          DO KC = KSTAL,KSTOL
CVM            DO JC = JSTAL,JSTOL
CVM
CVM              NCOUNT = NCOUNT + 1
CVM              PARRAY(NCOUNT) = STRUXL(JC,KC)
CVM              NCOUNT = NCOUNT + 1
CVM              PARRAY(NCOUNT) = STRVXL(JC,KC)
CVM              NCOUNT = NCOUNT + 1
CVM              PARRAY(NCOUNT) = STRWXL(JC,KC)
CVM              NCOUNT = NCOUNT + 1
CVM              PARRAY(NCOUNT) = DUDTXL(JC,KC)
CVM              NCOUNT = NCOUNT + 1
CVM              PARRAY(NCOUNT) = DVDTXL(JC,KC)
CVM              NCOUNT = NCOUNT + 1
CVM              PARRAY(NCOUNT) = DWDTXL(JC,KC)
CVM
CVM            ENDDO
CVM          ENDDO
CVM
CVM          PCOUNT = REAL(NCOUNT)
CVM          IRPROC = NPROCX(0)
CVM          IRTAG = IPROC*NPROC+IRPROC
CVM          CALL P_SEND(PCOUNT,1,1,IRPROC,IRTAG)
CVM          CALL P_SEND(PARRAY,NPARAY,NCOUNT,IRPROC,IRTAG)
CVM
CVM        ENDIF
CVM
CVM      ENDIF

C     GENERATING TURBULENT FIELD USING SYNTHETIC DIGITAL FILTERING
C     METHOD
C     VM: NXLPRM(1)=4 IMPLIES THAT THE VELOCITY SYTHETIC DIGITAL FILTERING
C     IS ON
      IF(NXLPRM(1).EQ.4)THEN

c        VFAC=MIN(1.0D0,0.1+ETIME/0.0001)
        VFAC=ONE

        COFLOW = RXLPRM(4)

        IF (VFAC.LT.1.0D0) THEN
          DVFDT=10000.0D0
        ELSE
          DVFDT=0.0D0
        END IF

        DO KC = KSTAL,KSTOL
          DO JC = JSTAL,JSTOL
            STRUXL(JC,KC) = (RXLPRM(1)*USTEAD(JC,KC)+
     +                      UINF2(JC,KC))*VFAC+COFLOW
            STRVXL(JC,KC) = VINF2(JC,KC)*VFAC
            STRWXL(JC,KC) = WINF2(JC,KC)*VFAC


            DUDTXL(JC,KC) = (UINF2(JC,KC)-UINF1(JC,KC))/TSTEP
     +             *VFAC+(RXLPRM(1)*USTEAD(JC,KC)+UINF2(JC,KC))
     +             *DVFDT
            DVDTXL(JC,KC) = (VINF2(JC,KC)-VINF1(JC,KC))/TSTEP
     +             *VFAC+VINF2(JC,KC)*DVFDT
            DWDTXL(JC,KC) = (WINF2(JC,KC)-WINF1(JC,KC))/TSTEP
     +             *VFAC+WINF2(JC,KC)*DVFDT


          ENDDO
        ENDDO

      ENDIF

C     =========================================================================


      RETURN
      END
