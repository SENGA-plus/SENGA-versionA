      SUBROUTINE INFLOW
 
C     *************************************************************************
C
C     INFLOW
C     ======
C
C     AUTHOR
C     ------
C     M.KLEIN  --   UNIVERSITÄT DER BUNDESWEHR MÜNCHEN
C
C     CHANGE RECORD
C     -------------
C     27-MAR-2015:  CREATED
C     10-APR-2019:  CHANGED FOR SENGA2 (M.A.) 
C
C     DESCRIPTION
C     -----------

C     *************************************************************************
c      IMPLICIT NONE


C     GLOBAL DATA
C     ===========
C     -------------------------------------------------------------------------
      INCLUDE 'mpif.h'
      INCLUDE 'com_senga2.h'
C     -------------------------------------------------------------------------


C     LOCAL DATA
C     ==========
      REAL    RAN1        ! Random number Generator
      DOUBLE PRECISION LENX,LENY,LENZ,THETA,UMEAN,DELX,PHI,TSTP,URMS

      DOUBLE PRECISION NORM,AY(-NFY:NFY),AZ(-NFZ:NFZ)  ! filter coefficients  
      
      LOGICAL PERIY,PERIZ!,TRADIT
      PARAMETER (PERIY=.TRUE.,PERIZ=.FALSE.)!,TRADIT=.TRUE.)

      INTEGER YINFL,YINFR,ZINFL,ZINFR,nylcl,nzlcl      ! rectangular inflow area
      PARAMETER (YINFL=1,YINFR=NYGLBL,
     +           ZINFL=1,ZINFR=NZGLBL)

      DOUBLE PRECISION UFILT(YINFL:NYSIZE,ZINFL:NZSIZE),
     +                 VFILT(YINFL:NYSIZE,ZINFL:NZSIZE),
     +                 WFILT(YINFL:NYSIZE,ZINFL:NZSIZE),  ! filtered inflow data
     +                 YFILT(YINFL:NYSIZE,ZINFL:NZSIZE,NSPEC)

      DOUBLE PRECISION UFOLD(YINFL:NYSIZE,ZINFL-NFZ:NZSIZE+NFZ),
     +                 VFOLD(YINFL:NYSIZE,ZINFL-NFZ:NZSIZE+NFZ),
     +                 WFOLD(YINFL:NYSIZE,ZINFL-NFZ:NZSIZE+NFZ),
     +                 YFOLD(YINFL:NYSIZE,ZINFL-NFZ:NZSIZE+NFZ,NSPEC)

      DOUBLE PRECISION 
     + URAND(YINFL-NFY:YINFR+NFY,ZINFL-NFZ:ZINFR+NFZ),
     + VRAND(YINFL-NFY:YINFR+NFY,ZINFL-NFZ:ZINFR+NFZ),
     + WRAND(YINFL-NFY:YINFR+NFY,ZINFL-NFZ:ZINFR+NFZ),
     + YRAND(YINFL-NFY:YINFR+NFY,ZINFL-NFZ:ZINFR+NFZ,NSPEC)

      DOUBLE PRECISION UMEANGLBL,SUMUMEANGLBL

      INTEGER IC,JC,KC,JC2,KC2,YP,ZP
      INTEGER JFILTSTART,KFILTSTART
      INTEGER TSPC(2),NSPC,ISPEC

C     VM: CALCULATING USTEAD
      DOUBLE PRECISION RAD,DENOM,SUMDENOM
      DOUBLE PRECISION ZDIST,YDIST
      DOUBLE PRECISION YREF(NSPEC)
      INTEGER KG,JG
      INTEGER YOFFSET,ZOFFSET,YMIDPNT,ZMIDPNT
C     VM: DEBUG
      INTEGER SNAP,ICPROC
      CHARACTER*40 FNAME
      CHARACTER*4 PSNAP

C     LNX, LNY, LNZ: LENGTH SCALES IN TERMS OF DELX
C     THETA: INDUCED TIME SCALE
C     PHI: REQUIRED WEIGHTING PARAMETER
C     NFY,NFZ: FILTER SIZE: NF=2*LN
      DELX = XGDLEN/REAL(NXGLBL-1)
      UMEAN = RXLPRM(1)
      TSTP = TSTEP
      URMS = RXLPRM(2)
      YRMS = 0.050d0
      LENX = DBLE(LNX)
      LENY = DBLE(LNY)
      LENZ = DBLE(LNZ)
      THETA=LENX*DELX/UMEAN
      PHI=1.0-TSTP/THETA

      ZP=INT(IPROC/(NXPROC*NYPROC))
      YP=INT((IPROC-ZP*NYPROC*NXPROC)/NXPROC)

      JFILTSTART=YP*NYSIZE
      KFILTSTART=ZP*NZSIZE

      TSPC = (/1,2/)
      NSPC = 2
      DO ISPEC=1,NSPEC
        YREF(ISPEC)=0.0D0
      ENDDO
      YREF(1)=1.0D0
      YREF(2)=0.233D0

C     ------------------------------------------------------------------
      
CC    URMS=MAX(5.0D0,3.0D0+FLOAT(ITIME)/10000.D0*1.0D0)
C     CHECK RAMP UP OF VELOCITY
C      IF (NXLPRM(2).EQ.1)THEN
C         TO BE DONE
C      END IF
      PI = FOUR*ATAN(1.0D0)
C     CALCULATE FILTER COEFFICIENTS
      NORM=0.0
      DO JC=-NFY,NFY
         AY(JC)=EXP(-PI*JC*JC/(2*LNY*LNY))
         NORM=NORM+AY(JC)**2
      ENDDO   
      NORM=sqrt(NORM)
      DO JC=-NFY,NFY
         AY(JC)=AY(JC)/NORM
      ENDDO         

      NORM=0.0
      DO KC=-NFZ,NFZ
         AZ(KC)=exp(-Pi*KC*KC/(2*LNZ*LNZ))
         NORM=NORM+AZ(KC)**2
      ENDDO   
      NORM=sqrt(NORM)
      DO KC=-NFZ,NFZ
         AZ(KC)=AZ(KC)/NORM
      ENDDO         

      INTRAN = ITIME

C     INITIALIZE RANDOM ARRAYS
      DO JC=YINFL-NFY,YINFR+NFY
         DO KC=ZINFL-NFZ,ZINFR+NFZ
            URAND(JC,KC)=RAN1(INTRAN)
            VRAND(JC,KC)=RAN1(INTRAN)
            WRAND(JC,KC)=RAN1(INTRAN)
         ENDDO         
      ENDDO         
      IF(NXLPRM(2)==1)THEN
        DO JC=YINFL-NFY,YINFR+NFY
           DO KC=ZINFL-NFZ,ZINFR+NFZ
              YRAND(JC,KC,:)=0.0d0
              DO ISPC = 1,NSPC
                ISPEC = TSPC(ISPC)
                YRAND(JC,KC,ISPEC)=RAN1(INTRAN)
              ENDDO
           ENDDO         
        ENDDO         
      ENDIF

      IF (PERIY) THEN   ! overwrite in a periodic manner
        DO JC=-NFY+1,0
           DO KC=ZINFL-NFZ,ZINFR+NFZ
              URAND(JC,KC)=URAND(JC+YINFR,KC)
              VRAND(JC,KC)=VRAND(JC+YINFR,KC)
              WRAND(JC,KC)=WRAND(JC+YINFR,KC)
           ENDDO         
        ENDDO         
        DO JC=YINFR+1,YINFR+NFY
           DO KC=ZINFL-NFZ,ZINFR+NFZ
              URAND(JC,KC)=URAND(JC-YINFR,KC)
              VRAND(JC,KC)=VRAND(JC-YINFR,KC)
              WRAND(JC,KC)=WRAND(JC-YINFR,KC)
           ENDDO         
        ENDDO         
        IF(NXLPRM(2)==1)THEN
          DO JC=-NFY+1,0
             DO KC=ZINFL-NFZ,ZINFR+NFZ
                DO ISPC = 1,NSPC
                  ISPEC = TSPC(ISPC)
                  YRAND(JC,KC,ISPEC)=YRAND(JC+YINFR,KC,ISPEC)
                ENDDO
             ENDDO         
          ENDDO         
          DO JC=YINFR+1,YINFR+NFY
             DO KC=ZINFL-NFZ,ZINFR+NFZ
                DO ISPC = 1,NSPC
                  ISPEC = TSPC(ISPC)
                  YRAND(JC,KC,ISPEC)=YRAND(JC-YINFR,KC,ISPEC)
                ENDDO
             ENDDO         
          ENDDO         
        ENDIF
      ENDIF

      IF (PERIZ) THEN   ! overwrite in a periodic manner
        DO KC=-NFZ+1,0
           DO JC=YINFL-NFY,YINFR+NFY
              URAND(JC,KC)=URAND(JC,KC+ZINFR)
              VRAND(JC,KC)=VRAND(JC,KC+ZINFR)
              WRAND(JC,KC)=WRAND(JC,KC+ZINFR)
           ENDDO         
        ENDDO         
        DO KC=ZINFR+1,ZINFR+NFZ
           DO JC=YINFL-NFY,YINFR+NFY
              URAND(JC,KC)=URAND(JC,KC-ZINFR)
              VRAND(JC,KC)=VRAND(JC,KC-ZINFR)
              WRAND(JC,KC)=WRAND(JC,KC-ZINFR)
           ENDDO         
        ENDDO         
        IF(NXLPRM(2)==1)THEN
          DO KC=-NFZ+1,0
             DO JC=YINFL-NFY,YINFR+NFY
                DO ISPC = 1,NSPC
                  ISPEC = TSPC(ISPC)
                  YRAND(JC,KC,ISPEC)=YRAND(JC,KC+ZINFR,ISPEC)
                ENDDO
             ENDDO         
          ENDDO         
          DO KC=ZINFR+1,ZINFR+NFZ
             DO JC=YINFL-NFY,YINFR+NFY
                DO ISPC = 1,NSPC
                  ISPEC = TSPC(ISPC)
                  YRAND(JC,KC,ISPEC)=YRAND(JC,KC-ZINFR,ISPEC)
                ENDDO
             ENDDO         
          ENDDO         
        ENDIF
      ENDIF

      UFILT(:,:)=0.0
      VFILT(:,:)=0.0
      WFILT(:,:)=0.0
      YFILT(:,:,:)=0.0
      UFOLD(:,:)=0.0
      VFOLD(:,:)=0.0
      WFOLD(:,:)=0.0
      YFOLD(:,:,:)=0.0

 
C     What is TRADIT?
      DO JC=JSTAL,JSTOL
         DO KC=KSTAL-NFZ,KSTOL+NFZ
            DO JC2=-NFY,NFY

              UFOLD(JC,KC)=UFOLD(JC,KC)+
     +          URAND(JC+JFILTSTART+JC2,KC+KFILTSTART)*AY(JC2)
              VFOLD(JC,KC)=VFOLD(JC,KC)+
     +          VRAND(JC+JFILTSTART+JC2,KC+KFILTSTART)*AY(JC2)
              WFOLD(JC,KC)=WFOLD(JC,KC)+
     +          WRAND(JC+JFILTSTART+JC2,KC+KFILTSTART)*AY(JC2)
            ENDDO         
         ENDDO         
      ENDDO         
      IF(NXLPRM(2)==1)THEN
        DO JC=JSTAL,JSTOL
           DO KC=KSTAL-NFZ,KSTOL+NFZ
              DO JC2=-NFY,NFY
                DO ISPC = 1,NSPC
                  ISPEC = TSPC(ISPC)
                  YFOLD(JC,KC,ISPEC)=YFOLD(JC,KC,ISPEC)+
     +              YRAND(JC+JFILTSTART+JC2,KC+KFILTSTART,ISPEC)*AY(JC2)
                ENDDO
              ENDDO         
           ENDDO         
        ENDDO         
      ENDIF

      DO JC=JSTAL,JSTOL
         DO KC=KSTAL,KSTOL
            DO KC2=-NFZ,NFZ
               UFILT(JC,KC)=UFILT(JC,KC)+
     +           UFOLD(JC,KC+KC2)*AZ(KC2)
               VFILT(JC,KC)=VFILT(JC,KC)+
     +           VFOLD(JC,KC+KC2)*AZ(KC2)
               WFILT(JC,KC)=WFILT(JC,KC)+
     +           WFOLD(JC,KC+KC2)*AZ(KC2)
            ENDDO         
         ENDDO         
      ENDDO
      IF(NXLPRM(2)==1)THEN
        DO JC=JSTAL,JSTOL
           DO KC=KSTAL,KSTOL
              DO KC2=-NFZ,NFZ
                 DO ISPC = 1,NSPC
                   ISPEC = TSPC(ISPC)
                   YFILT(JC,KC,ISPEC)=YFILT(JC,KC,ISPEC)+
     +               YFOLD(JC,KC+KC2,ISPEC)*AZ(KC2)
                 ENDDO
              ENDDO         
           ENDDO         
        ENDDO
      ENDIF

C     Turbulence intensity
      UFILT=UFILT*SQRT(URMS**TWO*(ONE-PHI**TWO))
      VFILT=VFILT*SQRT(URMS**TWO*(ONE-PHI**TWO))
      WFILT=WFILT*SQRT(URMS**TWO*(ONE-PHI**TWO))
      IF(NXLPRM(2)==1)THEN
        DO JC=JSTAL,JSTOL
          DO KC=KSTAL,KSTOL
            DO ISPC = 1,NSPC
              ISPEC = TSPC(ISPC)
              YFILT(JC,KC,ISPEC)=YFILT(JC,KC,ISPEC)*
     +          SQRT((YRMS*YREF(ISPEC))**TWO*(ONE-PHI**TWO))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      UINF1=UINF2
      VINF1=VINF2
      WINF1=WINF2
      IF(NXLPRM(2)==1)THEN
        DO JC=JSTAL,JSTOL
          DO KC=KSTAL,KSTOL
            DO ISPC = 1,NSPC
              ISPEC = TSPC(ISPC)
              YINF1(JC,KC,ISPEC)=YINF2(JC,KC,ISPEC)
            ENDDO
          ENDDO
        ENDDO
      ENDIF


      DO JC=JSTAL,JSTOL
         DO KC=KSTAL,KSTOL
            UINF2(JC,KC)=UFILT(JC,KC)+PHI*UINF1(JC,KC)
            VINF2(JC,KC)=VFILT(JC,KC)+PHI*VINF1(JC,KC)
            WINF2(JC,KC)=WFILT(JC,KC)+PHI*WINF1(JC,KC)
         ENDDO
      ENDDO
      IF(NXLPRM(2)==1)THEN
        DO JC=JSTAL,JSTOL
           DO KC=KSTAL,KSTOL
              DO ISPC = 1,NSPC
                ISPEC = TSPC(ISPC)
                YINF2(JC,KC,ISPEC)=YFILT(JC,KC,ISPEC)+
     +            PHI*YINF1(JC,KC,ISPEC)
              ENDDO
           ENDDO
        ENDDO
      ENDIF

C     VM: CALCULATING USTEAD

      YOFFSET = 0
      DO ICPROC = 0, IYPROC-1
        YOFFSET = YOFFSET + NPMAPY(ICPROC)
      ENDDO
      ZOFFSET = 0
      YMIDPNT = NYGLBL/2
      ZMIDPNT = 1
      if (NYGLBL==1) then
        DELTAY = 0.0d0
      else
        DELTAY = YGDLEN/DBLE(NYGLBL-1)
      end if
      if (NZGLBL==1) then
        DELTAZ = 0.0d0
      else
        DELTAZ = ZGDLEN/DBLE(NZGLBL-1)
      end if

      DO KC = KSTAL,KSTOL
        DO JC = JSTAL,JSTOL
          KG = ZOFFSET+KC
          JG = YOFFSET+JC
          YDIST = DELTAY*DBLE((JG-YMIDPNT))
          ZDIST = DELTAZ*DBLE((KG-ZMIDPNT))
          RAD=SQRT(YDIST**2+ZDIST**2)
          USTEAD(JC,KC)=1.0D0!0.5D0*(1.0D0-TANH((RAD-RXLPRM(3)*YGDLEN)
C     +                  /(0.001*YGDLEN)))
        END DO
      END DO

      DO KC=KSTAL,KSTOL
        DO JC=JSTAL,JSTOL
          UINF2(JC,KC)=UINF2(JC,KC)*USTEAD(JC,KC)
          VINF2(JC,KC)=VINF2(JC,KC)*USTEAD(JC,KC)
          WINF2(JC,KC)=0.0D0!WINF2(JC,KC)*USTEAD(JC,KC)
        END DO
      END DO
      IF(NXLPRM(2)==1)THEN
        DO KC=KSTAL,KSTOL
          DO JC=JSTAL,JSTOL
            DO ISPC = 1,NSPC
              ISPEC = TSPC(ISPC)
              YINF2(JC,KC,ISPEC)=YINF2(JC,KC,ISPEC)*USTEAD(JC,KC)
            ENDDO
          END DO
        END DO
      ENDIF
         
C Calculate the average speed in the entire domain, which is then
C subtracted from u_new
      SUMUMEANGLBL=ZERO
      DENOM = ZERO
      UMEANGLBL=ZERO
      DO JC=JSTAL,JSTOL
         DO KC=KSTAL,KSTOL
            UMEANGLBL=UMEANGLBL+UINF2(JC,KC)
            DENOM = DENOM+USTEAD(JC,KC)
         ENDDO         
      ENDDO         
C      UMEANGLBL=UMEANGLBL/DENOM!FLOAT(NYSIZE*NZSIZE)
      CALL MPI_ALLREDUCE(UMEANGLBL,SUMUMEANGLBL,1,
     +         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(DENOM,SUMDENOM,1,
     +         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      UINF2=UINF2-SUMUMEANGLBL/(SUMDENOM)

      SUMUMEANGLBL=ZERO
      DENOM = ZERO
      UMEANGLBL=ZERO
      DO JC=JSTAL,JSTOL
         DO KC=KSTAL,KSTOL
            UMEANGLBL=UMEANGLBL+VINF2(JC,KC)
            DENOM = DENOM+USTEAD(JC,KC)
         ENDDO         
      ENDDO         
C      UMEANGLBL=UMEANGLBL/DENOM!FLOAT(NYSIZE*NZSIZE)
      CALL MPI_ALLREDUCE(UMEANGLBL,SUMUMEANGLBL,1,
     +         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(DENOM,SUMDENOM,1,
     +         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      VINF2=VINF2-SUMUMEANGLBL/(SUMDENOM)

      SUMUMEANGLBL=ZERO
      DENOM = ZERO
      UMEANGLBL=ZERO
      DO JC=JSTAL,JSTOL
         DO KC=KSTAL,KSTOL
            UMEANGLBL=UMEANGLBL+WINF2(JC,KC)
            DENOM = DENOM+USTEAD(JC,KC)
         ENDDO         
      ENDDO         
C      UMEANGLBL=UMEANGLBL/DENOM!FLOAT(NYSIZE*NZSIZE)
      CALL MPI_ALLREDUCE(UMEANGLBL,SUMUMEANGLBL,1,
     +         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      CALL MPI_ALLREDUCE(DENOM,SUMDENOM,1,
     +         MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
      WINF2=WINF2-SUMUMEANGLBL/(SUMDENOM)

C     -------------------------------------------------------------------------

      END

C     ==================================================================
      FUNCTION RAN1(IDUM)
      INTEGER IDUM,IA,IM,IQ,IR,NTAB,NDIV
      REAL RAN1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER J,K,IV(NTAB),IY
      SAVE IV,IY
      DATA IV /NTAB*0/, IY /0/
      IF (IDUM.LE.0.OR.IY.EQ.0) THEN
        IDUM=MAX(-IDUM,1)
        DO 11 J=NTAB+8,1,-1
          K=IDUM/IQ
          IDUM=IA*(IDUM-K*IQ)-IR*K
          IF (IDUM.LT.0) IDUM=IDUM+IM
          IF (J.LE.NTAB) IV(J)=IDUM
11      CONTINUE
        IY=IV(1)
      ENDIF
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF (IDUM.LT.0) IDUM=IDUM+IM
      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RAN1=MIN(AM*IY,RNMX)
      RAN1=(RAN1*2.0-1.0)/0.577
      RETURN
      END
C     ==================================================================

