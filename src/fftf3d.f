      SUBROUTINE FFTF3D(CARRRE,CARRIM,NXPHYS,NYPHYS,NZPHYS,
     +                                NX,NY,NZ,IFORW)

C     *************************************************************************
C
C     FFTF3D
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     29-DEC-1998:  CREATED
C
C     DESCRIPTION
C     -----------
C     CARRIES OUT AN FFT IN 3D
C
C     REFERENCES
C     ----------
C     1) NUMERICAL RECIPES pp451-453
C
C     *************************************************************************


C     PARAMETERS
C     ==========
      INTEGER NFTMAX
      PARAMETER(NFTMAX=1024)


C     ARGUMENTS
C     =========
      INTEGER NXPHYS,NYPHYS,NZPHYS,NX,NY,NZ,IFORW
      DOUBLE PRECISION CARRRE(NXPHYS,NYPHYS,NZPHYS)
      DOUBLE PRECISION CARRIM(NXPHYS,NYPHYS,NZPHYS)


C     LOCAL DATA
C     ==========
      DOUBLE PRECISION TARRAY(NFTMAX)
      INTEGER IX,JX,KX


C     BEGIN
C     =====

C     FT ON I-INDEX
C     -------------
      CALL FFTGIN(NX,IFORW)

      DO KX = 1,NZ
        DO JX = 1,NY

          DO IX = 1,NX
            TARRAY(2*IX-1) = CARRRE(IX,JX,KX)
            TARRAY(2*IX)   = CARRIM(IX,JX,KX)
          ENDDO

          CALL FFTGEN(TARRAY,NX)

          DO IX = 1,NX
            CARRRE(IX,JX,KX) = TARRAY(2*IX-1)
            CARRIM(IX,JX,KX) = TARRAY(2*IX)
          ENDDO

        ENDDO
      ENDDO


C     FT ON J-INDEX
C     -------------
      CALL FFTGIN(NY,IFORW)

      DO KX = 1,NZ
        DO IX = 1,NX

          DO JX = 1,NY
            TARRAY(2*JX-1) = CARRRE(IX,JX,KX)
            TARRAY(2*JX)   = CARRIM(IX,JX,KX)
          ENDDO

          CALL FFTGEN(TARRAY,NY)

          DO JX = 1,NY
            CARRRE(IX,JX,KX) = TARRAY(2*JX-1)
            CARRIM(IX,JX,KX) = TARRAY(2*JX)
          ENDDO

        ENDDO
      ENDDO


C     FT ON K-INDEX
C     -------------
      CALL FFTGIN(NZ,IFORW)

      DO JX = 1,NY
        DO IX = 1,NX

          DO KX = 1,NZ
            TARRAY(2*KX-1) = CARRRE(IX,JX,KX)
            TARRAY(2*KX)   = CARRIM(IX,JX,KX)
          ENDDO

          CALL FFTGEN(TARRAY,NZ)

          DO KX = 1,NZ
            CARRRE(IX,JX,KX) = TARRAY(2*KX-1)
            CARRIM(IX,JX,KX) = TARRAY(2*KX)
          ENDDO

        ENDDO
      ENDDO


      RETURN
      END
