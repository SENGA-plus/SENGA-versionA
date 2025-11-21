C     *************************************************************************
C
C     mpif.h
C     ======
C
C     AUTHOR
C     ------
C     R.S.CANT  --  CAMBRIDGE UNIVERSITY ENGINEERING DEPARTMENT
C
C     CHANGE RECORD
C     -------------
C     07-OCT-2002:  CREATED
C
C     DESCRIPTION
C     -----------
C     STUB HEADER FILE FOR MPI
C
C     *************************************************************************


C     MPI PARAMETERS
C     ==============
      INTEGER MPI_DOUBLE_PRECISION
      PARAMETER(MPI_DOUBLE_PRECISION = 1)

      INTEGER MPI_STATUS_SIZE
      PARAMETER(MPI_STATUS_SIZE = 2)

      INTEGER MPI_SUM, MPI_MAX, MPI_MIN
      PARAMETER(MPI_SUM = 1, MPI_MAX = 2, MPI_MIN = 3)


C     MPI GLOBAL DATA
C     ===============
      INTEGER MPI_COMM_WORLD
      COMMON/MPICOM/MPI_COMM_WORLD


C     END OF FILE
C     ===========
