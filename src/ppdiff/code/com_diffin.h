C     DIFFIN-------------------------------------------------------------------

C     PARAMETERS
C     ==========
      INTEGER NSPCMX
      PARAMETER(NSPCMX = 100)

      INTEGER NDCFMX,NVCFMX,NCCFMX
      PARAMETER(NDCFMX = 4, NVCFMX = 4, NCCFMX = 4)

C     MOLECULAR TRANSPORT DATA
C     ========================
      DOUBLE PRECISION DIFFCO(NDCFMX,NSPCMX,NSPCMX)
      DOUBLE PRECISION TDRCCO(NDCFMX,NSPCMX,NSPCMX)
      DOUBLE PRECISION VISCCO(NVCFMX,NSPCMX)
      DOUBLE PRECISION CONDCO(NCCFMX,NSPCMX)
      DOUBLE PRECISION PREFGB,TREFGB

      INTEGER LENSYM(NSPCMX)
      INTEGER NSPEC
      INTEGER NCODIF,NCOTDR,NCOVIS,NCOCON
 
      CHARACTER*10 SPCSYM(NSPCMX)

      COMMON/CHEMIN/DIFFCO,TDRCCO,VISCCO,CONDCO,
     +              PREFGB,TREFGB,
     +              LENSYM,
     +              NSPEC,
     +              NCODIF,NCOTDR,NCOVIS,NCOCON,
     +              SPCSYM

C     DIFFIN-------------------------------------------------------------------
