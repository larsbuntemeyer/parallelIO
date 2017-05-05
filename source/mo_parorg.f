!
      MODULE MO_PARORG
!
      use mpi
      IMPLICIT NONE
!
!**********************************************************************
!
!     MODULE:    PARORG.F
!     CONTAINS ORGANIZATIONAL DATA OF THE PARALLEL (NODE-) PROGRAM
!
!**********************************************************************
!
!     DATE      : 02.02.94
!     AUTHOR    : ULRICH SCHAETTLER
!
!     11.11.15 - converted to fortran module (Lars Buntemeyer)
!                
!
!**********************************************************************
!      INCLUDE "mpif.h"

      INTEGER, PARAMETER :: MAXPRC=4096
!
      INTEGER :: IE, JE, KE, IEJE, IEKE, IEJEKE, KE1,
     &  MOIE,MOJE,MOKE,MOIEJE,MOIEKE,MOIEJEKE,MOKE1,
     &  MYID, MYPID, MYPOS(2), MYPOSGRD(8), NEIGHBOR(4),
     &  COMM,IERROR, SOURCE,TYPE,STATUS(MPI_STATUS_SIZE),COUNT,
     &  REQUEST(MAXPRC*100),REQUESTNR,DEST,
     &  REQUESTR(MAXPRC*100),REQUESTNRR,TAGCOUNT,TAGTABLE(MAXPRC),
     &  MYI, MYJ, NPROCX, NPROCY,IPRCIDS(MAXPRC),
     &  MYFFT_ILO, MYFFT_IUP, MYFFT_JLO, MYFFT_JUP,
     &  MYGAUSS_ILO, MYGAUSS_IUP, MYGAUSS_JLO, MYGAUSS_JUP,
     &  MPARTTREE(2,4),
     &  NFFTJLO(MAXPRC),NFFTJUP(MAXPRC),NGAUSSILO(MAXPRC),
     &  NGAUSSIUP(MAXPRC),ISUBPOS(MAXPRC,8),ISUBCOMP(MAXPRC,2),
     &  ISUBNEIGH(MAXPRC,4),IALLOGI,IALLOGJ,
     &  NPPMAX_SEND1,NPPMAX_RECV1,NPPINDEX1(MAXPRC),
     &  NPPMAX_SEND2,NPPMAX_RECV2,NPPINDEX2(MAXPRC),
     &  NPPMAX_SEND3,NPPMAX_RECV3,NPPINDEX3(MAXPRC),
     &  NPPMAX_SEND4,NPPMAX_RECV4,NPPINDEX4(MAXPRC),
     &  NPPMAX_SEND5,NPPMAX_RECV5,NPPINDEX5(MAXPRC),
     &  NPPMAX_SEND6,NPPMAX_RECV6,NPPINDEX6(MAXPRC)
      INTEGER :: NPPMAX_SEND7,NPPMAX_RECV7,NPPINDEX7(MAXPRC),
     &  NPPMAX_SEND8,NPPMAX_RECV8,NPPINDEX8(MAXPRC),
     &  IAMIN1,IEMAX1,JAMIN1,JEMAX1,
     &  IAMIN2,IEMAX2,JAMIN2,JEMAX2,
     &  IAMIN3,IEMAX3,JAMIN3,JEMAX3,
     &  IAMIN4,IEMAX4,JAMIN4,JEMAX4,
     &  IAMIN5,IEMAX5,JAMIN5,JEMAX5,
     &  IAMIN6,IEMAX6,JAMIN6,JEMAX6,
     &  IAMIN7,IEMAX7,JAMIN7,JEMAX7,
     &  IAMIN8,IEMAX8,JAMIN8,JEMAX8,
     &  IW5,IS5,IS6,IO6, IO7,IN7,IN8,IW8,P_REAL

       INTEGER, PARAMETER :: P_IO = 0
!
!**********************************************************************
!
!     DESCRIPTION OF THE VARIABLES:
!
!     MOIE        :   NUMBER OF GRIDPOINTS IN EAST-WEST-DIRECTION
!                     OF THE WHOLE MODEL
!     MOJE        :   NUMBER OF GRIDPOINTS IN NORTH-SOUTH-DIRECTION
!                     OF THE WHOLE MODEL
!     IE          :   NUMBER OF GRID-POINTS IN EAST-WEST-DIRECTION
!                     FOR ONE SUBDOMAIN
!     JE          :   NUMBER OF GRID-POINTS IN NORTH-SOUTH-DIRECTION
!                     FOR ONE SUBDOMAIN
!     KE          :   NUMBER OF GRID-POINTS IN VERTICAL DIRECTION     
!     IEJE        :   IE * JE
!     IEKE        :   IE * KE
!     IEJEKE      :   IE * JE * KE
!     KE1         :   KE + 1
!     NZLP        :   NUMBER OF GRID-POINTS BELONGING TO LAND-POINTS
!                     FOR ONE SUBDOMAIN
!
!     MYID        :   IDENTIFICATION OF THE PROCESS USING 1D NUMBERING
!     MYPID       :   PROCESS-IDENTIFICATION 
!     MYPOS       :   OWN POSITION IN THE PROCESSOR-TORUS IN X- AND
!                     Y-DIRECTION
!     MYPOSGRD    :   POSITION OF THE SUBDOMAIN IN THE ENTIRE GRID;
!                     CONTAINS THE COORDINATES OF THE LOWER LEFT AND
!                     THE UPPER RIGHT GRID-POINT IN THE WHOLE DOMAIN
!                     IN THE ORDER I_LL, J_LL, I_UR, J_UR.
!                     (BELONGING TO ISUBPOS IN *HPARORG*)
!     NEIGHBOR    :   CONTAINS THE PROCESS IDENTIFICATIONS OF THE 
!                     NEIGHBORS IN THE ORDER: WEST, NORTH, EAST, SOUTH
!                     (-1, IF NO NEIGHBORS ARE PRESENT)
!
!     MYI         :   NUMBER OF GRID-POINTS IN EAST-WEST-DIRECTION, FOR
!                     WHICH THE NEW TIME-LEVEL IS COMPUTED
!     MYJ         :   NUMBER OF GRID-POINTS IN NORTH-SOUTH-DIRECTION, FOR
!                     WHICH THE NEW TIME-LEVEL IS COMPUTED
!
!     NPROC       :   NUMBER OF PROCESSORS USED
!     NPROCX      :   NUMBER OF PROCESSORS IN EAST-WEST DIRECTION
!     NPROCY      :   NUMBER OF PROCESSORS IN NORTH-SOUTH DIRECTION
!
!     MYFFT_ILO   :   LOWER I-INDEX FOR THE COLUMNS BELONGING TO THIS 
!                     PROCESSOR FOR THE FFT
!     MYFFT_IUP   :   UPPER I-INDEX FOR THE COLUMNS BELONGING TO THIS 
!                     PROCESSOR FOR THE FFT
!     MYFFT_JLO   :   LOWER J-INDEX FOR THE ROWS BELONGING TO THIS 
!                     PROCESSOR FOR THE FFT
!     MYFFT_JUP   :   UPPER J-INDEX FOR THE ROWS BELONGING TO THIS 
!                     PROCESSOR FOR THE FFT
!     MYGAUSS_ILO :   LOWER I-INDEX FOR THE COLUMNS BELONGING TO THIS 
!                     PROCESSOR FOR THE GAUSS-ELIMINATION
!     MYGAUSS_IUP :   UPPER I-INDEX FOR THE COLUMNS BELONGING TO THIS 
!                     PROCESSOR FOR THE GAUSS-ELIMINATION
!     MYGAUSS_JLO :   LOWER J-INDEX FOR THE ROWS BELONGING TO THIS 
!                     PROCESSOR FOR THE GAUSS-ELIMINATION
!     MYGAUSS_JUP :   UPPER J-INDEX FOR THE ROWS BELONGING TO THIS 
!                     PROCESSOR FOR THE GAUSS-ELIMINATION
!
!**********************************************************************
      CONTAINS
      !
      !
      !
      SUBROUTINE P_ABORT
      IMPLICIT NONE
      !
      INTEGER :: P_ERROR
      CALL MPI_ABORT (MPI_COMM_WORLD, 1, P_ERROR)
      IF (p_error /= MPI_SUCCESS) THEN
         WRITE (*,'(a)') ' MPI_ABORT failed.'
         WRITE (*,'(a,i4)') ' Error =  ', P_ERROR
         STOP
      END IF
      WRITE (*,'(a)') 'mo_parorg: p_abort ...'
      !
      END SUBROUTINE P_ABORT
      !
      END MODULE MO_PARORG
