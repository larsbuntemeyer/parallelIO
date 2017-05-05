      SUBROUTINE DECOMP
      !
      USE MO_ORG
      USE MO_PARORG
      USE grid_
      !
      IMPLICIT NONE
      !
C**********************************************************************
C     PURPOSE   : COMPUTE THE VALUES FOR THE DOMAIN DECOMPOSITION
C     DATE      : 16.09.94
C
C     EXTERNALS :              -----
C
C     PARAMATERLIST INPUT:     -----
C
C     PARAMETERLIST OUTPUT:    -----
C
C
C     METHOD:
C
C     TREATMENT OF ERRORS:     -----
C
C**********************************************************************

C     STATT DER VERWENDETEN VARIANTEN VON SEND UND RECEIVE FUER DIE VERSCHIE-
C     DENEN DATENTYPEN KOENNTE MAN AUCH ALLES MIT DEM DATENTYP MPI_BYTE
C     BEHANDELN. DANN IST ABER DIE PORTABILITAET NICHT MEHR GARANTIERT. FALLS
C     DIES TROTZDEM ERWUENSCHT IST, MUSS DANN IN DEN ZEILEN MIT "COUNT=..."
C     DIE "1" DURCH DIE ANZAHL DER BYTES FUER EIN ELEMENT DES JEWEILIGEN
C     DATENTYPS ERSETZT WERDEN
C     LOCAL VARIABLES
      INTEGER ISUB,JSUB, IX1,IX2,IX2L,ICOMP,JCOMP,
     &                   JY1,JY2,JY2L,I,K,
     &                   IX,IY,I1D,
     &                   JR1,IR1,JX1,IY1,IALLOI,IALLOJ,
     &                   ISUBSIZE(MAXPRC,2)

C     COMPUTATION OF THE DECOMPOSITION
      MOIE = NX_GLOBAL
      MOJE = NY_GLOBAL
      MOKE = NLAYER  

C     BOUNDARY LINES OF THE WHOLE DOMAIN ARE NOT CONSIDERED
      ICOMP = MOIE - 2
      JCOMP = MOJE - 2

C     NUMBER OF GRID-POINTS A SUBDOMAIN GETS AT LEAST: ISUB, JSUB
      ISUB = ICOMP / NPROCX
      JSUB = JCOMP / NPROCY

C     DETERMINE, HOW MANY SUBDOMAINS WILL GET *ISUB* (IX1) AND HOW
C     MANY WILL GET *ISUB+1* (IX2) GRID-POINTS
      IX2 = ICOMP - NPROCX * ISUB
      IX1 = NPROCX - IX2

C     DETERMINE, HOW MANY SUBDOMAINS WILL GET *JSUB* (JY1) AND HOW
C     MANY WILL GET *JSUB+1* (JY2) GRID-POINTS
      JY2 = JCOMP - NPROCY * JSUB
      JY1 = NPROCY - JY2

C     DETERMINE THE DISTRIBUTION OF THE SUBDOMAINS WITH DIFFERENT
C     SIZES; THE ONES WITH MORE GRID-POINTS ARE PLACED TO THE
C     LEFT AND RIGHT (LOWER AND UPPER) BOUNDARY, THE ONES WITH
C     LESS GRID-POINTS ARE PLACED TO THE MIDDLE OF THE DOMAIN.
      IX2L = IX2 / 2
      JY2L = JY2 / 2

C     COMPUTATION OF ISUBCOMP, ISUBPOS AND ISUBSIZE
C     TO COMPUTE ISUBNEIGH YOU STILL NEED IPRCIDS
      DO IY=1,NPROCY
         DO IX=1,NPROCX

C           1D-NUMBERING
            I1D = (IY-1)*NPROCX + IX

C           COMPUTATION OF ISUBCOMP AND ISUBPOS
            IF ( (1 .LE. IX) .AND. (IX .LE. IX2L) ) THEN
               ISUBCOMP(I1D,1) = ISUB + 1
               ISUBPOS(I1D,1)  = (IX-1)*(ISUB+1) + 2
               ISUBPOS(I1D,3)  = IX * (ISUB+1) + 1
            ELSE IF ( (IX2L+1 .LE. IX) .AND. (IX .LE. IX2L+IX1) ) THEN
               ISUBCOMP(I1D,1) = ISUB
               ISUBPOS(I1D,1)  = (IX-1)*ISUB + IX2L + 2
               ISUBPOS(I1D,3)  = IX * ISUB + IX2L + 1
            ELSE IF ( (IX2L+IX1+1 .LE. IX).AND.(IX .LE. NPROCX) ) THEN
               ISUBCOMP(I1D,1) = ISUB + 1
               ISUBPOS(I1D,1)  = (IX-1)*(ISUB+1) - IX1 + 2
               ISUBPOS(I1D,3)  = IX * (ISUB+1) - IX1 + 1
            ENDIF

C           COMPUTATION OF ISUBCOMP AND ISUBPOS
            IF ( (1 .LE. IY) .AND. (IY .LE. JY2L) ) THEN
               ISUBCOMP(I1D,2) = JSUB + 1
               ISUBPOS(I1D,2)  = (IY-1)*(JSUB+1) + 2
               ISUBPOS(I1D,4)  = IY * (JSUB+1) + 1
            ELSE IF ( (JY2L+1 .LE. IY) .AND. (IY .LE. JY2L+JY1) ) THEN
               ISUBCOMP(I1D,2) = JSUB
               ISUBPOS(I1D,2)  = (IY-1)*JSUB + JY2L + 2
               ISUBPOS(I1D,4)  = IY * JSUB + JY2L + 1
            ELSE IF ( (JY2L+JY1+1 .LE. IY).AND.(IY .LE. NPROCY) ) THEN
               ISUBCOMP(I1D,2) = JSUB + 1
               ISUBPOS(I1D,2)  = (IY-1)*(JSUB+1) - JY1 + 2
               ISUBPOS(I1D,4)  = IY * (JSUB+1) - JY1 + 1
            ENDIF

C           SIZE OF THE WHOLE DOMAIN FOR ISUBPOS
            ISUBPOS(I1D,5) = 1
            ISUBPOS(I1D,6) = 1
            ISUBPOS(I1D,7) =MOIE
            ISUBPOS(I1D,8) =MOJE

C           COMPUTATION OF ISUBSIZE
            ISUBSIZE(I1D,1) = ISUBCOMP(I1D,1) + 4
            ISUBSIZE(I1D,2) = ISUBCOMP(I1D,2) + 4

C           DECOMPOSITION FOR THE FFT
            JX1  = ISUBCOMP(I1D,2) / NPROCX
            JR1  = ISUBCOMP(I1D,2) - NPROCX * JX1
            IF (IX .LE. JR1) THEN
               NFFTJLO(I1D) = ISUBPOS(I1D,2)     + (IX-1) * (JX1+1)
               NFFTJUP(I1D) = ISUBPOS(I1D,2) - 1 +  IX    * (JX1+1)
            ELSE
               NFFTJLO(I1D) = ISUBPOS(I1D,2)     + (IX-1-JR1) *  JX1
     &                                           +       JR1  * (JX1+1)
               NFFTJUP(I1D) = ISUBPOS(I1D,2) - 1 + (IX  -JR1) *  JX1
     &                                           +       JR1  * (JX1+1)
            ENDIF

C           DECOMPOSITION FOR THE GAUSSELIMINATION
            IY1  = ISUBCOMP(I1D,1) / NPROCY
            IR1  = ISUBCOMP(I1D,1) - NPROCY * IY1
            IF (IY .LE. IR1) THEN
               NGAUSSILO(I1D) = ISUBPOS(I1D,1)   + (IY-1) * (IY1+1)
               NGAUSSIUP(I1D) = ISUBPOS(I1D,1)-1 +  IY    * (IY1+1)
            ELSE
               NGAUSSILO(I1D) = ISUBPOS(I1D,1)   + (IY-1-IR1) *  IY1
     &                                           +       IR1  * (IY1+1)
               NGAUSSIUP(I1D) = ISUBPOS(I1D,1)-1 + (IY  -IR1) *  IY1
     &                                           +       IR1  * (IY1+1)
            ENDIF

         ENDDO
      ENDDO

C     AB HIER NUR NOCH EIGENE GROESSEN

      DO IY=1,NPROCY
         DO IX=1,NPROCX
C         (* 1D NUMBERING *)
            I1D = (IY-1) * NPROCX + IX
            IF (MYID+1 .EQ. I1D) THEN
               MYPOS(1)=IX
               MYPOS(2)=IY
               MYI = ISUBCOMP(I1D,1)
               MYJ = ISUBCOMP(I1D,2)
               DO I=1,8
                  MYPOSGRD(I)=ISUBPOS(I1D,I)
               ENDDO
               EXIT
            ENDIF
         ENDDO
      ENDDO

C     MAXIMALE LOKALE DIMENSIONEN (MIT RAND UND UEBERLAPP)
      IALLOGI=0
      IALLOGJ=0
      DO K=1,NPROC
         IALLOI=ISUBCOMP(K,1)+4
         IALLOJ=ISUBCOMP(K,2)+4
         IF(IALLOI.GT.IALLOGI)IALLOGI=IALLOI
         IF(IALLOJ.GT.IALLOGJ)IALLOGJ=IALLOJ
      ENDDO

C*********************************************************************

C     BESTIMMUNG DER EINDIMENSIONALEN ZEILENAUFTEILUNG AUF DIE
C     PROZESSOREN FUER DIE FFT UND DER EINDIMENSIONALEN SPALTEN-
C     AUFTEILUNG FUER DIE GAUSSELIMINATION

C     FFT:   AUFTEILUNG DER ZEILEN
C     DIE ZEILEN WERDEN INNERHALB EINER PROZESSORENZEILE
C     VON J_LO = MYPOSGRD(2) BIS J_UP = MYPOSGRD(4)
C     AN DIE PROZESSOREN VON LINKS NACH RECHTS VERTEILT. DER LINKE
C     PROZESSOR IN DER PROZESSORENZEILE (MYPOS(1) = 1) BEKOMMT DIE
C     UNTEREN ZEILEN, DER RECHTE PROZESSOR (MYPOS(1) = NPROCX)
C     BEKOMMT DIE OBEREN ZEILEN. KOENNEN DIE ZEILEN NICHT
C     GLEICHMAESSIG AUFGETEILT WERDEN, ERHALTEN DIE ERSTEN PROZESSOREN
C     MEHR ZEILEN.

      JX1 = MYJ / NPROCX
      JR1 = MYJ - NPROCX * JX1
      IF (MYPOS(1) .LE. JR1) THEN
         MYFFT_JLO = MYPOSGRD(2)     + (MYPOS(1)-1) * (JX1+1)
         MYFFT_JUP = MYPOSGRD(2) - 1 +  MYPOS(1)    * (JX1+1)
      ELSE
         MYFFT_JLO = MYPOSGRD(2)     + (MYPOS(1)-1-JR1) *  JX1
     &                               + JR1              * (JX1+1)
         MYFFT_JUP = MYPOSGRD(2) - 1 + (MYPOS(1) - JR1) *  JX1
     &                               + JR1              * (JX1+1)
      ENDIF

C**** WIRD NICHT VERWENDET********
      MYFFT_ILO = MYPOSGRD(5) + 1
      MYFFT_IUP = MYPOSGRD(7) - 1
C********************************

C*********************************************************************

C     GAUSSELIMINATION:   AUFTEILUNG DER SPALTEN
C     DIE SPALTEN WERDEN INNERHALB EINER PROZESSORENSPALTE
C     VON I_LO = MYPOSGRD(1) BIS I_UP = MYPOSGRD(3)
C     AN DIE PROZESSOREN VON UNTEN NACH OBEN VERTEILT. DER UNTERE
C     PROZESSOR IN DER PROZESSORENSPALTE (MYPOS(2) = 1) BEKOMMT DIE
C     LINKEN SPALTEN, DER OBERSTE PROZESSOR (MYPOS(2) = NPROCY)
C     BEKOMMT DIE RECHTEN SPALTEN. KOENNEN DIE SPALTEN NICHT
C     GLEICHMAESSIG AUFGETEILT WERDEN, ERHALTEN DIE ERSTEN PROZESSOREN
C     MEHR SPALTEN

      IY1  = MYI / NPROCY
      IR1 = MYI - NPROCY * IY1
      IF (MYPOS(2) .LE. IR1) THEN
         MYGAUSS_ILO = MYPOSGRD(1)     + (MYPOS(2)-1) * (IY1+1)
         MYGAUSS_IUP = MYPOSGRD(1) - 1 +  MYPOS(2)    * (IY1+1)
      ELSE
         MYGAUSS_ILO = MYPOSGRD(1)     + (MYPOS(2)-1-IR1) *  IY1
     &                               + IR1              * (IY1+1)
         MYGAUSS_IUP = MYPOSGRD(1) - 1 + (MYPOS(2) - IR1) *  IY1
     &                               + IR1              * (IY1+1)
      ENDIF
      MYGAUSS_JLO = MYPOSGRD(6) + 1
      MYGAUSS_JUP = MYPOSGRD(8) - 1


C     COMPUTE THE NEIGHBORS OF THE SUBDOMAINS
      DO IY=1,NPROCY
         DO IX=1,NPROCX

C         1D NUMBERING
            I1D = (IY-1)*NPROCX + IX

            IF (IX .EQ. 1) THEN
C     NO NEIGHBOR TO THE WEST (LEFT)
               ISUBNEIGH(I1D,1) = -1
            ELSE
               ISUBNEIGH(I1D,1) = I1D-2
            ENDIF

            IF (IY .EQ. NPROCY) THEN
C     NO NEIGHBOR TO THE NORTH (UP)
               ISUBNEIGH(I1D,2) = -1
            ELSE
               ISUBNEIGH(I1D,2) = I1D+NPROCX-1
            ENDIF

            IF (IX .EQ. NPROCX) THEN
C     NO NEIGHBOR TO THE EAST (RIGHT)
               ISUBNEIGH(I1D,3) = -1
            ELSE
               ISUBNEIGH(I1D,3) = I1D+1-1
            ENDIF

            IF (IY .EQ. 1) THEN
C     NO NEIGHBOR TO THE SOUTH (DOWN)
               ISUBNEIGH(I1D,4) = -1
            ELSE
               ISUBNEIGH(I1D,4) = I1D-NPROCX-1
            ENDIF

         ENDDO
      ENDDO
      DO I=1,4
         NEIGHBOR(I)=ISUBNEIGH(MYID+1,I)
      ENDDO
      IE = ISUBSIZE(MYID+1,1)
      JE = ISUBSIZE(MYID+1,2)
      IEJE=IE*JE
      KE=MOKE
      IEKE=IE*KE
      IEJEKE=IEJE*KE
      KE1=KE+1

      CALL EXRANDASSSETUP

      P_REAL = MPI_DOUBLE_PRECISION

CKS
C 9911 FORMAT ('Verwendete Prozessoren:  ')
C 9912 FORMAT ('Maximal verwendbare Prozessoren:  ',I5)
C 9913 FORMAT ('Anzahl NPROC = ',I5,';   Aufteilung:    NPROCX = ',
C     +         I5,'     NPROCY = ',I5)
C
C 9901 FORMAT ('Knoten:  ',I5)
C 9902 FORMAT ('Groesse des Teilgebietes:   MYie = ',I5,
C     +           ';     MYje = ',I5)
C 9903 FORMAT ('Anzahl der Gitterpunkte, fuer die der neue Zeitschritt ',
C     +         'berechnet wird:')
C 9904 FORMAT (T29,'MYI  = ',I5,';     MYJ  = ',I5)
C 9905 FORMAT ('Position des Teilgebietes im Gesamtgebiet:')
C 9906 FORMAT ('       i(links unten) = ',I5,'    i(rechts oben) = ',I5)
C 9907 FORMAT ('       j(links unten) = ',I5,'    j(rechts oben) = ',I5)
CKS
      RETURN
      !
      END SUBROUTINE DECOMP
