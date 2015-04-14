PROGRAM MICROMATTER
    !#################################################################################
    !#SPECGEN CALCULATES THE FOLLOWING DATA:                                         #
    !#      -TUBE SPECTRUM                                                           #
    !#      -SECONDARY TARGET SPECTRUM                                               #
    !#      -INTENSITY OF CHARACTERISTIC LINES OF A SAMPLE                           #
    !#                                                                               #
    !#CHRISTOPHER GEERLINGS 07/04/2015                                               #
    !#################################################################################
    USE     :: ISO_FORTRAN_ENV
    USE     :: xraylib
    USE     :: XRLDATA
    USE     :: CFGDATA
    USE     :: ANODE
    USE     :: CONSTANTS
    USE     :: SECCOMP
    USE     :: MICROMATTER

    IMPLICIT NONE
    INTEGER :: CNT, N, CNT2, NELEMENT
    CHARACTER(LEN=16)   :: ARG0
    REAL(16) :: EI, I, K, EC, EA
    REAL(16)    :: KAPPA, PI
    REAL(16)    :: ISAM1, ISAM2
    REAL        :: START, FINISH

    PI = 2.D0*DASIN(1.D0)

    OPEN(UNIT=999, FILE='DEBUG.LOG', ACCESS='APPEND', STATUS='REPLACE')
    OPEN(UNIT=100, STATUS='SCRATCH')    !BUFFER

    OPEN(UNIT=101, STATUS='SCRATCH')    !ANODE CONTINUUM
    OPEN(UNIT=102, STATUS='SCRATCH')    !ANODE CHARACTERISTIC LINES
    OPEN(UNIT=103, STATUS='SCRATCH')    !TUBE SPECTRUM
    OPEN(UNIT=111, STATUS='SCRATCH')    !SECONDARY TARGET CONTINUUM
    OPEN(UNIT=112, STATUS='SCRATCH')    !SECONDARY TARGET CHARACTERISTIC LINES
    OPEN(UNIT=113, STATUS='SCRATCH')    !SECONDARY TARGET SPECTRUM
    OPEN(UNIT=131, STATUS='SCRATCH')    !SAMPLE CHARACTERISTIC LINES

    OPEN(UNIT=121,FILE='OUTPUT.DAT', ACCESS='APPEND', STATUS='REPLACE')

    !SET UP CONFIG-FILES
    CALL INITCFG()
    !READ SELECTOR FROM COMMAND LINE
    !   = 0 THEN LOAD CONFIG FROM FILES AND DISPLAY MENU
    !   = 1 THEN READ PARAMETERS FROM COMMAND LINE
    CALL GET_COMMAND_ARGUMENT(1, ARG0)
    ARG0 = ADJUSTL(ARG0)
    READ(ARG0,'(I1)') CNT
    IF (CNT.EQ. 0) THEN
        CALL LOADCFG()
    ELSEIF (CNT.EQ. 1) THEN
        CALL READCFG()
    ENDIF

    CALL DEBUGCFG()
    CALL DISPCFG()

    !CALCULATE A VALUE FOR THE PROPORTIONAL CONSTANT KAPPA
    !WRITE (6,*) 'CALCULATING KAPPA'
    !CALL CHECKKAPPA(32._16, KAPPA)
    KAPPA = 1E4

    !WRITING LIST OF ENERGY VALUES IN CONTINUUM
    WRITE (6,*) 'INITIALIZING ENERGY VALUES'
    DO CNT = 1, NSTEP
        WRITE (6,'(I4,1H/,I4,A1,$)',ADVANCE='NO') CNT, NSTEP, CHAR(13)
        WRITE (100,200) (EMIN + ESTEP*DBLE(CNT))
    END DO
    REWIND(100)

    DO CNT2 = 1, CP_SAM%NELEMENTS
        DO CNT = 3, 3!1, SIZE(LINE)
            CALL CPU_TIME(START)
            NELEMENT = CP_SAM%ELEMENTS(CNT2)
            EI = LineEnergy(NELEMENT, LINE(CNT))
            IF (EI.EQ. 0) CYCLE
            IF (EI.LE. EMIN) CYCLE
            !WRITE (6,'(I3,1H/,I3,A1,$)',ADVANCE='NO') CNT, SIZE(LINE), CHAR(13)
            ISAM1 = MM1(CNT, NELEMENT, KAPPA)
            ISAM1 = ISAM1 + MM2(CNT, NELEMENT, KAPPA)
            ISAM1 = ISAM1 + MM3(CNT, NELEMENT, KAPPA)
            ISAM1 = ISAM1*DETEFF(EI)
            ISAM1 = ISAM1/(ITUBE*250._16)

            WRITE (131,202) CP_SAM%ELEMENTS(1), CNT, LineEnergy(NELEMENT, LINE(CNT)), CONC, ISAM1
            CALL CPU_TIME(FINISH)
            !write (6,'(1H[,A16,2H](,A3,1H), I32, 2ES32.20E3)') 'MAIN','OUT', CNT2, LineEnergy(NELEMENT, LINE(CNT)), isam1
        END DO
  END DO

    !WRITING OUTPUT TO FILE
    WRITE (6,*) 'WRITING OUTPUT TO FILE'
    REWIND(131)
    DO
        READ (131,202,END=999) NELEMENT, CNT, EI, CONC, ISAM1
        WRITE (121,202) NELEMENT, CNT, EI, CONC, ISAM1
    END DO
200 FORMAT(ES32.20E3)
201 FORMAT(ES32.20E3, 2X, ES32.20E4)
202 FORMAT(I3, I3, ES32.20E3, 2X, 2ES32.20E4)
999 CLOSE(100)
    CLOSE(101)
    CLOSE(102)
    CLOSE(103)
    CLOSE(104)
    CLOSE(111)
    CLOSE(112)
    CLOSE(121)
    CLOSE(131)
    CLOSE(999)
END PROGRAM MICROMATTER

    SUBROUTINE CHECKKAPPA(DELTA, KAPPA)
        USE     :: xraylib
        USE     :: XRLDATA
        USE     :: CFGDATA
        USE     :: CONSTANTS
        !#################################################################################
        !#THIS FUNCTION CALCULATES A VALUE FOR THE PROPORTIONALITY CONSTANT KAPPA        #
        !#BY CALCULATING CALC_KR() FOR EACH LINE AND INCREASING KAPPA BY 2^DELTA EACH    #
        !#LOOP. THIS LOOP CONTINUES UNTIL THE CHARACTERISTIC LINE INTENSITY IS GREATER   #
        !#THAN THE CONTINUUM.                                                            #
        !#################################################################################
        USE     :: ISO_FORTRAN_ENV
        IMPLICIT NONE
        INTEGER(INT64) :: N, BASE
        INTEGER :: CNT
        REAL(16) :: EI, EC
        REAL(16) :: TAU
        REAL(16) :: R, PI
        REAL(16), INTENT(IN)  ::DELTA
        REAL(16), INTENT(OUT) :: KAPPA
        REAL(16)    :: CON, K, NR

        DATA BASE/2/

        PI = 2.D0*DASIN(1.D0)
        K = 1.D0
        KAPPA = K
        DO N = 1, SIZE(LINE)
        WRITE (6,'(I4,1H/,I4,A1,$)',ADVANCE='NO') N, SIZE(LINE), CHAR(13)
        NR = N
        EC = ANINT(LineEnergy(Z_ANODE, LINE(N))/ESTEP)*ESTEP
        TAU = (((VTUBE/EC)*LOG((VTUBE/EC)/((VTUBE/EC)-1)))-1)
1       R = CALC_KR(Z_ANODE, EC, INT(N), K)*TAU
        WRITE (6,'(ES32.20,A1,$)',ADVANCE='NO') R, CHAR(13)
        IF (R.EQ. 0) CYCLE
        IF (R.LT. 1) THEN
            CON = 2_16**DELTA
            K = K + CON
            GO TO 1
        ENDIF
        IF (K.GE. KAPPA) KAPPA = K
        END DO
        RETURN
    END SUBROUTINE CHECKKAPPA
