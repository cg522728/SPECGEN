PROGRAM STSPEC
    USE     :: ISO_FORTRAN_ENV
    USE     :: xraylib
    USE     :: TYPES
    USE     :: XRLDATA
    USE     :: CFGDATA
    USE     :: ANODE
    USE     :: SECCOMP

    IMPLICIT NONE
    INTEGER :: CNT, CNT2, N
    CHARACTER(LEN=16)   :: ARG0
    CHARACTER(LEN=16)    :: LABEL
    CHARACTER(LEN=3)    :: NLABEL
    REAL(DP) :: EI, EC, EA
    REAL(QP) :: I, ITMP
    REAL(QP)    :: PI
    REAL(QP)    :: TMP, I_CHAR
    REAL(QP), DIMENSION(:,:), ALLOCATABLE   :: I_ST_CHAR
    REAL(QP), DIMENSION(:), ALLOCATABLE   :: I_ST_CONT

    PI = 2.D0*DASIN(1.D0)

    OPEN(UNIT=999, FILE='DEBUG.LOG', ACCESS='APPEND', STATUS='REPLACE')
    OPEN(UNIT=100, STATUS='SCRATCH')    !BUFFER

    OPEN(UNIT=111, STATUS='SCRATCH')    !SECONDARY TARGET CONTINUUM
    OPEN(UNIT=112, STATUS='SCRATCH')    !SECONDARY TARGET CHARACTERISTIC LINES
    OPEN(UNIT=113, STATUS='SCRATCH')    !SECONDARY TARGET SPECTRUM

    OPEN(UNIT=121,FILE='OUTPUT_ST.DAT', ACCESS='APPEND', STATUS='REPLACE')

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

    CALL DISPCFG()
    CALL LOAD_NIST()

    WRITE (6,*) 'INITIALIZING ENERGY VALUES'
    DO CNT = 1, NSTEP
        WRITE (6,'(I4,1H/,I4,A1,$)',ADVANCE='NO') CNT, NSTEP, CHAR(13)
        WRITE (100,200) (EMIN + ESTEP*DBLE(CNT))
    END DO
    REWIND(100)

!    !CALCULATING INTENSITY OF SCATTERED CONTINUUM OF X-RAY TUBE
!    WRITE (6,*) 'CALCULATING INTENSITY OF SCATTERED CONTINUUM OF X-RAY TUBE'
!    REWIND(100)
!    DO CNT = 1,NSTEP
!        WRITE (6,'(I5,1H/,I5,A1,$)',ADVANCE='NO') CNT, NSTEP, CHAR(13)
!        READ (100,200) EI
!        WRITE (111,201) EI, SECTCONT(EI)
!    END DO
!
!    !CALCULATING INTENSITIES OF CHARACTERISTIC LINES
!    WRITE (6,*) 'CALCULATING INTENSITIES OF CHARACTERISTIC LINES'
!        DO N = 1, SIZE(LINE)
!            WRITE (6,'(I10,1H/,I10,3H-->,I3,1H/,I3,A1,$)',ADVANCE='NO'), CNT, NSTEP, N, SIZE(LINE), CHAR(13)
!            EA = LineEnergy(Z_ANODE, LINE(N))
!            IF (EA.EQ.0) CYCLE
!            IF (ISNAN(EA)) CYCLE
!            WRITE (112,201) EA, SECTRAYL(N)
!            EC = ComptonEnergy(DBLE(EA), DBLE((PI/2)-A_ST_AZIM_OUT))
!            WRITE (112,201) EC, SECTCOMP(N)
!            !$OMP PARALLEL DO SHARED(ITMP) PRIVATE(CNT2, EC)
!            DO CNT2 = 1, CP_ST%NELEMENTS
!                EC = LineEnergy(CP_ST%ELEMENTS(CNT2), LINE(N))
!                IF (EC.EQ.0) CYCLE
!                WRITE (112,201) EC, SECCOMPCHAR(N, CP_ST%ELEMENTS(CNT2))
!            END DO
!            !$OMP END PARALLEL DO
!        END DO
!    REWIND(112)
!    REWIND(111)
!    REWIND(100)

    !CALCULATING SECONDARY TARGET SPECTRUM USING ENERGY VALUES
    !CALCULATED EARLIER. IF THERE IS CHARACTERISTIC LINE AT THAT
    !POSITION, ADD SPACING OF ESTEP/10 BEFORE AND AFTER
    !CHARACTERISTIC LINE WITH THE INTENSITY OF THE SCATTERED
    !CONTINUUM. THEN PLOT A POINT WITH THE INTENSITY OF THE
    !CHARACTERISTIC LINE I NBETWEEN.
    !IF THERE IS NO CHARACTERISTIC LINE AT THAT POSITION,
    !PLOT INTENSITY OF SCATTERED CONTINUUM.
    !THIS IS REPEATED FOR THE SCATTERED CHARACTERISTIC LINES
    !OF THE ANODE (BOTH RAYLEIGH AND COMPTON SCATTERED).
    WRITE (6,*) 'CALCULATING SECONDARY TARGET SPECTRUM'
    REWIND(100)

    ALLOCATE(I_ST_CHAR(SIZE(LINE),CP_ST%NELEMENTS))
    ALLOCATE(I_ST_CONT(NSTEP))

    DO CNT = 1, NSTEP
        EI = 0_DP
        READ (100,200) EI
        I_ST_CONT(CNT) = ST_SCAT_CONT(EI)
!        IF (ISNAN(I_ST_CONT(CNT))) THEN
!            WRITE (6,*) 'ERROR:', EI, I_ST_CONT(CNT)
!        ENDIF
        WRITE (6,'(I10,1H/,I10,A1,$)',ADVANCE='NO'), CNT, NSTEP, CHAR(13)
    END DO
    REWIND(100)
    DO CNT= 1, SIZE(LINE)
        DO CNT2 = 1, CP_ST%NELEMENTS
            I_ST_CHAR(CNT, CNT2) = ST_CHAR(CP_ST%ELEMENTS(CNT2), CNT)
!            IF (ISNAN(I_ST_CHAR(CNT, CNT2))) THEN
!                WRITE (6,*) 'ERROR:', EI, I_ST_CHAR(CNT, CNT2)
!            ENDIF
        END DO
        WRITE (6,'(I10,1H/,I10,A1,$)',ADVANCE='NO'), CNT, SIZE(LINE), CHAR(13)
    END DO

    DO CNT = 1, NSTEP
        READ (100,200) EI
        TMP = 0_QP
        I_CHAR = 0_QP
        ITMP = 0!I_ST_CONT(CNT)
        LABEL = ''
        DO N = 1, SIZE(LINE)
            WRITE (6,'(I10,1H/,I10,3H-->,I3,1H/,I3,A1,$)',ADVANCE='NO'), CNT, NSTEP, N, SIZE(LINE), CHAR(13)
            EA = LineEnergy(Z_ANODE, LINE(N))
            IF (EA.EQ.0) CYCLE
            IF (ISNAN(EA)) CYCLE
            IF (EA.GE.EI .AND. EA.LT.(EI+ESTEP)) THEN
                TMP = ST_SCAT_R(N)
                ITMP = ITMP + TMP
                TMP = 0_QP
                LABEL = TRIM(LABEL)//'1'
            ENDIF
            EC = ComptonEnergy(DBLE(EA), DBLE((PI/2)-A_ST_TAKE_OFF))
            IF (EC.GE.EI .AND. EC.LT.(EI+ESTEP)) THEN
                TMP = ST_SCAT_C(N)
                ITMP = ITMP + TMP
                TMP = 0_QP
                LABEL = TRIM(LABEL)//'C'
            ENDIF
            DO CNT2 = 1, CP_ST%NELEMENTS
                EC = LineEnergy(CP_ST%ELEMENTS(CNT2), LINE(N))
                IF (EC.EQ.0) CYCLE
                IF (EC.LT.EI) CYCLE
                IF (EC.GE.EI .AND. EC.LT.(EI+ESTEP)) THEN
                    I_CHAR = I_ST_CHAR(N, CNT2)
                    ITMP = ITMP + I_CHAR
                    WRITE (NLABEL,'(I3)') N
                    NLABEL = TRIM(NLABEL)
                    !LABEL = TRIM(LABEL)//TRIM(NLABEL)
                ENDIF
            END DO
        END DO
        WRITE (104,201) EI, ITMP, LABEL
        IF (ISNAN(ITMP)) THEN
            WRITE (6,*) 'ERROR:', EI, ITMP
        ENDIF
    END DO
    WRITE (6,*) 'WRITING OUTPUT TO FILE'
    REWIND(104)
    DO
        READ (104,201,END=999) EI, I, LABEL
        WRITE (121,201) EI, I, LABEL
    END DO
    REWIND(104)

200 FORMAT(ES32.20E3)
201 FORMAT(ES32.20E3, 2X, ES32.20E4, 2X, A16)
999 CLOSE(100)
    CLOSE(104)
    CLOSE(111)
    CLOSE(112)
    CLOSE(121)
    CLOSE(999)
END PROGRAM STSPEC
