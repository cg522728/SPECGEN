PROGRAM MMSENS
    USE     :: ISO_FORTRAN_ENV
    USE     :: xraylib
    USE     :: TYPES
    USE     :: XRLDATA
    USE     :: CFGDATA
    USE     :: ANODE
    USE     :: SECCOMP
    USE     :: MICROMATTER

    IMPLICIT NONE
    INTEGER :: CNT
    INTEGER :: CNT2
    INTEGER :: NELEMENT = 1
    CHARACTER(LEN=16)   :: ARG0
    REAL(WP)    :: ISAM = 0_QP
    REAL(WP)    :: ISAM1 = 0_QP
    REAL(WP)    :: ISAM2 = 0_QP
    REAL(WP)    :: ISAM3 = 0_QP
    REAL(DP)    :: EI = 0_DP
    REAL(WP), DIMENSION(:,:), ALLOCATABLE :: I_CHAR

    OPEN(UNIT=121,FILE='OUTPUT_MM.DAT', ACCESS='APPEND', STATUS='REPLACE')

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
    CALL LOAD_NIST()

    ALLOCATE(I_CHAR(CP_ST%NELEMENTS,SIZE(LINE)))
    I_CHAR = 0_QP
    DO CNT2 = 1, CP_ST%NELEMENTS
        DO CNT = 1, SIZE(LINE)
            WRITE (6,'(1H[, A16, 2H]>,I4,1H/,I4,A1)',ADVANCE='NO') 'MAIN', CNT, SIZE(LINE), CHAR(13)
            EI = LINE_ENERGY(CP_ST%ELEMENTS(CNT2), CNT)
            IF (EI.EQ.0 .OR. EI.LT.EMIN) THEN
                I_CHAR(CNT2,CNT) = 0_QP
            ELSE
                I_CHAR(CNT2,CNT) = ST_CHAR(CP_ST%ELEMENTS(CNT2), CNT)
            ENDIF
        END DO
    END DO
    DO CNT2 = 10, 92!CP_SAM%NELEMENTS
        NELEMENT = CNT2!CP_SAM%ELEMENTS(CNT2)
        DO CNT = 3, 3!1, SIZE(LINE)
            EI = LINE_ENERGY(NELEMENT, CNT)
            IF (EI.EQ.0) CYCLE
            IF (EI.LT.EMIN) CYCLE
            ISAM = MM_SENS(NELEMENT, CNT, I_CHAR)
            WRITE (6,'(I3,1H/,I3,A1,$)',ADVANCE='NO') CNT2, 92, CHAR(13)
        END DO
        WRITE (121,201) NELEMENT, ISAM
        ISAM = 0_WP
    END DO
201 FORMAT(I4, 2X, 3ES32.20E4)
202 FORMAT(I3, I3, ES32.20E3, 2X, 2ES32.20E4)
    CLOSE(121)
END PROGRAM MMSENS
