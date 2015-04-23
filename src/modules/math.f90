MODULE MATHCHG
    USE :: TYPES
    implicit none
CONTAINS
    FUNCTION INTEGRATE(FUNC, XMIN, XMAX, NSTEP, Z, N)
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        REAL(QP)                        :: INTEGRATE
        REAL(QP)                        :: FUNC
        EXTERNAL                        :: FUNC
        REAL(DP), INTENT(IN)            :: XMIN
        REAL(DP), INTENT(IN)            :: XMAX
        INTEGER, INTENT(IN)             :: NSTEP
        INTEGER, INTENT(IN), OPTIONAL   :: Z
        INTEGER, INTENT(IN), OPTIONAL   :: N
        REAL(QP)    :: FUNCVAL = 0_QP
        REAL(QP)    :: TMP = 0_QP
        REAL(QP)    :: FI = 0_QP
        REAL(QP)    :: FF = 0_QP
        REAL(DP)    :: XVAL = 0_DP
        REAL(DP)    :: XSTEP = 0_DP
        INTEGER     :: CNT
        LOGICAL     :: N_PRESENT = .FALSE.
        LOGICAL     :: Z_PRESENT = .FALSE.

        N_PRESENT = PRESENT(N)
        Z_PRESENT = PRESENT(Z)
        XSTEP = (XMAX-XMIN)*(DBLE(NSTEP)**(-1))
        FUNCVAL = 0_QP
        IF (Z_PRESENT) THEN
            IF (N_PRESENT) THEN
                FI = FUNC(XMIN, Z, N)/2
                FF = FUNC(XMAX, Z, N)/2
                DO CNT = 1, (NSTEP-1)
                !WRITE (6,'(1H[, A16, 2H]*,I4,1H/,I4,A1,$)',ADVANCE='NO') 'INTEGRATE', CNT, NSTEP-1, CHAR(13)
                    XVAL = XMIN + CNT*XSTEP
                    TMP = FUNC(XVAL, Z, N)
                    FUNCVAL = FUNCVAL + TMP
                END DO
                INTEGRATE = XSTEP*(FI + FUNCVAL + FF)
                RETURN
            ELSE
                FI = FUNC(XMIN, Z)/2
                FF = FUNC(XMAX, Z)/2
                DO CNT = 1, (NSTEP-1)
                 !WRITE (6,'(1H[, A16, 2H]*,I4,1H/,I4,A1,$)',ADVANCE='NO') 'INTEGRATE', CNT, NSTEP-1, CHAR(13)
                    XVAL = XMIN + CNT*XSTEP
                    TMP = FUNC(XVAL, Z)
                    FUNCVAL = FUNCVAL + TMP
                END DO
                INTEGRATE = XSTEP*(FI + FUNCVAL + FF)
                RETURN
            ENDIF
        ELSE
            IF (N_PRESENT) THEN
                FI = FUNC(XMIN, N)/2
                FF = FUNC(XMAX, N)/2
                DO CNT = 1, (NSTEP-1)
                !WRITE (6,'(1H[, A16, 2H]*,I4,1H/,I4,A1,$)',ADVANCE='NO') 'INTEGRATE', CNT, NSTEP-1, CHAR(13)
                    XVAL = XMIN + CNT*XSTEP
                    TMP = FUNC(XVAL, N)
                    FUNCVAL = FUNCVAL + TMP
                END DO
                INTEGRATE = XSTEP*(FI + FUNCVAL + FF)
                RETURN
            ELSE
                FI = FUNC(XMIN)/2
                FF = FUNC(XMAX)/2
                DO CNT = 1, (NSTEP-1)
                 !WRITE (6,'(1H[, A16, 2H]*,I4,1H/,I4,A1,$)',ADVANCE='NO') 'INTEGRATE', CNT, NSTEP-1, CHAR(13)
                    XVAL = XMIN + CNT*XSTEP
                    TMP = FUNC(XVAL)
                    FUNCVAL = FUNCVAL + TMP
                END DO
                INTEGRATE = XSTEP*(FI + FUNCVAL + FF)
                RETURN
            ENDIF
        ENDIF
        RETURN
    END FUNCTION INTEGRATE
END MODULE MATHCHG
