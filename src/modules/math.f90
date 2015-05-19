MODULE MATHCHG
    USE :: TYPES
    implicit none
CONTAINS
    FUNCTION INTEGRATE(FUNC, XMIN, XMAX, NSTEP, Z, N)
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        REAL(WP)                        :: INTEGRATE
        EXTERNAL                        :: FUNC
        REAL(WP)                        :: FUNC
        REAL(DP), INTENT(IN)            :: XMIN
        REAL(DP), INTENT(IN)            :: XMAX
        INTEGER, INTENT(IN)             :: NSTEP
        INTEGER, INTENT(IN), OPTIONAL   :: Z
        INTEGER, INTENT(IN), OPTIONAL   :: N
        REAL(WP)    :: FUNCVAL
        REAL(WP)    :: TMP
        REAL(WP)    :: FI
        REAL(WP)    :: FF
        REAL(DP)    :: XVAL
        REAL(DP)    :: XSTEP
        INTEGER     :: CNT
        LOGICAL     :: N_PRESENT = .FALSE.
        LOGICAL     :: Z_PRESENT = .FALSE.

        N_PRESENT = PRESENT(N)
        Z_PRESENT = PRESENT(Z)
        XSTEP = (XMAX-XMIN)*(DBLE(NSTEP)**(-1))
        FUNCVAL = 0_QP
        TMP = 0_QP
        IF (Z_PRESENT .AND. N_PRESENT) THEN
            FI = FUNC(XMIN, Z, N)/2
            FF = FUNC(XMAX, Z, N)/2
            DO CNT = 1, (NSTEP-1)
                TMP = 0_QP
                XVAL = XMIN + CNT*XSTEP
                TMP = FUNC(XVAL, Z, N)
                FUNCVAL = FUNCVAL + TMP
            END DO
            INTEGRATE = XSTEP*(FI + FUNCVAL + FF)
            FUNCVAL = 0_QP
            RETURN
        ENDIF
        IF (Z_PRESENT .AND. .NOT.N_PRESENT) THEN
            FI = FUNC(XMIN, Z)/2
            FF = FUNC(XMAX, Z)/2
            DO CNT = 1, (NSTEP-1)
                TMP = 0_QP
                XVAL = XMIN + CNT*XSTEP
                TMP = FUNC(XVAL, Z)
                FUNCVAL = FUNCVAL + TMP
            END DO
            INTEGRATE = XSTEP*(FI + FUNCVAL + FF)
            FUNCVAL = 0_QP
            RETURN
        ENDIF
        IF (.NOT.Z_PRESENT .AND. .NOT.N_PRESENT) THEN
            FI = FUNC(XMIN)/2
            FF = FUNC(XMAX)/2
            DO CNT = 1, (NSTEP-1)
                TMP = 0_QP
                XVAL = XMIN + CNT*XSTEP
                TMP = FUNC(XVAL)
                FUNCVAL = FUNCVAL + TMP
            END DO
            INTEGRATE = XSTEP*(FI + FUNCVAL + FF)
            FUNCVAL = 0_QP
            RETURN
        ENDIF
        RETURN
    END FUNCTION INTEGRATE
END MODULE MATHCHG
