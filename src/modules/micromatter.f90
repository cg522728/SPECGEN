MODULE MICROMATTER
    USE :: xraylib
    USE :: TYPES
    USE :: XRLDATA
    USE :: CFGDATA
    USE :: ANODE
    USE :: SECCOMP
    USE :: MATHCHG
    implicit none
CONTAINS
    FUNCTION MM_SENS(Z, N, I_ST_CHAR) RESULT(I)
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP), DIMENSION(:,:), INTENT(IN) :: I_ST_CHAR
        REAL(WP) :: I

        REAL(DP) :: E_SAM_CHAR
        REAL(WP) :: TMP

        E_SAM_CHAR = LINE_ENERGY(Z, N)

        TMP = I_X_ST_CHAR(Z, N, I_ST_CHAR)&
            +I_X_AN_SCAT_CHAR(Z, N)&
            +I_X_AN_SCAT_CONT(Z ,N)
        I  = TMP*DETEFF(E_SAM_CHAR)
        RETURN
    END FUNCTION MM_SENS

    FUNCTION I_X_ST_CHAR(Z, N, I_ST_CHAR) RESULT (I)
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP), DIMENSION(:,:), INTENT(IN) :: I_ST_CHAR
        REAL(WP) :: I

        REAL(DP) :: E_ST_CHAR
        REAL(WP) :: TMP
        REAL(WP) :: I_TMP

        INTEGER :: CNT
        INTEGER :: CNT2

        DO CNT2 = 1, CP_ST%NELEMENTS
            DO CNT = 1, SIZE(LINE)
                E_ST_CHAR = LINE_ENERGY(CP_ST%ELEMENTS(CNT2), CNT)
                IF (E_ST_CHAR.EQ.0) CYCLE
                TMP = I_ST_CHAR(CNT2, CNT)&
                        *CS_FLUOR_CHG(Z, N, E_ST_CHAR)
                I_TMP = I_TMP + TMP
            END DO
        END DO
        I = I_TMP
        I_TMP = 0_WP
        RETURN
    END FUNCTION I_X_ST_CHAR

    FUNCTION I_X_AN_SCAT_CHAR(Z, N) RESULT(I)
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: I

        REAL(DP) :: E_AN_CHAR
        REAL(DP) :: E_AN_COMPT
        REAL(WP) :: I_R
        REAL(WP) :: I_C
        REAL(WP) :: TMP
        REAL(WP) :: I_TMP

        INTEGER :: CNT

        DO CNT = 1, SIZE(LINE)
            E_AN_CHAR = LINE_ENERGY(Z_ANODE, CNT)
            IF (E_AN_CHAR.EQ.0) CYCLE
            E_AN_COMPT = COMPTON_ENERGY(E_AN_CHAR, A_ST_POL)

            I_R = ST_SCAT_R(CNT)*CS_FLUOR_CHG(Z, N, E_AN_CHAR)
            I_C = ST_SCAT_C(CNT)*CS_FLUOR_CHG(Z, N, E_AN_COMPT)
            TMP = I_R + I_C
            I_TMP = I_TMP + TMP
        END DO
        I = I_TMP
        RETURN
    END FUNCTION I_X_AN_SCAT_CHAR

    FUNCTION I_X_AN_SCAT_CONT(Z ,N) RESULT(I)
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: I

        REAL(DP) :: E_SAM_EDGE

        E_SAM_EDGE = EDGE_ENERGY(Z, N)

        I = INTEGRATE(DERIV_I_X_AN_SCAT_CONT, E_SAM_EDGE, VTUBE, INT(100), Z, N)
        RETURN
    END FUNCTION I_X_AN_SCAT_CONT

    FUNCTION DERIV_I_X_AN_SCAT_CONT(E, Z, N) RESULT(I)
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: I

        I = ST_SCAT_CONT(E)*CS_FLUOR_CHG(Z, N, E)
        RETURN
    END FUNCTION DERIV_I_X_AN_SCAT_CONT

    FUNCTION DETEFF(EI)
        REAL(WP)   :: DETEFF
        REAL(DP), INTENT(IN)    :: EI

        REAL(DP)    :: EC
        REAL(WP)    :: F_PE
        REAL(WP)    :: CON1
        REAL(WP)    :: CON2
        REAL(WP)    :: ETA
        REAL(WP)    :: K1
        REAL(WP)    :: K2

        CON1 = CS_PHOTO_CHG(Z_DET_BODY, EI)/MAC(Z_DET_BODY, EI)
        EC = LINE_ENERGY(Z_DET_BODY, 3)
        CON2 = MAC(Z_DET_BODY, EC)/MAC(Z_DET_BODY, EI)
        ETA = 1-(0.5_WP*FLUOR_YIELD(Z_DET_BODY, 3)&
                *CON1&
                *(1-CON2*LOG(1+(1/CON2))))

        K1 = EXP(-(CS_PHOTO_CHG(Z_DET_WINDOW, EI)*D_DET_WINDOW*1E-4&
            +CS_PHOTO_CHG(Z_DET_DL, EI)*D_DET_DL*1E-4&
            +CS_PHOTO_CHG(Z_DET_GAP, EI)*D_DET_GAP*1E-4))
        K2 = 1-EXP(-MAC(Z_DET_BODY, EI)*(D_DET_BODY*1E-4))
        DETEFF = ETA*K1*K2
        RETURN
    END FUNCTION
END MODULE
