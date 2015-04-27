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
    FUNCTION MM1(N, NELEMENT, I_CHAR)
        IMPLICIT NONE
        REAL(QP) :: MM1
        INTEGER, INTENT(IN) :: NELEMENT
        INTEGER, INTENT(IN) :: N
        REAL(QP), DIMENSION(:,:), INTENT(IN) :: I_CHAR

        REAL(QP) :: ITMP
        REAL(QP) :: TMP

        INTEGER     :: CNT
        INTEGER     :: CNT2
        REAL(QP)    :: PI
        REAL(QP)    :: INTEN
        REAL(DP)    :: EI
        REAL(WP)    :: CS


        TMP = 0_QP
        ITMP = 0_QP
        PI = 2.D0*DASIN(1.D0)

        DO CNT2 = 1, CP_ST%NELEMENTS
            DO CNT = 1, SIZE(LINE)
                WRITE (6,'(I4,1H/,I4,A1,$)',ADVANCE='NO') CNT, SIZE(LINE), CHAR(13)
                INTEN = I_CHAR(CNT2,CNT)
                EI = LINE_ENERGY(CP_ST%ELEMENTS(CNT2), CNT)
                IF (EI.LE. 0) CYCLE
                    CS = CS_FLUOR_CHG(NELEMENT, N, EI)
                    TMP = TMP + (INTEN*CS)
            END DO
            ITMP = ITMP + TMP*CP_ST%massFractions(CNT2)
            TMP = 0_WP
        END DO
        MM1 = ITMP*((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_INCID)))
        WRITE (6,'(1H[,A16,2H](,A3,1H),ES32.20E3)') 'MM1', 'OUT', MM1
        ITMP = 0_WP
        RETURN
    END FUNCTION
    FUNCTION MM2(N, NELEMENT)
        IMPLICIT NONE
        REAL(QP) :: MM2
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: NELEMENT

        INTEGER :: CNT
        REAL(QP) :: ITMP = 0_QP
        REAL(QP) :: TMPR = 0_QP
        REAL(QP)    :: TMPC = 0_QP
        REAL(QP), DIMENSION(:), ALLOCATABLE :: TMP
        REAL(DP) :: EI = 0_DP
        REAL(QP) :: PI = 0_QP
        REAL(DP)    :: EC = 0_DP

        ALLOCATE(TMP(SIZE(LINE)))
        TMP = 0._QP

        PI = 2.D0*DASIN(1.D0)

        DO CNT = 1, SIZE(LINE)
            WRITE (6,'(I4,1H/,I4,A1,$)',ADVANCE='NO') CNT, SIZE(LINE), CHAR(13)
            TMPR = 0_QP
            TMPC = 0_QP
            EI = LINE_ENERGY(Z_ANODE, CNT)
            IF (EI.LE. EDGE_ENERGY(NELEMENT, N).OR. EI.EQ. 0) CYCLE
                TMPR = ST_SCAT_R(CNT)&
                        *CS_FLUOR_CHG(NELEMENT, N, EI)
                EC = COMPTON_ENERGY(EI, A_ST_TAKE_OFF)
                IF (EC.LT.EDGE_ENERGY(NELEMENT, N) .OR. EC.EQ.0) THEN
                    TMPC = 0_QP
                ELSE
                    TMPC = ST_SCAT_C(CNT)&
                            *CS_FLUOR_CHG(NELEMENT, N, EC)
                ENDIF
            TMP(CNT) = TMPR + TMPC
        END DO
        ITMP = SUM(TMP)
        MM2 = ITMP*((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_INCID)))
        WRITE (6,'(1H[,A16,2H](,A3,1H),ES32.20E3)') 'MM2', 'OUT', MM2
        DEALLOCATE(TMP)
        RETURN
    END FUNCTION
    FUNCTION MM3(N, NELEMENT)
        IMPLICIT NONE
        REAL(QP)    :: MM3
        INTEGER,    INTENT(IN)  :: N
        INTEGER,    INTENT(IN)  :: NELEMENT
        REAL(DP)    :: E_ST_EDGE = 0_DP
        REAL(DP)    :: ETUBE = 0_DP
        REAL(QP)    :: PI = 0_QP
        REAL(QP)    :: TMP = 0_QP

        PI = 2.D0*DASIN(1.D0)
        ETUBE = VTUBE
        E_ST_EDGE = EDGE_ENERGY(NELEMENT, N)
        TMP = INTEGRATE(DERIV_MM3, E_ST_EDGE, ETUBE, INT(10), NELEMENT, N)
        MM3 = TMP*((SA_ST_IN*CONC)/(4*PI*SIN(A_ST_INCID)))
        WRITE (6,'(1H[,A16,2H](,A3,1H),ES32.20E3)') 'MM3', 'OUT', MM3
        RETURN
    END FUNCTION
    FUNCTION DERIV_MM3(EI, Z, N)
        implicit none
        REAL(QP)    :: DERIV_MM3
        REAL(DP),   INTENT(IN)     :: EI
        INTEGER,    INTENT(IN)     :: Z
        INTEGER,    INTENT(IN)     :: N
        DERIV_MM3 = ST_SCAT_CONT(EI)&
                *CS_FLUOR_CHG(Z, N, EI)
    END FUNCTION DERIV_MM3

    FUNCTION DETEFF(EI)
        REAL(WP)   :: DETEFF
        REAL(DP), INTENT(IN)    :: EI
        REAL(WP)    :: K1
        REAL(WP)    :: K2
        K1 = MAC(Z_DET_WINDOW, EI)*D_DET_WINDOW*1E-4&
            +MAC(Z_DET_DL, EI)*D_DET_DL*1E-4&
            +MAC(Z_DET_GAP, EI)*D_DET_GAP*1E-4
        K2 = MAC(Z_DET_BODY, EI)*D_DET_BODY*1E-4
        DETEFF = EXP(-K1)*(1-EXP(-K2))
        RETURN
    END FUNCTION
END MODULE
