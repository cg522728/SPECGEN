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

        REAL(QP), DIMENSION(:), ALLOCATABLE :: ITMP
        REAL(QP), DIMENSION(:), ALLOCATABLE :: TMP

        INTEGER     :: CNT
        INTEGER     :: CNT2
        REAL(QP)    :: PI = 0_QP
        REAL(QP)    :: INTEN = 0_QP
        REAL(DP)    :: EI = 0_DP

        ALLOCATE(TMP(SIZE(LINE)))
        ALLOCATE(ITMP(CP_ST%NELEMENTS))

        TMP = 0_QP
        ITMP = 0_QP
        PI = 2.D0*DASIN(1.D0)

        DO CNT2 = 1, CP_ST%NELEMENTS
            DO CNT = 1, SIZE(LINE)
                WRITE (6,'(I4,1H/,I4,A1,$)',ADVANCE='NO') CNT, SIZE(LINE), CHAR(13)
                INTEN = I_CHAR(CNT2,CNT)
                EI = LineEnergy(CP_ST%ELEMENTS(CNT2), LINE(CNT))
                IF (EI.LE. 0) CYCLE
                    TMP(CNT) = INTEN&
                        *CS_FLUOR_CHG(NELEMENT, N, EI)
            END DO
            ITMP(CNT2) = SUM(TMP)*CP_ST%massFractions(CNT2)
            TMP = 0_QP
        END DO
        MM1 = SUM(ITMP)*((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))
        WRITE (6,'(1H[,A16,2H](,A3,1H),ES32.20E3)') 'MM1', 'OUT', MM1
        DEALLOCATE(TMP)
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
            EI = LineEnergy(Z_ANODE, LINE(CNT))
            IF (EI.LE. EdgeEnergy(NELEMENT, SHELL(N)).OR. EI.EQ. 0) CYCLE
                TMPR = AN_SCAT_RAYL(CNT)&
                        *CS_FLUOR_CHG(NELEMENT, N, EI)
                EC = ComptonEnergy(EI, DBLE(A_ST_AZIM_OUT))
                IF (EC.LT.EdgeEnergy(NELEMENT, SHELL(N)) .OR. EC.EQ.0) THEN
                    TMPC = 0_QP
                ELSE
                    TMPC = AN_SCAT_COMP(CNT)&
                            *CS_FLUOR_CHG(NELEMENT, N, EC)
                ENDIF
            TMP(CNT) = TMPR + TMPC
        END DO
        ITMP = SUM(TMP)
        MM2 = ITMP*((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))
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
        ETUBE = DBLE(VTUBE)
        E_ST_EDGE = EdgeEnergy(NELEMENT, SHELL(N))
        TMP = INTEGRATE(DERIV_MM3, E_ST_EDGE, ETUBE, INT(10), NELEMENT, N)
        MM3 = TMP*((SA_ST_IN*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))
        RETURN
    END FUNCTION
    FUNCTION DERIV_MM3(EI, Z, N)
        implicit none
        REAL(QP)    :: DERIV_MM3
        REAL(DP),   INTENT(IN)     :: EI
        INTEGER,    INTENT(IN)     :: Z
        INTEGER,    INTENT(IN)     :: N
        DERIV_MM3 = SECT_CONT(EI)&
                *CS_FLUOR_CHG(Z, N, EI)
    END FUNCTION DERIV_MM3
END MODULE
