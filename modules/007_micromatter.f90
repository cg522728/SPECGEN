MODULE MICROMATTER
    USE :: xraylib
    USE :: XRLDATA
    USE :: CFGDATA
    USE :: ANODE
    USE :: SECCOMP
CONTAINS
    FUNCTION MM1(N, NELEMENT, KAPPA)
        IMPLICIT NONE
        REAL(16) :: MM1
        INTEGER :: CNT, CNT2, NELEMENT
        INTEGER :: N
        REAL(16) :: ITMP, TMP, EI
        REAL(16) :: PI, KAPPA, I_CHAR
        DATA ITMP/0/
        PI = 2.D0*DASIN(1.D0)
        DO CNT2 = 1, CP_ST%NELEMENTS
            I_CHAR = SECCOMPCHAR(N, CP_ST%ELEMENTS(CNT2), KAPPA)
            DO CNT = 1, SIZE(LINE)
                EI = LineEnergy(CP_ST%ELEMENTS(CNT2), LINE(CNT))
                IF (EI.LE. 0) CYCLE
                IF (EI.LE. LineEnergy(NELEMENT, LINE(N))) CYCLE
                TMP = I_CHAR&
                    *CS_Total_CP(STR_SAMPLE, LineEnergy(CP_ST%ELEMENTS(CNT2), LINE(CNT)))
1               ITMP = ITMP + TMP
            END DO
                ITMP = ITMP&
            *((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))&
            *FluorYield(NELEMENT, SHELL(N))&
            *RadRate(NELEMENT, LINE(N))
        END DO
        MM1 = ITMP
    END FUNCTION
    FUNCTION MM2(N, NELEMENT, KAPPA)
        IMPLICIT NONE
        REAL(16) :: MM2
        INTEGER :: CNT, CNT2, NELEMENT
        INTEGER :: LINES, SHELLS, N
        REAL(16) :: ITMP, TMP
        REAL(16) :: PI, KAPPA, EI, EC
        DATA ITMP/0/
        PI = 2.D0*DASIN(1.D0)
            DO CNT = 1, SIZE(LINE)
                EI = LineEnergy(Z_ANODE, LINE(CNT))
                IF (EI.LE. 0) CYCLE
                IF (EI.LE. LineEnergy(NELEMENT, LINE(N))) CYCLE
            TMP = SECTRAYL(KAPPA, N)&
                *CS_Total_CP(STR_SAMPLE, DBLE(EI))
            EC = ComptonEnergy(DBLE(EI), DBLE(A_ST_AZIM_OUT))
            TMP = TMP + SECTCOMP(KAPPA, N)&
                *CS_Total_CP(STR_SAMPLE, DBLE(EC))
            ITMP = ITMP + TMP
        END DO
        ITMP = ITMP&
            *((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))&
            *FluorYield(NELEMENT, SHELL(N))&
            *RadRate(NELEMENT, LINE(N))
        MM2 = ITMP
    END FUNCTION
    FUNCTION MM3(N, NELEMENT, KAPPA)
        IMPLICIT NONE
        REAL(16) :: MM3
        INTEGER :: CNT, CNT2, Z, INTSTEP
        INTEGER :: LINES, SHELLS, N, ECNT, NELEMENT
        REAL(16) :: ITMP, TMP
        REAL(16) :: PI, EI, KAPPA
        DATA ITMP/0/
        PI = 2.D0*DASIN(1.D0)
        INTSTEP = INT((VTUBE-EdgeEnergy(NELEMENT, SHELL(N)))/ESTEP)
        DO ECNT = 0, INTSTEP
            EI = EdgeEnergy(NELEMENT, SHELL(N)) + (DBLE(ECNT)*ESTEP)
            TMP = SECTCONT(EI)&
                *CS_Total_CP(STR_SAMPLE, DBLE(EI))
            ITMP = ITMP + TMP
        END DO
        ITMP = ITMP&
            *((SA_ST_IN*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))&
            *FluorYield(NELEMENT, SHELL(N))&
            *RadRate(NELEMENT, LINE(N))
        MM3 = ITMP
    END FUNCTION
END MODULE
