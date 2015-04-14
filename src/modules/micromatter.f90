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
        INTEGER :: CNT, CNT2
        INTEGER, INTENT(IN) :: NELEMENT
        INTEGER, INTENT(IN) :: N
        REAL(16) :: ITMP = 0
        REAL(16) :: PI
        REAL(16), INTENT(IN) :: KAPPA
        REAL(16)    :: I_CHAR
        REAL(16), DIMENSION(:), ALLOCATABLE :: TMP, EI
        ALLOCATE(TMP(SIZE(LINE)))
        ALLOCATE(EI(SIZE(LINE)))
        DATA ITMP/0/
        TMP = 0_16
        EI = 0._16
        PI = 2.D0*DASIN(1.D0)
        !WRITE (6,'(1H[,A16,2H](,A3,1H),I32,I32,ES32.20E3)') 'MM1', 'IN', N, NELEMENT, KAPPA
        DO CNT2 = 1, CP_ST%NELEMENTS
            I_CHAR = SECCOMPCHAR(N, CP_ST%ELEMENTS(CNT2), KAPPA)
            !WRITE (6,'(1H[,A16,2H](,A3,1H),6HICHAR=,ES32.20E3)') 'MM1', 'DEB', I_CHAR
            DO CNT = 2, 3!SIZE(LINE)
                EI(CNT) = LineEnergy(CP_ST%ELEMENTS(CNT2), LINE(CNT))
                !WRITE (6,'(1H[,A16,2H](,A3,1H),3HEI=ES32.20E3)') 'MM1', 'DEB', EI(CNT)
            END DO
            DO CNT = 1, SIZE(LINE)
                IF (EI(CNT).LE. 0) CYCLE
                IF (EI(CNT).LE. LineEnergy(NELEMENT, LINE(N))) CYCLE
                    TMP(CNT) = I_CHAR&
                        *CS_Photo_CP(STR_SAMPLE, DBLE(EI(CNT)))
            END DO
            ITMP = ITMP + SUM(TMP)&
                *((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))&
                *FLUORYIELD_CHG(NELEMENT, SHELL(N))&
                *RadRate(NELEMENT, LINE(N))
!WRITE (6,'(1H[,A16,2H](,A3,1H),5HITMP=,ES32.20E3)') 'MM1', 'DEB', ITMP
        END DO
        MM1 = ITMP
WRITE (6,'(1H[,A16,2H](,A3,1H),ES32.20E3)') 'MM1', 'OUT', MM1
        DEALLOCATE(TMP)
        DEALLOCATE(EI)
        RETURN
    END FUNCTION
    FUNCTION MM2(N, NELEMENT, KAPPA)
        IMPLICIT NONE
        REAL(16) :: MM2
        INTEGER :: CNT, CNT2, NELEMENT
        INTEGER :: LINES, SHELLS, N
        REAL(16) :: ITMP, TMPR, TMPC
        REAL(16), DIMENSION(:), ALLOCATABLE :: TMP, EI
        REAL(16) :: PI, KAPPA, EC
        DATA ITMP/0/
        ALLOCATE(TMP(SIZE(LINE)))
        ALLOCATE(EI(SIZE(LINE)))
        TMP = 0._16
        EI = 0._16
        PI = 2.D0*DASIN(1.D0)
        DO CNT = 1, SIZE(LINE)
            EI(CNT) = LineEnergy(Z_ANODE, LINE(CNT))
        END DO
        DO CNT = 2, 3!SIZE(LINE)
            IF (EI(CNT).LE. 0) CYCLE
            IF (EI(CNT).LE. LineEnergy(NELEMENT, LINE(N))) CYCLE
                TMPR = SECTRAYL(KAPPA, N)&
                    *CS_Photo_CP(STR_SAMPLE, DBLE(EI(CNT)))
                EC = ComptonEnergy(DBLE(EI(CNT)), DBLE(A_ST_AZIM_OUT))
                TMPC = SECTCOMP(KAPPA, N)&
                    *CS_Photo_CP(STR_SAMPLE, DBLE(EC))
            TMP(CNT) = TMPR + TMPC
        END DO
        ITMP = SUM(TMP)&
            *((SA_ST_OUT*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))&
            *FLUORYIELD_CHG(NELEMENT, SHELL(N))&
            *RadRate(NELEMENT, LINE(N))
        MM2 = ITMP
        DEALLOCATE(TMP)
        DEALLOCATE(EI)
        RETURN
    END FUNCTION
    FUNCTION MM3(N, NELEMENT, KAPPA)
        IMPLICIT NONE
        REAL(16) :: MM3
        INTEGER :: CNT, CNT2, Z, INTSTEP
        INTEGER :: LINES, SHELLS, N, ECNT, NELEMENT
        REAL(16) :: ITMP
        REAL(16) :: PI, KAPPA
        REAL(16), DIMENSION(:), ALLOCATABLE :: TMP, EI
        DATA ITMP/0/
        PI = 2.D0*DASIN(1.D0)
        INTSTEP = INT((VTUBE-EdgeEnergy(NELEMENT, SHELL(N)))/ESTEP)
        ALLOCATE(EI(0:INTSTEP))
        ALLOCATE(TMP(0:INTSTEP))
        EI = 0._16
        TMP = 0._16
        DO ECNT = 0, INTSTEP
            EI(ECNT) = EdgeEnergy(NELEMENT, SHELL(N)) + (DBLE(ECNT)*ESTEP)
        END DO
        DO ECNT = 0, INTSTEP
            TMP(ECNT) = SECTCONT(EI(ECNT))&
                *CS_Photo_CP(STR_SAMPLE, DBLE(EI(ECNT)))*ESTEP
        END DO
        ITMP = SUM(TMP)&
            *((SA_ST_IN*CONC)/(4*PI*SIN(A_ST_AZIM_IN)))&
            *FLUORYIELD_CHG(NELEMENT, SHELL(N))&
            *RadRate(NELEMENT, LINE(N))
        MM3 = ITMP
        DEALLOCATE(EI)
        DEALLOCATE(TMP)
    END FUNCTION
END MODULE
