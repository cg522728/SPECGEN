MODULE SECCOMP
    USE :: xraylib
    USE :: XRLDATA
    USE :: CFGDATA
    USE :: ANODE
    CONTAINS
    FUNCTION SECCOMPCHAR(N, ZELEMENT, KAPPA)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        IMPLICIT NONE
        REAL(16) :: ITMP, ITMP2, ITMP3
        REAL(16) :: EI, EST, E1, E2
        INTEGER :: CNT, ECNT, CMIN, N
        INTEGER :: ZELEMENT
        REAL(16) :: SECCOMPCHAR, CONST
        REAL(16) :: PI, KAPPA, TAU, MAC1, MAC2, FLUOR, RADR
        REAL(16), DIMENSION(:), ALLOCATABLE :: TMP1, TMP2, IA, IA2, EA
        PI = 2.D0*DASIN(1.D0)
        DATA ITMP/0/
        DATA ITMP2/0/
        DATA ITMP3/0/

        ALLOCATE(TMP1(SIZE(LINE)))
        ALLOCATE(TMP2(NSTEP))
        ALLOCATE(IA(SIZE(LINE)))
        ALLOCATE(EA(SIZE(LINE)))
        ALLOCATE(IA2(NSTEP))
        TMP1 = 0._16
        TMP2 = 0._16
        IA = 0._16
        EA = 0._16
        IA2 = 0._16

        EST = LineEnergy(ZELEMENT, LINE(N))
        IF (EST.EQ.0) THEN
            SECCOMPCHAR = SECTCONT(EST)
            RETURN
        ENDIF
        DO  CNT=1,SIZE(LINE)
            EA(CNT) = LineEnergy(Z_ANODE, LINE(CNT))
            IA(CNT) = ANODECHAR(CNT, KAPPA)
        END DO
        CONST = CS_Photo_CP(STR_SECTARGET, DBLE(EST))&
                    *FluorYield(ZELEMENT, SHELL(N))&
                    *RadRate(ZELEMENT, LINE(N))
        MAC1 = MAC_COMP(CP_ST, EST)/SIN(A_ST_AZIM_OUT)
        !$OMP PARALLEL SHARED(EA, IA, TMP1, TMP2)
        !$OMP DO
        DO  CNT=1,SIZE(LINE)
            IF (EA(CNT).EQ.0) CYCLE
            IF (EA(CNT).LT. EST) CYCLE
            IF (ISNAN(IA(CNT))) CYCLE
            TMP1(CNT) = IA(CNT)&
                /((MAC_COMP(CP_ST, EA(CNT))&
                /SIN(A_ST_AZIM_IN))&
                +(MAC1))
        END DO
        !$OMP END DO
        CONST = FLUORYIELD_CHG(ZELEMENT, SHELL(N))&
                   *RadRate(Z_ANODE, LINE(N))
        !$OMP DO
        DO ECNT= 1,NSTEP
            EI = EMIN+ESTEP*DBLE(ECNT)
            E1 = EI - (ESTEP/2)
            E2 = EI + (ESTEP/2)
            IF (ISNAN(ANODECONT(E1))) CYCLE
            IF (ISNAN(ANODECONT(E2))) CYCLE
            TMP2(ECNT) = ((ANODECONT(E1)&
                *CS_Photo_CP(STR_SECTARGET, DBLE(E1))&
                *CONST&
                /((MAC_COMP(CP_ST, E1)&
                /SIN(A_ST_AZIM_IN))&
                +MAC1))&
                +(ANODECONT(E2)&
                *CS_Photo_CP(STR_SECTARGET, DBLE(E2))&
                *CONST&
                /((MAC_COMP(CP_ST, E2)&
                /SIN(A_ST_AZIM_IN))&
                +MAC1)))*ESTEP
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        ITMP = SUM(TMP1)*(SA_ST_IN/(4*PI*SIN(A_ST_AZIM_IN)))
        ITMP3 = SUM(TMP2)*ESTEP*(SA_ST_IN/(4*PI*SIN(A_ST_AZIM_IN)))
        SECCOMPCHAR = (ITMP + ITMP3)
        RETURN
    END FUNCTION
    FUNCTION SECTCONT(EI)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        USE :: ANODE
        IMPLICIT NONE
        REAL(16) :: SECTCONT
        REAL(16) :: EI, EINCOH
        REAL(16) :: IRAY, ICOMPT, ITMP
        REAL(16) :: CRAY, CCOMPT
        REAL(16) :: PI

        DATA    IRAY/0/
        DATA    ITMP/0/
        PI = 2.D0*DASIN(1.D0)
        IRAY = (ANODECONT(EI)&
            *DCSP_Rayl_CP(STR_SECTARGET, DBLE(EI), DBLE(A_ST_POL), DBLE((PI/2)-A_ST_AZIM_OUT)))
        CRAY = ((MAC_COMP(CP_ST, EI))/SIN((PI/2)-A_ST_AZIM_IN))&
            +(MAC_COMP(CP_ST, EI)/SIN((PI/2)-A_ST_AZIM_OUT))
        IRAY = (IRAY/CRAY)
        EINCOH = ComptonEnergy(DBLE(EI), DBLE((PI/2)-A_ST_AZIM_OUT))
        ICOMPT = (ANODECONT(EI)&
            *DCSP_Compt_CP(STR_SECTARGET, DBLE(EINCOH), DBLE(A_ST_POL), DBLE((PI/2)-A_ST_AZIM_OUT)))
        CCOMPT= ((MAC_COMP(CP_ST, EINCOH))/SIN((PI/2)-A_ST_AZIM_IN))&
            +(MAC_COMP(CP_ST, EI)/SIN((PI/2)-A_ST_AZIM_OUT))
        ICOMPT = (ICOMPT/CCOMPT)
        DATA ITMP/0/
        ITMP = (IRAY + ICOMPT)&
            *(SA_ST_IN/SIN((PI/2)-A_ST_AZIM_IN))
        SECTCONT = ITMP*ESTEP
    END FUNCTION SECTCONT
    FUNCTION SECTRAYL(KAPPA, N)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        USE :: ANODE
        IMPLICIT NONE
        REAL(16)    :: SECTRAYL
        REAL(16)    :: PI, KAPPA
        REAL(16)    :: I_AN, I_RAYL, CON
        REAL(16)    :: EA
        INTEGER     :: N

        PI = 2.D0*DASIN(1.D0)

        EA = LineEnergy(Z_ANODE, LINE(N))
        I_AN = ANODECHAR(N, KAPPA)
        I_RAYL = I_AN*DCSP_Rayl_CP(STR_SECTARGET, DBLE(EA), DBLE(A_ST_POL), DBLE((PI/2)-A_ST_AZIM_OUT))
        !I_RAYL = I_AN*DCS_Rayl_CP(STR_SECTARGET, DBLE(EA), DBLE(A_ST_POL))
        CON = ((MAC_COMP(CP_ST, EA))/SIN((PI/2)-A_ST_AZIM_IN))&
            +(MAC_COMP(CP_ST, EA)/SIN((PI/2)-A_ST_AZIM_OUT))
        I_RAYL = (I_RAYL/CON)*(SA_ST_IN/SIN((PI/2)-A_ST_AZIM_IN))
        SECTRAYL = I_RAYL
        RETURN
    END FUNCTION SECTRAYL
    FUNCTION SECTCOMP(KAPPA, N)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        USE :: ANODE
        IMPLICIT NONE
        REAL(16)    :: SECTCOMP
        REAL(16)    :: PI, KAPPA
        REAL(16)    :: I_AN, I_COMP, CON
        REAL(16)    :: EI, EA
        INTEGER     :: N

        PI = 2.D0*DASIN(1.D0)

        EA = LineEnergy(Z_ANODE, LINE(N))
        EI = ComptonEnergy(DBLE(EA), DBLE((PI/2)-A_ST_AZIM_OUT))
        I_AN = ANODECHAR(N, KAPPA)
        I_COMP = I_AN*DCSP_Compt_CP(STR_SECTARGET, DBLE(EI), DBLE(A_ST_POL), DBLE((PI/2)-A_ST_AZIM_OUT))
        !I_COMP = I_AN*DCS_Compt_CP(STR_SECTARGET, DBLE(EI), DBLE(A_ST_POL))
        CON = ((MAC_COMP(CP_ST, EA))/SIN((PI/2)-A_ST_AZIM_IN))&
            +(MAC_COMP(CP_ST, EI)/SIN((PI/2)-A_ST_AZIM_OUT))
        I_COMP = (I_COMP/CON)*(SA_ST_IN/SIN((PI/2)-A_ST_AZIM_IN))
        SECTCOMP = I_COMP
        RETURN
    END FUNCTION SECTCOMP
END MODULE SECCOMP
