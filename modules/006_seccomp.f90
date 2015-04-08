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
        REAL(16) :: EI, EST, EA
        INTEGER :: CNT, ECNT, CMIN, N
        INTEGER :: ZELEMENT
        REAL(16) :: SECCOMPCHAR
        REAL(16) :: PI, KAPPA, TAU, MAC1, MAC2, FLUOR, RADR
        PI = 2.D0*DASIN(1.D0)
        DATA ITMP/0/
        DATA ITMP2/0/
        DATA ITMP3/0/

        EST = LineEnergy(ZELEMENT, LINE(N))
        IF (EST.EQ.0) THEN
            SECCOMPCHAR = 0.
            RETURN
        ENDIF
        DO  CNT=1,SIZE(LINE)
            EA = LineEnergy(Z_ANODE, LINE(CNT))
            IF (EA.EQ.0) CYCLE
            IF (EA.LT. EST) CYCLE
            ITMP = ANODECHAR(CNT, KAPPA)&
                *CS_Photo_CP(STR_SECTARGET, DBLE(EST))&
                *FluorYield(ZELEMENT, SHELL(N))&
                *RadRate(ZELEMENT, LINE(N))&
                /((MAC_COMP(CP_ST, EA)&
                /SIN(A_ST_AZIM_IN))&
                +(MAC_COMP(CP_ST, EST)&
                /SIN(A_ST_AZIM_OUT))) + ITMP
        END DO
        ITMP = ITMP*(SA_ST_IN/(4*PI*SIN(A_ST_AZIM_IN)))
        DO ECNT= 1,NSTEP
            EI = EMIN+ESTEP*DBLE(ECNT)
            ITMP2 = ANODECONT(EI)&
                *CS_Photo_CP(STR_SECTARGET, DBLE(EI))&
                *FluorYield(ZELEMENT, SHELL(N))&
                *RadRate(Z_ANODE, LINE(N))&
                /((MAC_COMP(CP_ST, EI)&
                /SIN(A_ST_AZIM_IN))&
                +(MAC_COMP(CP_ST, EST)&
                /SIN(A_ST_AZIM_OUT)))
            ITMP2 = ITMP2 + ITMP3
            ITMP3 = ITMP2
2       END DO
        ITMP3 = ITMP3*ESTEP*(SA_ST_IN/(4*PI*SIN(A_ST_AZIM_IN)))
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
