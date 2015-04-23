MODULE SECCOMP
    USE :: xraylib
    USE :: TYPES
    USE :: XRLDATA
    USE :: CFGDATA
    USE :: ANODE
    USE :: MATHCHG
    implicit none
    PRIVATE
    PUBLIC SECT_CHAR, SECT_CONT, AN_SCAT_RAYL, AN_SCAT_COMP
    CONTAINS
    FUNCTION SECT_CHAR(N, ZELEMENT)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        IMPLICIT NONE
        REAL(QP) :: SECT_CHAR
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: ZELEMENT

        REAL(QP)    :: ITMP1 = 0_QP
        REAL(QP)    :: ITMP2 = 0_QP
        REAL(DP)    :: E_CHAR_ST = 0_DP
        REAL(DP)    :: E_EDGE_ST = 0_DP
        REAL(DP)    :: E1 = 0_DP
        REAL(DP)    :: E2 = 0_DP
        INTEGER     :: CNT
        REAL(QP)    :: PI = 0_QP
        REAL(QP)    :: MAC1 = 0_QP
        REAL(QP), DIMENSION(:), ALLOCATABLE :: TMP
        REAL(QP), DIMENSION(:), ALLOCATABLE :: IA
        REAL(QP), DIMENSION(:), ALLOCATABLE :: MAC2
        REAL(QP), DIMENSION(:), ALLOCATABLE :: CONST
        REAL(DP), DIMENSION(:), ALLOCATABLE :: EA

        ALLOCATE(TMP(SIZE(LINE)))
        ALLOCATE(IA(SIZE(LINE)))
        ALLOCATE(EA(SIZE(LINE)))
        ALLOCATE(CONST(SIZE(LINE)))
        ALLOCATE(MAC2(SIZE(LINE)))

        PI = 2.D0*DASIN(1.D0)
        TMP = 0_QP
        IA = 0_QP
        EA = 0_DP
        CONST = 0_QP
        MAC2 = 0_QP

        E_CHAR_ST = LineEnergy(ZELEMENT, LINE(N))
        E_EDGE_ST = EdgeEnergy(ZELEMENT, SHELL(N))
        IF (E_CHAR_ST.EQ.0) THEN
            SECT_CHAR = 0_QP
            RETURN
        ENDIF
        MAC1 = MAC_COMP(CP_ST, E_CHAR_ST)
        DO  CNT=1,SIZE(LINE)
            EA(CNT) = LineEnergy(Z_ANODE, LINE(CNT))
            IF (EA(CNT) .EQ. 0) CYCLE
            IF (EA(CNT).LT.EMIN .OR. EA(CNT).LT.E_EDGE_ST) THEN
                IA(CNT) = 0_QP
            ELSE
                IA(CNT) = ANODE_CHAR(CNT, Z_ANODE)
            ENDIF
            CONST(CNT) = CS_FLUOR_CP_CHG(CP_ST, N, EA(CNT))
            MAC2(CNT) = MAC_COMP(CP_ST, EA(CNT))
        END DO

        DO  CNT=1,SIZE(LINE)
            IF (EA(CNT).EQ.0 .OR. EA(CNT).LT. E_CHAR_ST .OR. EA(CNT).LT. E_EDGE_ST) CYCLE
            TMP(CNT) = IA(CNT)*CONST(CNT)&
                /((MAC2(CNT)/SIN(A_ST_AZIM_IN))&
                +(MAC1/SIN(A_ST_AZIM_OUT)))
        END DO
        ITMP1 = SUM(TMP)*(SA_ST_IN/(4*PI*SIN(A_ST_AZIM_IN)))

        DEALLOCATE(EA)
        DEALLOCATE(IA)
        DEALLOCATE(TMP)
        DEALLOCATE(MAC2)
        DEALLOCATE(CONST)

        E1 = E_EDGE_ST
        E2 = DBLE(VTUBE)
        ITMP2 = INTEGRATE(DERIV_SECT_CHAR, E1, E2, NSTEP, ZELEMENT, N)&
                    *(SA_ST_IN/(4*PI*SIN(A_ST_AZIM_IN)))
        SECT_CHAR = (ITMP1 + ITMP2)
        RETURN
    END FUNCTION
    FUNCTION DERIV_SECT_CHAR(EI, Z, N)
        IMPLICIT NONE
        REAL(QP)    :: DERIV_SECT_CHAR
        REAL(DP), INTENT(IN)    :: EI
        INTEGER, INTENT(IN)     :: Z
        INTEGER, INTENT(IN)     :: N

        REAL(DP)    ::  E_EDGE = 0_DP
        REAL(DP)    ::  E_CHAR = 0_DP
        REAL(QP)    ::  N_CHAR = 0_QP

            E_EDGE = EdgeEnergy(Z, SHELL(N))
            E_CHAR = LineEnergy(Z, LINE(N))
            N_CHAR = ANODE_CONT(EI)&
                *CS_FLUOR_CP_CHG(CP_ST, N, EI)&
                /((MAC_COMP(CP_ST, EI)/SIN(A_ST_AZIM_IN))&
                +(MAC_COMP(CP_ST, E_CHAR)/SIN(A_ST_AZIM_OUT)))
            DERIV_SECT_CHAR = N_CHAR
            RETURN
    END FUNCTION DERIV_SECT_CHAR
    FUNCTION SECT_CONT(EI)
        IMPLICIT NONE
        REAL(QP)    :: SECT_CONT
        REAL(DP), INTENT(IN)    :: EI

        REAL(QP)    :: N_CONT = 0_QP
        REAL(QP)    :: TMP = 0_QP
        REAL(DP)    :: E1 = 0_DP
        REAL(DP)    :: E2 = 0_DP

        E1 = EI
        E2 = EI+ESTEP
        TMP = INTEGRATE(DERIV_SECT_CONT, E1, E2, INT(5))
        N_CONT = TMP*(SA_ST_IN/SIN(A_ST_AZIM_IN))
        SECT_CONT = N_CONT
        RETURN
    END FUNCTION SECT_CONT
    FUNCTION DERIV_SECT_CONT(EI)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        USE :: ANODE
        IMPLICIT NONE
        REAL(QP)    :: DERIV_SECT_CONT
        REAL(DP), INTENT(IN)    :: EI

        REAL(DP)    :: EINCOH = 0_DP
        REAL(QP)    :: IRAY = 0_QP
        REAL(QP)    :: ICOMPT = 0_QP
        REAL(QP)    :: ITMP = 0_QP
        REAL(QP)    :: CRAY = 0_QP
        REAL(QP)    :: CCOMPT = 0_QP
        REAL(QP)    :: PI = 0_QP
        REAL(QP)    :: DCSR = 0_QP
        REAL(QP)    :: DCSC = 0_QP
        INTEGER     :: CNT

        DATA    IRAY/0/
        DATA    ITMP/0/
        IF (EI.EQ.0) THEN
            DERIV_SECT_CONT = 0_QP
            RETURN
        ENDIF
        IF (EI.GT. VTUBE) THEN
            DERIV_SECT_CONT = 0_QP
            RETURN
        ENDIF

        DO CNT = 1, CP_ST%NELEMENTS
            DCSR = DCSR + DCSP_Rayl(CP_ST%ELEMENTS(CNT), EI, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))&
                    *CP_ST%massFractions(CNT)
        END DO

        PI = 2.D0*DASIN(1.D0)
        IRAY = ANODE_CONT(EI)*DCSR!DCSP_Rayl(Z_ANODE, EI, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))
        CRAY = ((MAC_COMP(CP_ST, EI))/SIN(A_ST_AZIM_IN))&
            +(MAC_COMP(CP_ST, EI)/SIN(A_ST_AZIM_OUT))
        IRAY = (IRAY/CRAY)
        EINCOH = ComptonEnergy(EI, DBLE(A_ST_AZIM_OUT))
        DO CNT = 1, CP_ST%NELEMENTS
            DCSC = DCSC + DCSP_CompT(CP_ST%ELEMENTS(CNT), EINCOH, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))&
                    *CP_ST%massFractions(CNT)
        END DO
        ICOMPT = ANODE_CONT(EI)*DCSC!DCSP_Compt(Z_ANODE, EINCOH, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))
        CCOMPT= ((MAC_COMP(CP_ST, EINCOH))/SIN(A_ST_AZIM_IN))&
            +(MAC_COMP(CP_ST, EI)/SIN(A_ST_AZIM_OUT))
        ICOMPT = (ICOMPT/CCOMPT)
        DERIV_SECT_CONT = (IRAY + ICOMPT)
        RETURN
    END FUNCTION DERIV_SECT_CONT
    FUNCTION AN_SCAT_RAYL(N)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        USE :: ANODE
        IMPLICIT NONE
        REAL(QP)    :: AN_SCAT_RAYL
        INTEGER, INTENT(IN)     :: N

        REAL(QP)    :: PI = 0_QP
        REAL(QP)    :: I_AN = 0_QP
        REAL(QP)    :: I_RAYL = 0_QP
        REAL(QP)    :: CON = 0_QP
        REAL(DP)    :: EA = 0_DP
        REAL(QP)    :: DCS = 0_QP
        INTEGER     :: CNT

        PI = 2.D0*DASIN(1.D0)

        EA = LineEnergy(Z_ANODE, LINE(N))
        I_AN = ANODE_CHAR(N, Z_ANODE)

        DO CNT = 1, CP_ST%NELEMENTS
            DCS = DCS + DCSP_Rayl(CP_ST%ELEMENTS(CNT), EA, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))&
                    *CP_ST%massFractions(CNT)
        END DO

        I_RAYL = I_AN*DCS!DCSP_Rayl(Z_ANODE, EA, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))
        CON = ((MAC_COMP(CP_ST, EA))/SIN(A_ST_AZIM_IN))&
            +(MAC_COMP(CP_ST, EA)/SIN(A_ST_AZIM_OUT))
        I_RAYL = (I_RAYL/CON)*(SA_ST_IN/SIN(A_ST_AZIM_IN))
        AN_SCAT_RAYL = I_RAYL
        RETURN
    END FUNCTION AN_SCAT_RAYL
    FUNCTION AN_SCAT_COMP(N)
        USE :: xraylib
        USE :: CFGDATA
        USE :: XRLDATA
        USE :: ANODE
        IMPLICIT NONE
        REAL(QP)    :: AN_SCAT_COMP
        INTEGER, INTENT(IN)     :: N

        REAL(QP)    :: PI = 0_QP
        REAL(QP)    :: I_AN = 0_QP
        REAL(QP)    :: I_COMP = 0_QP
        REAL(QP)    :: CON = 0_QP
        REAL(DP)    :: EI = 0_DP
        REAL(DP)    :: EA = 0_DP
        REAL(QP)    :: MAC1 = 0_QP
        REAL(QP)    :: MAC2 = 0_QP
        REAL(QP)    :: DCS = 0_QP
        INTEGER     :: CNT

        PI = 2.D0*DASIN(1.D0)

        EA = LineEnergy(Z_ANODE, LINE(N))
        EI = ComptonEnergy(EA, DBLE(A_ST_AZIM_OUT))
        I_AN = ANODE_CHAR(N, Z_ANODE)

        DO CNT = 1, CP_ST%NELEMENTS
            DCS = DCS + DCSP_Compt(CP_ST%ELEMENTS(CNT), EA, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))&
                    *CP_ST%massFractions(CNT)
        END DO

        I_COMP = I_AN*DCS!DCSP_Compt(Z_ANODE, EA, DBLE(A_ST_POL), DBLE(A_ST_AZIM_OUT))
        MAC1 = MAC_COMP(CP_ST, EA)
        MAC2 = MAC_COMP(CP_ST, EI)
        CON = (MAC1/SIN(A_ST_AZIM_IN))&
            +(MAC2/SIN(A_ST_AZIM_OUT))
        I_COMP = (I_COMP/CON)*(SA_ST_IN/SIN(A_ST_AZIM_IN))
        AN_SCAT_COMP = I_COMP
        RETURN
    END FUNCTION AN_SCAT_COMP
END MODULE SECCOMP
