MODULE ANODE
    USE :: xraylib
    USE :: TYPES
    USE :: CFGDATA
    USE :: XRLDATA
    USE :: CONSTANTS
    USE :: MATHCHG
    IMPLICIT NONE
    PRIVATE
    PUBLIC ANODE_CONT, ANODE_CHAR, TUBE_ATTEN
CONTAINS
    FUNCTION ANODE_CONT(EI)
        REAL(QP)    :: ANODE_CONT
        REAL(DP), INTENT(IN)    :: EI

        REAL(DP)    :: E1 = 0_DP
        REAL(DP)    :: E2 = 0_DP

        E1 = EI - (ESTEP/2)
        E2 = EI + (ESTEP/2)
        ANODE_CONT = INTEGRATE(DERIV_ANODE_CONT, E1, E2, INT(25), Z_ANODE)
        RETURN
    END FUNCTION ANODE_CONT
    FUNCTION ANODE_CHAR(N, Z_INT)
        IMPLICIT NONE
        REAL(QP) :: ANODE_CHAR
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN):: Z_INT
        
        REAL(DP) :: E_CHAR = 0_DP
        REAL(DP)    :: E_EDGE = 0_DP
        REAL(QP) :: R_BACKSCATTER = 0_QP
        REAL(DP) :: UZ = 0_DP
        REAL(DP)    :: ETUBE = 0_DP
        REAL(QP), DIMENSION(5)   :: RCON
        REAL(DP)    :: Z = 0_DP

        DATA RCON/1, 0.0081517, 3.613E-5, 0.009583, 0.001141/

        Z = DBLE(Z_INT)
        E_CHAR = LineEnergy(Z_INT, LINE(N))
        E_EDGE = EdgeEnergy(Z_INT, SHELL(N))
        IF (E_CHAR.EQ. 0) THEN
            ANODE_CHAR = 0
            RETURN
        ENDIF
        ETUBE = DBLE(VTUBE)
        UZ = ETUBE/E_EDGE
        R_BACKSCATTER = RCON(1)&
                            - RCON(2)*Z&
                            + RCON(3)*(Z**2)&
                            +RCON(4)*Z*EXP(-UZ)&
                            +RCON(5)*ETUBE
        ANODE_CHAR = SA_ANODE_OUT&
                    *ITUBE&
                    *EBEL_CONST(N, Z_INT)&
                    *R_BACKSCATTER&
                    *STOPPINGFACTOR(N, Z_INT)&
                    *CALC_F(E_EDGE, Z_INT, N)&
                    *TRANSPROB(Z_INT, LINE(N))&
                    *FLUORYIELD_CHG(Z_INT, SHELL(N))!&
                    !*FLUORYIELD(Z_INT, SHELL(N))
        RETURN
    END FUNCTION ANODE_CHAR
    FUNCTION DERIV_ANODE_CONT(EI, Z)
        IMPLICIT NONE
        REAL(QP) :: DERIV_ANODE_CONT
        REAL(DP), INTENT(IN)    :: EI
        INTEGER , INTENT(IN)    :: Z
        
        REAL(QP) :: CONST = 0_QP
        REAL(QP)    :: EBEL = 1.37E9_QP
        REAL(QP) :: X = 0_QP
        REAL(QP) :: N_CONT = 0_QP
        REAL(DP)    :: UZ = 0_DP
        REAL(DP)    :: ETUBE

        ETUBE = DBLE(VTUBE)
        UZ = ETUBE/EI
        IF (EI.GE.ETUBE) THEN
            DERIV_ANODE_CONT = 0_QP
            RETURN
        ENDIF
        CONST = SA_ANODE_OUT*ITUBE*EBEL*Z
        X = 1.0314-0.0032*Z+0.0047*VTUBE
        N_CONT = CONST*((UZ-1)**X)*CALC_F(EI, Z)
        DERIV_ANODE_CONT = N_CONT
        RETURN
    END FUNCTION DERIV_ANODE_CONT
    FUNCTION TUBE_ATTEN(EI)
        IMPLICIT NONE
        REAL(QP)    :: TUBE_ATTEN
        REAL(DP), INTENT(IN)    :: EI

        REAL(QP)    :: ATTEN_WINDOW = 0_QP
        REAL(QP)    :: ATTEN_FILTER = 0_QP

        ATTEN_WINDOW = EXP(-MAC(Z_WINDOW, EI)*D_WINDOW*ElementDensity(Z_WINDOW)*1E-4)
        IF (D_FILTER.GT. 0) THEN
            ATTEN_FILTER = EXP(-MAC(Z_FILTER, EI)*D_FILTER*ElementDensity(Z_FILTER)*1E-4)
        ELSE
            ATTEN_FILTER = 1_QP
        ENDIF
        TUBE_ATTEN = ATTEN_WINDOW*ATTEN_FILTER
    END FUNCTION TUBE_ATTEN
    FUNCTION CALC_F(EI, Z, N)
        IMPLICIT NONE
        REAL(QP)    :: CALC_F
        REAL(DP), INTENT(IN)    :: EI
        INTEGER, INTENT(IN)     :: Z
        INTEGER, INTENT(IN), OPTIONAL   :: N

        REAL(QP)    :: TAU = 0_QP
        REAL(QP)    :: TMP = 0_QP

        IF (PRESENT(N)) TAU = CS_Photo(Z, EI)
        IF (.NOT.PRESENT(N)) TAU = CS_TOTAL(Z, EI)
        TMP = TAU*2*DDF(EI, Z)*(SIN(A_INCID)/SIN(A_TAKE_OFF))
        CALC_F = (1-EXP(-TMP))/TMP
        RETURN
    END FUNCTION CALC_F
        FUNCTION DDF(EI, Z_INT)
        REAL(QP)    ::DDF
        REAL(DP), INTENT(IN)    :: EI
        INTEGER, INTENT(IN)     :: Z_INT

        REAL(QP)    :: ETA = 0_QP
        REAL(QP)    :: TMP1 = 0_QP
        REAL(QP)    :: TMP2 = 0_QP
        REAL(QP)    :: TMP3 = 0_QP
        REAL(QP)    :: M = 0_QP
        REAL(QP)    :: UZ = 0_QP
        REAL(QP)    :: J = 0_QP
        REAL(QP)    :: EC = 0_QP
        REAL(QP)    :: Z = 0_QP
        REAL(QP)    :: A = 0_QP
        REAL(DP)    :: ETUBE = 0_DP
        REAL(QP)    :: LOGZ = 0_QP
        REAL(QP), DIMENSION(2)    :: MCON
        REAL(QP), DIMENSION(4)    :: ECCON
        REAL(QP), DIMENSION(8)    :: DCON

        DATA MCON/0.1382, 0.9211/
        DATA ECCON/0.1904, 0.2236, 0.1292, 0.0149/
        DATA DCON/0.787E-5, 0.735E-6,&
                 0.49269, 1.0987, 0.78557,&
                 0.70256, 1.09865, 1.0046/

        Z = DBLE(Z_INT)
        A = AtomicWeight(Z_INT)
        ETUBE = DBLE(VTUBE)
        UZ = ETUBE/EI
        M = MCON(1) - (MCON(2)*(Z**(-0.5)))
        J = 0.0135_QP*Z
        LOGZ = LOG(Z)
        EC = ECCON(1)&
                - ECCON(2)*LOGZ&
                +ECCON(3)*(LOGZ**2)&
                -ECCON(4)*(LOGZ**3)
        ETA = (ETUBE**M)*EC
        TMP1 = (A/Z)*(DCON(1)*SQRT(J)*(ETUBE**(1.5_QP))+DCON(2)*(ETUBE**(2_DP)))
        TMP2 = DCON(3)-DCON(4)*ETA+DCON(5)*ETA**2_DP
        TMP3 = DCON(6)-DCON(7)*ETA+DCON(8)*(ETA**2_DP)+LOG(UZ)
        DDF = TMP1*(TMP2/TMP3)*LOG(UZ)
        RETURN
    END FUNCTION DDF
    FUNCTION STOPPINGFACTOR(N, Z)
        IMPLICIT NONE
        REAL(QP)    :: STOPPINGFACTOR
        INTEGER, INTENT(IN)     :: N
        INTEGER, INTENT(IN)     :: Z

        REAL(QP)    :: UZ = 0_QP
        REAL(QP)    :: ZCON = 0_QP
        REAL(QP)    :: BCON = 0_QP
        REAL(QP)    :: J = 0_QP
        REAL(QP)    :: S = 0_QP
        REAL(QP)    :: TMP = 0_QP
        REAL(QP)    :: TMP2 = 0_QP
        REAL(QP)    :: TMP3 = 0_QP
        REAL(DP)    :: ETUBE
        REAL(DP)    :: E_EDGE

        IF (SHELL(N).EQ. K_SHELL) THEN
            ZCON = 2_QP
            BCON = 0.35_QP
        ENDIF
        IF (SHELL(N).EQ. L3_SHELL&
                .OR. SHELL(N).EQ. L2_SHELL&
                .OR. SHELL(N).EQ. L1_SHELL) THEN
            ZCON = 8_QP
            BCON = 0.25_QP
        ENDIF
                IF (SHELL(N).EQ. M5_SHELL&
                .OR. SHELL(N).EQ. M4_SHELL&
                .OR. SHELL(N).EQ. M3_SHELL&
                .OR. SHELL(N).EQ. M2_SHELL&
                .OR. SHELL(N).EQ. M1_SHELL) THEN
            ZCON = 1.2_QP
            BCON = 0.7_QP
        ENDIF
        ETUBE = DBLE(VTUBE)
        E_EDGE = EdgeEnergy(Z, SHELL(N))
        UZ = ETUBE/E_EDGE
        J = 0.0135_QP*DBLE(Z)
        TMP = UZ*LOG(UZ)+1-UZ
        TMP2 = SQRT(UZ)*LOG(UZ)+2*(1-SQRT(UZ))
        TMP3 = 16.05_QP*SQRT(J/E_EDGE)*(TMP2/TMP)
        S = ((ZCON*BCON)/DBLE(Z))*TMP*(1+TMP3)
        STOPPINGFACTOR = S
        RETURN
    END FUNCTION STOPPINGFACTOR
    FUNCTION EBEL_CONST(N, Z)
        REAL(QP)    :: EBEL_CONST
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: Z

        REAL(QP)    :: FCOR = 0_QP
        REAL(DP)    :: UZ = 0_DP
        REAL(DP)    :: E_EDGE = 0_DP
        REAL(DP)    :: E_TUBE = 0_DP
        REAL(QP)    :: TMP1 = 0_QP
        REAL(QP)    :: TMP2 = 0_QP
        REAL(QP), DIMENSION(4)    :: KSHELL
        REAL(QP)    :: L1SHELL
        REAL(QP)    :: L2SHELL
        REAL(QP)    :: L3SHELL
        REAL(QP)    :: M1SHELL
        REAL(QP)    :: M2SHELL
        REAL(QP)    :: M3SHELL
        REAL(QP)    :: M4SHELL
        REAL(QP)    :: M5SHELL
        REAL(QP), DIMENSION(3)    :: FCON

        DATA KSHELL/4.697, 0.134, 0.00268, 1E13/
        DATA L1SHELL/0.72E13/
        DATA L2SHELL/2.70E13/
        DATA L3SHELL/4.94E13/
        DATA M1SHELL/0/
        DATA M2SHELL/0/
        DATA M3SHELL/2.32E13/
        DATA M4SHELL/15.8E13/
        DATA M5SHELL/20.5E13/
        DATA FCON/0.4814, 0.03781, 2.413E-4/

        E_TUBE = DBLE(VTUBE)
        E_EDGE = EdgeEnergy(Z, SHELL(N))
        UZ = E_TUBE/E_EDGE
        IF (Z.LT.80) THEN
            FCOR = -FCON(1)&
                    + FCON(2)*DBLE(Z)&
                    - FCON(3)*(DBLE(Z)**2)
        ELSE
            FCOR = 1_QP
        ENDIF

        IF (SHELL(N).EQ. K_SHELL) THEN
            FCOR = 1
            TMP1 = KSHELL(1) + KSHELL(2)*UZ - KSHELL(3)*(UZ**2)
            TMP2 = 1-EXP(-7_QP*(UZ-1_QP))
            EBEL_CONST = KSHELL(4)*(TMP1/TMP2)
            RETURN
        ELSEIF (SHELL(N).EQ.L1_SHELL) THEN
            EBEL_CONST = FCOR*L1SHELL
            RETURN
        ELSEIF (SHELL(N).EQ.L2_SHELL) THEN
            EBEL_CONST = FCOR*L2SHELL
            RETURN
        ELSEIF (SHELL(N).EQ.L3_SHELL) THEN
            EBEL_CONST = L3SHELL
            RETURN
        ELSEIF (SHELL(N).EQ.M1_SHELL) THEN
            EBEL_CONST = M1SHELL
            RETURN
        ELSEIF (SHELL(N).EQ.M2_SHELL) THEN
            EBEL_CONST = M2SHELL
            RETURN
        ELSEIF (SHELL(N).EQ.M3_SHELL) THEN
            EBEL_CONST = M3SHELL
            RETURN
        ELSEIF (SHELL(N).EQ.M4_SHELL) THEN
            EBEL_CONST = M4SHELL
            RETURN
        ELSEIF (SHELL(N).EQ.M5_SHELL) THEN
            EBEL_CONST = M5SHELL
            RETURN
        ELSE
            EBEL_CONST = 0_QP
            RETURN
        ENDIF
    END FUNCTION EBEL_CONST
END MODULE ANODE
