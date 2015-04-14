MODULE CONSTANTS
    USE     :: ISO_FORTRAN_ENV
    USE :: xraylib
    USE :: XRLDATA
    USE :: CFGDATA

CONTAINS
    FUNCTION PELLA_CALC_F(EI)
        IMPLICIT NONE
        REAL(16) :: PELLA_CALC_F
        REAL(16) :: CONST
        REAL(16) :: EI
        REAL(16) :: F
        REAL(16) :: I
        REAL(16) :: ETA, ETACONST
        REAL(16) :: C, C1, C2, C3

        ETA = VTUBE**(1.65)-EI**(1.65)
        ETACONST = MAC(Z_ANODE, EI)/(SIN(A_TAKE_OFF)*63.67)
        ETA = ETA*ETACONST
        C1 = 1+(2.56E-3)*Z_ANODE**2
        C2 = 1+3.17E4*(VTUBE**(-1))*(Z_ANODE**(-2))
        C3 = 0.25*ETA+1E4
        C = (1+C1**(-1))*(C2**(-1)*C3**(-1))
        F = (1+C*ETA)**(-2)
        PELLA_CALC_F = F
    END FUNCTION PELLA_CALC_F
    FUNCTION CALC_KR(Z, EI, SHELLS, KAPPA)
    !#################################################################################
    !#THIS FUNCTION CALCULATES THE K_r FACTOR IN THE PELLA ALGORITHM, WHICH WAS      #
    !#SIMPLIFIED IN THE ORIGINAL ALGORITHM BY A POLYNOME-FIT OF THE ACTUAL VALUES    #
    !#INPUTS ARE:                                                                    #
    !#      -Z      (INT)       ATOMIC NUMBER                                        #
    !#      -EI     (REAL16)    ENERGY IN keV                                        #
    !#      -SHELLS (INT)       INDEX OF SHELL-ARRAY IN MODULE XRLDATA               #
    !#      -KAPPA  (REAL16)    PROPORTIONAL CONSTANT CALCULATED WITH CHECKKAPPA()   #
    !#################################################################################
        IMPLICIT NONE
        REAL(16) :: CALC_KR
        INTEGER :: Z, SHELLS
        REAL(16) :: EI, UZ
        REAL(16) :: K, PI, P, OMEGAQ, R, A, C, LAMBDAL
        REAL(16)    :: KAPPA
        REAL(16), DIMENSION(5)   :: RCON

        DATA K/2.72E-6/
        DATA C/4.4E5/
        DATA RCON/1, 0.008157, 3.613E-5, 0.009583, 0.001141/


        UZ = VTUBE/EI
        LAMBDAL = KEV2ANGST/LineEnergy(Z, LINE(SHELLS))
        A = AtomicWeight(Z)
        R = RCON(1) - RCON(2)*Z + RCON(3)*(Z**2)+RCON(4)*Z*EXP(-UZ)+RCON(5)*VTUBE
        P = 1.62E-13*(Z-2)**2*Z*A*C*R**(-1)
        OMEGAQ = FLUORYIELD_CHG(Z, SHELL(SHELLS))
        PI = RadRate(Z, LINE(SHELLS))
        CALC_KR = (KAPPA/K)&
                    *PI&
                    *(1+P)&
                    *OMEGAQ&
                    *R&
                    *(LAMBDAL**2)&
                    *((A*C*Z)**(-1))
        RETURN
    END FUNCTION
    FUNCTION DETEFF(EI)
        REAL(16)   :: DETEFF
        REAL(16), INTENT(IN)    :: EI
        REAL(16)    :: K1, K2
        K1 = MAC(Z_DET_WINDOW, EI)*D_DET_WINDOW*1E-4&
                +MAC(Z_DET_DL, EI)*D_DET_DL*1E-4&
                +MAC(Z_DET_GAP, EI)*D_DET_GAP*1E-4
        K2 = MAC(Z_DET_BODY, EI)*D_DET_BODY*1E-4
        DETEFF = EXP(-K1)*(1-EXP(-K2))
        RETURN
    END FUNCTION
    FUNCTION STOPPINGFACTOR(EI, N)
        IMPLICIT NONE
        REAL(16)    :: STOPPINGFACTOR
        REAL(16)    :: EI, UZ
        INTEGER     :: N
        REAL(16)    :: Z, B, J, S, TMP, TMP2, TMP3

        IF (SHELL(N).EQ. K_SHELL) THEN
            Z = 2
            B = 0.35_16
        ENDIF
        IF (SHELL(N).EQ. L3_SHELL&
                .OR. SHELL(N).EQ. L2_SHELL&
                .OR. SHELL(N).EQ. L1_SHELL) THEN
            Z = 8
            B = 0.25_16
        ENDIF
                IF (SHELL(N).EQ. M5_SHELL&
                .OR. SHELL(N).EQ. M4_SHELL&
                .OR. SHELL(N).EQ. M3_SHELL&
                .OR. SHELL(N).EQ. M2_SHELL&
                .OR. SHELL(N).EQ. M1_SHELL) THEN
            Z = 16
            B = 0.15_16
        ENDIF
        UZ = VTUBE/EdgeEnergy(Z_ANODE, SHELL(N))
        J = 0.0135*Z_ANODE
        TMP = UZ*LOG(UZ)+1-UZ
        TMP2 = SQRT(UZ)*LOG(UZ)+2*(1-SQRT(UZ))
        TMP3 = 16.05*SQRT(J/EdgeEnergy(Z_ANODE, SHELL(N)))*(TMP2/TMP)
        S = ((Z*B)/Z_ANODE)*TMP*(1+TMP3)
        STOPPINGFACTOR = 1/S
        RETURN
    END FUNCTION STOPPINGFACTOR
END MODULE CONSTANTS
