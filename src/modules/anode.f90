MODULE ANODE
    USE :: xraylib
    USE :: TYPES
    USE :: XRLDATA
    USE :: CFGDATA
    USE :: MATHCHG
    IMPLICIT NONE
    PRIVATE
    PUBLIC ANODE_CONT, ANODE_CHAR, TUBE_ATTEN, FILTER_ATTEN
CONTAINS
    FUNCTION ANODE_CHAR(Z, N, OMEGA, I) RESULT(N_CHAR)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE NUMBER OF COUNTS FOR A CHARACTERISTIC LINE        #
        !#AS DEFINED IN:                                                                 #
        !#  EBEL, Horst. X‐ray tube spectra. X‐Ray Spectrometry, 1999, 28.4: 255-266.    #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)       ATOMIC NUMBER                   [-]                  #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#      -OMEGA  (WP)        SOLID ANGLE OF EXIT             [sr]                 #
        !#      -I      (WP)        TUBE CURRENT                    [mA]                 #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP), INTENT(IN) :: OMEGA
        REAL(WP), INTENT(IN) :: I
        REAL(WP) :: N_CHAR

        REAL(DP) :: E_LINE
        REAL(DP) :: E_EDGE
        REAL(WP) :: CONST
        REAL(WP) :: SPF
        REAL(WP) :: R
        REAL(WP) :: F

        E_LINE = LINE_ENERGY(Z, N)
        E_EDGE = EDGE_ENERGY(Z, N)
        CONST = EBEL_CONST(N, Z)
        SPF = STOPPING_POWER_FACTOR(E_LINE, Z, N)
        R = BACKSCAT_FACTOR(E_LINE ,Z)
        F = ABSORB_TERM(E_LINE, Z)

        N_CHAR = OMEGA&
                    *I&
                    *CONST&
                    *SPF&
                    *R&
                    *FLUOR_YIELD(Z, N)&
                    *TRANSPROB(Z, N)&
                    *F
        RETURN
    END FUNCTION ANODE_CHAR

    FUNCTION ANODE_CONT(Z, E, OMEGA, I) RESULT(N_CONT)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE NUMBER OF COUNTS AT ENERGY E IN THE CONTINUUM     #
        !#AS DEFINED IN:                                                                 #
        !#  EBEL, Horst. X‐ray tube spectra. X‐Ray Spectrometry, 1999, 28.4: 255-266.    #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)       ATOMIC NUMBER                   [-]                  #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -OMEGA  (WP)        SOLID ANGLE OF EXIT             [sr]                 #
        !#      -I      (WP)        TUBE CURRENT                    [mA]                 #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        REAL(DP), INTENT(IN) :: E
        REAL(WP), INTENT(IN) :: OMEGA
        REAL(WP), INTENT(IN) :: I
        REAL(WP) :: N_CONT

        REAL(WP) :: INTEGRAL
        REAL(DP) :: E_INIT
        REAL(DP) :: E_FINAL

        E_INIT = E
        E_FINAL = E + ESTEP
        INTEGRAL = INTEGRATE(DERIV_CONT, E_INIT, E_FINAL, INT(10), Z)
        N_CONT = OMEGA&
                    *I&
                    *INTEGRAL
        RETURN
    END FUNCTION ANODE_CONT

    FUNCTION DERIV_CONT(E, Z) RESULT(N_CONT)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE DERIVATIVE OF THE NUMBER OF COUNTS AT ENERGY      #
        !#E IN THE CONTINUUM AS DEFINED IN:                                              #
        !#  EBEL, Horst. X‐ray tube spectra. X‐Ray Spectrometry, 1999, 28.4: 255-266.    #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER                   [-]                  #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        REAL(DP), INTENT(IN) :: E
        REAL(WP) :: N_CONT

        INTEGER :: N = 0
        REAL(WP) :: CS_SR
        REAL(WP) :: F_ABS

        F_ABS = ABSORB_TERM(E, Z)
        CS_SR = CS_SMITH_REED(E, Z, N)
        N_CONT = CS_SR&
                    *F_ABS
        RETURN
    END FUNCTION DERIV_CONT

   FUNCTION BACKSCAT_COEF_GRAD(Z) RESULT(G)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE GRADIENT OF THE BACKSCATTER COEFFICIENT           #
        !#AS DEFINED IN:                                                                 #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D: Applied     #
        !#  Physics, 1978, 11.10: 1369.                                                  #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: G

        REAL(WP), DIMENSION(3) :: CON

        DATA CON/1112.8E-4, 30.289E-4, 0.15498E-4/

        G = -CON(1)&
                +CON(2)*Z&
                -CON(3)*(Z**2)
        RETURN
    END FUNCTION BACKSCAT_COEF_GRAD

    FUNCTION BACKSCAT_COEF_INTERCEPT(Z) RESULT(I)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTERCEPT OF THE BACKSCATTER COEFFICIENT          #
        !#AS DEFINED IN:                                                                 #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D: Applied     #
        !#  Physics, 1978, 11.10: 1369.                                                  #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: I

        REAL(WP), DIMENSION(4) :: CON

        DATA CON/52.3791E-4, 150.48371E-4, 1.67373E-4, 0.00716E-4/

        I = -CON(1)&
                +CON(2)*Z&
                -CON(3)*(Z**2)&
                +CON(4)*(Z**3)
        RETURN
    END FUNCTION BACKSCAT_COEF_INTERCEPT

    FUNCTION BACKSCAT_COEF(Z) RESULT(ETA)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE THE BACKSCATTER COEFFICIENT                       #
        !#AS DEFINED IN:                                                                 #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D: Applied     #
        !#  Physics, 1978, 11.10: 1369.                                                  #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: ETA

        REAL(WP) :: ETA20
        REAL(WP) :: GETA20
        REAL(WP) :: ECALC = 20_WP

        ETA20 = BACKSCAT_COEF_INTERCEPT(Z)
        GETA20 = BACKSCAT_COEF_GRAD(Z)

        ETA = ETA20*(1+GETA20*LOG(VTUBE/ECALC))
        RETURN
    END FUNCTION BACKSCAT_COEF

    FUNCTION MAX_PATH_LENGTH(Z) RESULT(RSM)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE THE MAXIMUM PATH LENGTH OF AN ELECTRON IN         #
        !#A TARGET AS DEFINED IN:                                                        #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D: Applied     #
        !#  Physics, 1978, 11.10: 1369.                                                  #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: RSM

        REAL(WP) :: IONPOT
        REAL(WP) :: ZARATIO

        REAL(WP), DIMENSION(2) :: CON

        DATA CON/0.773E-5, 0.735E-6/

        IONPOT = NIST_ION_POT(Z)
        ZARATIO = NIST_ZA_RATIO(Z)

        RSM = (CON(1)*SQRT(IONPOT*VTUBE)&
                    *VTUBE&
                +CON(2)*(VTUBE**2))&
                /ZARATIO
        RETURN
    END FUNCTION MAX_PATH_LENGTH

    FUNCTION MEAN_MASS_DEPTH(E, Z) RESULT(RS)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE MEAN MASS DEPTH OF X-RAY GENERATION               #
        !#AS DEFINED IN:                                                                 #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D: Applied     #
        !#  Physics, 1978, 11.10: 1369.                                                  #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E     (DBLE)   ENERGY                              [keV]                #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: RS

        REAL(WP) :: TMP1
        REAL(WP) :: TMP2
        REAL(DP) :: UZ
        REAL(WP) :: ETA
        REAL(WP) :: RSM

        REAL(WP), DIMENSION(3) :: CON1
        REAL(WP), DIMENSION(3) :: CON2

        DATA CON1/0.49269, 1.09870, 0.78557/
        DATA CON2/0.70256, 1.09865, 1.00460/

        ETA = BACKSCAT_COEF(Z)
        RSM = MAX_PATH_LENGTH(Z)
        UZ = VTUBE/E
        TMP1 = (CON1(1)&
                -CON1(2)*ETA&
                +CON1(3)*(ETA**2))&
                *LOG(UZ)

        TMP2 = (CON2(1)&
                -CON2(2)*ETA&
                +CON2(3)*(ETA**2))&
                +LOG(UZ)

        RS = RSM*(TMP1/TMP2)
        RETURN
    END FUNCTION MEAN_MASS_DEPTH

    FUNCTION ABSORB_TERM(E, Z) RESULT(F)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE ABSORPTION CORRECTION TERM                        #
        !#AS DEFINED IN:                                                                 #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D: Applied     #
        !#  Physics, 1978, 11.10: 1369.                                                  #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E     (DBLE)   ENERGY                              [keV]                #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: F

        REAL(WP) :: CHI
        REAL(WP) :: RS

        CHI = MAC(Z, E)/SIN(A_TAKE_OFF)
        RS = 2*MEAN_MASS_DEPTH(E, Z)

        F = (1-EXP(-RS*CHI))/(RS*CHI)
        RETURN
    END FUNCTION ABSORB_TERM

    FUNCTION CS_SMITH_REED(E, Z, N) RESULT(CS)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE SMITH AND REED CROSS SECTION                      #
        !#AS DEFINED IN:                                                                 #
        !#  EBEL, Horst. X‐ray tube spectra. X‐Ray Spectrometry, 1999, 28.4: 255-266.    #                                                 #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E     (DBLE)   ENERGY                              [keV]                #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#      -N     (INT)    LINE SELECTOR                       [-]                  #
        !#THE LINE SELECTOR IS SET TO ZERO TU USE THE CROSS SECTION FOR THE              #
        !#CONTINUUM.                                                                     #
        !#################################################################################
        INTEGER, INTENT(IN) :: Z
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: CS

        REAL(WP), DIMENSION(3) :: CON

        REAL(WP) :: X
        REAL(WP) :: CONST

        DATA CON/1.109, 0.00435, 0.00175/

        CONST = EBEL_CONST(N, Z)
        X = CON(1)&
                -CON(2)*Z&
                +CON(3)*VTUBE
        CS = CONST*Z*((VTUBE/E)-1)**X
        RETURN
    END FUNCTION CS_SMITH_REED

    FUNCTION EBEL_CONST(N, Z)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE CONSTANT USED INT SMITH AND REED CROSS SECTION    #
        !#AS DEFINED IN:                                                                 #
        !#  EBEL, Horst. X‐ray tube spectra. X‐Ray Spectrometry, 1999, 28.4: 255-266.    #
        !#                                                                               #
        !#THE CONSTANTS THEMSELVES ARE EXPERIMENTALLY DETERMINED IN:                     #
        !#  EBEL, Horst. Lι, Lα1, 2, Lη, Lβ1, 2, 3, 4, 5, 6 and Lγ1, 2, 3 spectra        #
        !3  of x‐ray tubes. X‐Ray Spectrometry, 2003, 32.1: 46-51.                       #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -N     (INT)    LINE SELECTOR                       [-]                  #
        !#      -Z     (INT)    ATOMIC NUMBER                       [-]                  #
        !#THE LINE SELECTOR IS SET TO ZERO TU USE THE CROSS SECTION FOR THE              #
        !#CONTINUUM.                                                                     #
        !#################################################################################
        REAL(QP)    :: EBEL_CONST
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: Z

        REAL(QP)    :: FCOR = 0_QP
        REAL(DP)    :: UZ = 0_DP
        REAL(DP)    :: E_EDGE = 0_DP
        REAL(DP)    :: E_TUBE = 0_DP
        REAL(QP)    :: TMP1 = 0_QP
        REAL(QP)    :: TMP2 = 0_QP

        REAL(QP)    :: L1SHELL
        REAL(QP)    :: L2SHELL
        REAL(QP)    :: L3SHELL
        REAL(QP)    :: M1SHELL
        REAL(QP)    :: M2SHELL
        REAL(QP)    :: M3SHELL
        REAL(QP)    :: M4SHELL
        REAL(QP)    :: M5SHELL
        REAL(WP) :: CONT

        REAL(QP), DIMENSION(4)    :: KSHELL
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
        DATA CONT/1.35E9/
        DATA FCON/0.4814, 0.03781, 2.413E-4/

        E_TUBE = DBLE(VTUBE)
        E_EDGE = EdgeEnergy(Z, SHELL(N))
        UZ = E_TUBE/E_EDGE

        IF (N.EQ.0) THEN
            EBEL_CONST = CONT
            RETURN
        ENDIF

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

    FUNCTION BACKSCAT_FACTOR_INTERCEPT(UZ) RESULT(I)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTERCEPT OF THE BACKSCATTERING FACTOR            #
        !#AS DESCRIBED IN:                                                               #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D:             #
        !#  Applied Physics, 1978, 11.10: 1369.                                          #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -UZ     (DBLE)      OVERVOLTAGE RATIO               [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: UZ
        REAL(WP) :: I

        REAL(WP), DIMENSION(4) :: CON
        REAL(WP) :: LN_UZ

        DATA CON/0.33148, 0.05596, 0.06339, 0.00947/

        LN_UZ = LOG(UZ)

        I = CON(1)*LN_UZ&
                +CON(2)*(LN_UZ**2)&
                -CON(3)*(LN_UZ**3)&
                +CON(4)*(LN_UZ**4)
        RETURN
    END FUNCTION BACKSCAT_FACTOR_INTERCEPT

    FUNCTION BACKSCAT_FACTOR_GRAD(UZ) RESULT(G)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE GRADIENT OF THE BACKSCATTERING FACTOR             #
        !#AS DESCRIBED IN:                                                               #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D:             #
        !#  Applied Physics, 1978, 11.10: 1369.                                          #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -UZ     (DBLE)      OVERVOLTAGE RATIO               [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) ::UZ
        REAL(WP) :: G

        REAL(WP), DIMENSION(4) :: CON
        REAL(WP) :: LN_UZ

        DATA CON/2.87898, 1.51307, 0.81312, 0.08241/

        LN_UZ = LOG(UZ)

        G = (1/UZ)&
                *(CON(1)*LN_UZ&
                    -CON(2)*(LN_UZ**2)&
                    +CON(3)*(LN_UZ**3)&
                    -CON(4)*(LN_UZ**4))
        RETURN
    END FUNCTION BACKSCAT_FACTOR_GRAD

    FUNCTION BACKSCAT_FACTOR(E, Z) RESULT(R)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE BACKSCATTERING FACTOR AS DESCRIBED IN:            #
        !#  LOVE, G.; SCOTT, V. D. Evaluation of a new correction procedure for          #
        !#  quantitative electron probe microanalysis. Journal of Physics D:             #
        !#  Applied Physics, 1978, 11.10: 1369.                                          #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER OF WINDOW         [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: R

        REAL(WP) :: I
        REAL(WP) :: G
        REAL(WP) :: ETA
        REAL(DP) :: UZ

        UZ = VTUBE/E
        I = BACKSCAT_FACTOR_INTERCEPT(UZ)
        G = BACKSCAT_FACTOR_GRAD(UZ)
        ETA = BACKSCAT_COEF(Z)

        R = 1 - ETA*((I+ETA*G)**1.67)
        RETURN
    END FUNCTION BACKSCAT_FACTOR

    FUNCTION CS_ION(E, Z, N) RESULT(CS)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE IONIZATION CROSS SECTION AS DESCRIBED IN:         #
        !#  GREEN, Martin; COSSLETT, V. E. The efficiency of production of               #
        !#  characteristic X-radiation in thick targets of a pure element.               #
        !#  Proceedings of the Physical Society, 1961, 78.6: 1206.                       #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER OF WINDOW         [-]                  #
        !#      -N      (INT)       SHELL SELECTOR                  [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: CS

        REAL(DP) :: E_EDGE
        REAL(DP) :: UZ
        REAL(QP) :: PI
        REAL(WP) :: CON0
        REAL(WP) :: CON1
        REAL(WP) :: CON2
        REAL(WP) :: CON3
        REAL(WP), DIMENSION(2) :: CON4

        DATA CON4/1.65, 2.35/

        PI = 2.D0*DASIN(1.D0)
        E_EDGE = EDGE_ENERGY(Z, N)
        UZ = E/E_EDGE
        CON0 = 4
        CON1 = 2*PI*EXP(CON0)
        CON2 = 0.35
        CON3 = E_EDGE*(CON4(1)+CON4(2)*EXP(1-UZ))
        CS = (CON1/(E*E_EDGE))*CON2*LOG((4*E)/CON3)
        RETURN
    END FUNCTION CS_ION

    FUNCTION STOPPING_FACTOR(E, Z) RESULT(F)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE STOPPING FACTOR AS DESCRIBED IN                   #
        !#  LOVE, G.; COX, M. G.; SCOTT, V. D. A versatile atomic number correction for  #
        !#  electron-probe microanalysis. Journal of Physics D: Applied Physics, 1978,   #
        !#  11.1: 7.                                                                     #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER OF WINDOW         [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        REAL(WP) :: F

        REAL(WP), DIMENSION(2) :: CON
        REAL(WP) :: IONPOT
        REAL(WP) :: ZA_RATIO
        REAL(WP) :: V
        REAL(WP) :: F_V

        DATA CON/1.18E-5, 1.47E-6/
        ZA_RATIO = NIST_ZA_RATIO(Z)
        IONPOT = NIST_ION_POT(Z)
        V = E/IONPOT

        F_V = CON(1)*SQRT(V)+CON(2)*V

        !F = -(1/IONPOT)*ZA_RATIO*(1/F_V)
        F = -78500*(ZA_RATIO/E)*LOG((1.166*E)/IONPOT)
        RETURN
    END FUNCTION STOPPING_FACTOR

    FUNCTION DERIV_SPF(E, Z, N) RESULT(S)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE DERIVATIVE OF THESTOPPING POWER FACTOR            #
        !#AS DESCRIBED IN                                                                #
        !#  LOVE, G.; COX, M. G.; SCOTT, V. D. A versatile atomic number correction for  #
        !#  electron-probe microanalysis. Journal of Physics D: Applied Physics, 1978,   #
        !#  11.1: 7.                                                                     #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER OF WINDOW         [-]                  #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: S

        REAL(WP) :: CS
        REAL(WP) :: SF

        CS = CS_ION(E, Z, N)
        SF = STOPPING_FACTOR(E, Z)

        S = CS/SF
        RETURN
    END FUNCTION DERIV_SPF

    FUNCTION STOPPING_POWER_FACTOR(E, Z, N) RESULT(S)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE STOPPING POWER FACTOR AS DESCRIBED IN             #
        !#  LOVE, G.; COX, M. G.; SCOTT, V. D. A versatile atomic number correction for  #
        !#  electron-probe microanalysis. Journal of Physics D: Applied Physics, 1978,   #
        !#  11.1: 7.                                                                     #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER OF WINDOW         [-]                  #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: S

        REAL(DP) :: E_INIT
        REAL(DP) :: E_FINAL

        E_INIT = VTUBE
        E_FINAL = E

        S = INTEGRATE(DERIV_SPF, E_INIT, E_FINAL, INT(10), Z, N)
        RETURN
    END FUNCTION STOPPING_POWER_FACTOR

    FUNCTION TUBE_ATTEN(E, Z ,D) RESULT(W)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE ATTENUATION OF THE INTENSITY CAUSED BY            #
        !#THE TUBE-WINDOW                                                                #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER OF WINDOW         [-]                  #
        !#      -D      (WP)        THICKNESS OF WINDOW             [um]                 #
        !#################################################################################
        IMPLICIT NONE
        REAL(DP), INTENT(IN)    :: E
        INTEGER, INTENT(IN)     :: Z
        REAL(WP), INTENT(IN)    :: D
        REAL(WP)    :: W

        W = EXP(-MAC(Z, E)*D*NIST_DENSITY(Z)*1E-4)
        RETURN
    END FUNCTION TUBE_ATTEN

    FUNCTION FILTER_ATTEN(E, Z ,D) RESULT(W)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE ATTENUATION OF THE INTENSITY CAUSED BY            #
        !#THE FILTER                                                                     #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#      -Z      (INT)       ATOMIC NUMBER OF FILTER         [-]                  #
        !#      -D      (WP)        THICKNESS OF FILTER             [um]                 #
        !#################################################################################
        IMPLICIT NONE
        REAL(DP), INTENT(IN)    :: E
        INTEGER, INTENT(IN)     :: Z
        REAL(WP), INTENT(IN)    :: D
        REAL(WP)    :: W

        IF (D.GT. 0) THEN
            W = EXP(-MAC(Z, E)*D*NIST_DENSITY(Z)*1E-4)
            RETURN
        ELSE
            W = 1_QP
            RETURN
        ENDIF
        RETURN
    END FUNCTION FILTER_ATTEN
END MODULE ANODE
