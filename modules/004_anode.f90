MODULE ANODE
    USE :: xraylib
    USE :: CFGDATA
    USE :: XRLDATA
    IMPLICIT NONE

CONTAINS
    FUNCTION ANODECONT(EI)
    !#################################################################################
    !#THIS FUNCTION CALCULATES THE INTENSITY OF THE CONTINUUM OF THE GIVEN           #
    !#X-RAY TUBE AT THE GIVEN ENERGY. THE INPUTS ARE:                                #
    !#      -EI     (DBLE)  :ENERGY IN keV                                           #
    !#THE CALCULATION ALGORITHM IS BASED ON THE PELLA ALGORITHM (PELLA ET AL, 1985)  #
    !#FROM 'X-RAY SPECTROMETRY'.                                                     #
    !#################################################################################
        USE :: xraylib
        USE :: XRLDATA
        USE :: CFGDATA
        IMPLICIT NONE

        REAL(16) :: ANODECONT
        REAL(16) :: CONST
        REAL(16) :: EI
        REAL(16) :: X
        REAL(16) :: F
        REAL(16) :: I
        REAL(16) :: ETA, ETACONST
        REAL(16) :: C, C1, C2, C3

        DATA    CONST/1.35E9/

        X = 1.109-0.00435*Z_ANODE+0.00175*VTUBE
        ETA = VTUBE**(1.65)-EI**(1.65)
        ETACONST = MAC(Z_ANODE, EI)/(SIN(A_TAKE_OFF)*63.67)
        ETA = ETA*ETACONST
        C1 = 1+2.56E-3*Z_ANODE**2
        C2 = 1+3.17E4*VTUBE**(-1)*Z_ANODE**(-2)
        C3 = 0.25*ETA+1E4
        C = (1+C1**(-1))*C2**(-1)*C3**(-1)
        F = (1+C*ETA)**(-2)
        I = SA_ANODE_OUT&
                *ITUBE&
                *CONST&
                *Z_ANODE&
                *(((VTUBE/EI)-1)**X)&
                *F
        I = I*EXP(-MAC(Z_WINDOW, EI)*D_WINDOW*1E-4)   !Be-venster
        IF (D_WINDOW.GT. 0) THEN
            I = I*EXP(-MAC(Z_FILTER, EI)*D_FILTER*1E-4)
        ENDIF
        ANODECONT = I*ESTEP
        RETURN
    END FUNCTION ANODECONT
    FUNCTION ANODECHAR(N, KAPPA)
    !#################################################################################
    !#THIS FUNCTION CALCULATES THE INTENSITY OF THE CHARACETRISTIC LINES OF          #
    !#OF THE GIVEN X-RAY TUBE AT THE GIVEN LINE. THE INPUTS ARE:                     #
    !#      -N     (INT)  :IDENTIFIER OF THE CHARACTERISTIC LINE AS DEFINED IN       #
    !#                          XRAYLIB                                              #
    !#THE CALCULATION ALGORITHM IS BASED ON THE PELLA ALGORITHM (PELLA ET AL, 1985)  #
    !#FROM 'X-RAY SPECTROMETRY'. HOWEVER, THE ALGORITHM HAS BEEN ADJUSTED TO WORK    #
    !#FOR HIGHER ATOMIC NUMBERS THAN THE ORIGINAL ALGORITHM. TO DO THIS WE HAVE      #
    !#REMOVED A SIMPLIFICATION WHERE A FIT WAS USED TO SIMULATE A CONSTANT COMPOSED  #
    !#OF A NUMBER OF PHYSICAL CONSTANTS. THIS HAS BEEN IMPLEMENTED                   #
    !#IN THE FUNCTION CALCKR(). A PROPORTIONALITY FACTOR KAPPA REMAINS.              #
    !#THIS KAPPA HAS AN INFLUENCE ON THE RATIO BETWEEN THE CHARACTERISTIC LINE       #
    !#INTENSITY AND THE INTENSITY OF THE BACKGROUND.                                 #
    !#################################################################################
        USE :: xraylib
        USE :: XRLDATA
        USE :: CFGDATA
        IMPLICIT NONE
        REAL(16) :: ANODECHAR
        INTEGER :: N
        REAL(16) :: EI, EC
        REAL(16) :: TAU
        REAL(16) :: R
        REAL(16) :: KAPPA
        REAL(16) :: PI, K

        PI = 2.D0*DASIN(1.D0)
        EC = LineEnergy(Z_ANODE, LINE(N))
        IF (EC.EQ.0) THEN
            ANODECHAR = 0._16
            RETURN
        ENDIF
        TAU = (((VTUBE/EC)*LOG((VTUBE/EC)/((VTUBE/EC)-1)))-1)
        R = CALC_KR(Z_ANODE, EC, N, KAPPA)*TAU
        IF (R.EQ. 0) THEN
            ANODECHAR = ANODECONT(EC)
            RETURN
        ENDIF
        ANODECHAR = R*ANODECONT(EC)*((EC**2)/KEV2ANGST)/(4*PI)
        RETURN
    END FUNCTION ANODECHAR
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
        USE :: xraylib
        USE :: XRLDATA
        USE :: CFGDATA
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
        OMEGAQ = FluorYield(Z, SHELL(SHELLS))
        PI = RadRate(Z, LINE(SHELLS))
        CALC_KR = (KAPPA/K)*PI*(1+P)*OMEGAQ*R*(LAMBDAL**2)*((A*C*Z)**(-1))
        RETURN
    END FUNCTION
END MODULE ANODE
