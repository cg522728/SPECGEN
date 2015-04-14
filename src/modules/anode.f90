MODULE ANODE
    USE :: xraylib
    USE :: CFGDATA
    USE :: XRLDATA
    USE :: CONSTANTS
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
        REAL(16) :: CONST, EBEL
        REAL(16) :: EI, E1, E2
        REAL(16) :: X
        REAL(16) :: F
        REAL(16) :: I, TMP1, TMP2
        REAL(16) :: ETA, ETACONST
        REAL(16) :: C, C1, C2, C3

        DATA    EBEL/1.37E9/

        E1 = EI - (ESTEP/2)
        E2 = EI + (ESTEP/2)
        CONST = SA_ANODE_OUT*ITUBE*EBEL*Z_ANODE
        X = 1.109-0.00435*Z_ANODE+0.00175*VTUBE
        TMP1 = CONST*(((VTUBE/E1)-1)**X)*CALC_F(E1)
        TMP2 = CONST*(((VTUBE/(E2))-1)**X)*CALC_F(E2)
        I = ESTEP*(TMP1+TMP2)/2
        I = I*EXP(-MAC(Z_WINDOW, EI)*D_WINDOW*ElementDensity(Z_WINDOW)*1E-4)   !Be-venster
        IF (D_FILTER.GT. 0) THEN
            I = I*EXP(-MAC(Z_FILTER, EI)*D_FILTER*ElementDensity(Z_FILTER)*1E-4)
        ENDIF
        ANODECONT = I
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
        REAL(16) :: PI, K, CONST, UZ, EBEL
        REAL(16), DIMENSION(5)   :: RCON

        DATA EBEL/6E13/
        DATA RCON/1, 0.008157, 3.613E-5, 0.009583, 0.001141/

        PI = 2.D0*DASIN(1.D0)
        EC = LineEnergy(Z_ANODE, LINE(N))
        UZ = VTUBE/EC
        R = RCON(1) - RCON(2)*Z_ANODE + RCON(3)*(Z_ANODE**2)+RCON(4)*Z_ANODE*EXP(-UZ)+RCON(5)*VTUBE
        IF (EC.EQ.0 .OR. EC.LT.EMIN) THEN
            ANODECHAR = 0._16
            RETURN
        ENDIF
        ANODECHAR = SA_ANODE_OUT&
                    *ITUBE&
                    *EBEL&
                    *R&
                    *STOPPINGFACTOR(EC, N)&
                    *CALC_F(EC)&
                    *RadRate(Z_ANODE, LINE(N))&
                    *FLUORYIELD_CHG(Z_ANODE, SHELL(N))
        IF (ISNAN(ANODECHAR)) ANODECHAR = ANODECONT(EC)
        RETURN
    END FUNCTION ANODECHAR
    FUNCTION CALC_F(EI)
        USE :: xraylib
        REAL(16)    :: CALC_F
        REAL(16)    :: EI
        REAL(16)    :: TMP, TMP1
        TMP = CS_Photo(Z_ANODE, DBLE(EI))*2*DDF(EI)*(SIN(A_INCID)/SIN(A_TAKE_OFF))
        TMP1 = 1-EXP(-TMP)
        CALC_F = TMP1/TMP
        RETURN
    END FUNCTION CALC_F
        FUNCTION DDF(EI)
        USE :: xraylib
        REAL(16)    ::DDF
        REAL(16)    :: EI
        REAL(16)    :: ETA
        REAL(16)    :: TMP, M, UZ, J, TMP2, TMP3
        UZ = VTUBE/EI
        M = 0.1382-(0.9211/SQRT(DBLE(Z_ANODE)))
        J = 0.0135*Z_ANODE
        ETA = (VTUBE**M)*(0.1904-0.2236*LOG(DBLE(Z_ANODE))+0.1292*(LOG(DBLE(Z_ANODE))**2)&
                -0.0149*(LOG(DBLE(Z_ANODE))**3))
        TMP = (AtomicWeight(Z_ANODE)/DBLE(Z_ANODE))&
            *(0.787E-5*SQRT(J*VTUBE**3)+0.735E-6*VTUBE**2)
        TMP2 = 0.49269-1.0987*ETA+0.78557*ETA**2
        TMP3 = 0.70256-1.09865*ETA+1.0046*ETA**2+LOG(UZ)
        DDF = TMP*(TMP2/TMP3)*LOG(UZ)
        RETURN
    END FUNCTION DDF
END MODULE ANODE
