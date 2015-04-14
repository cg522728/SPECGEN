MODULE PELLA
        USE :: xraylib
        USE :: XRLDATA
        USE :: CFGDATA
        USE :: CONSTANTS
CONTAINS
    FUNCTION PELLA_ANODECONT(EI)
    !#################################################################################
    !#THIS FUNCTION CALCULATES THE INTENSITY OF THE CONTINUUM OF THE GIVEN           #
    !#X-RAY TUBE AT THE GIVEN ENERGY. THE INPUTS ARE:                                #
    !#      -EI     (DBLE)  :ENERGY IN keV                                           #
    !#THE CALCULATION ALGORITHM IS BASED ON THE PELLA ALGORITHM (PELLA ET AL, 1985)  #
    !#FROM 'X-RAY SPECTROMETRY'.                                                     #
    !#################################################################################
        IMPLICIT NONE

        REAL(16) :: PELLA_ANODECONT
        REAL(16) :: CONST
        REAL(16) :: EI, UZ
        REAL(16) :: X
        REAL(16) :: F
        REAL(16) :: I

        DATA    CONST/1.37E9/

        X = 1.0314-0.00032*Z_ANODE+0.0047*VTUBE
        UZ = VTUBE/EI
        I = SA_ANODE_OUT&
                *ITUBE&
                *CONST&
                *Z_ANODE&
                *(((UZ)-1)**X)&
                *PELLA_CALC_F(EI)
        I = I*EXP(-MAC(Z_WINDOW, EI)&
                    *D_WINDOW&
                    *ElementDensity(Z_WINDOW)&
                    *1E-4)
        IF (D_FILTER.GT. 0) THEN
            I = I*EXP(-MAC(Z_FILTER, EI)&
                        *D_FILTER&
                        *ElementDensity(Z_FILTER)&
                        *1E-4)
        ENDIF
        PELLA_ANODECONT = I*ESTEP
        RETURN
    END FUNCTION PELLA_ANODECONT
    FUNCTION PELLA_ANODECHAR(N, KAPPA)
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
        IMPLICIT NONE
        REAL(16) :: PELLA_ANODECHAR
        INTEGER :: N
        REAL(16) :: EI, EC, UZ
        REAL(16) :: TAU
        REAL(16) :: R
        REAL(16) :: KAPPA
        REAL(16) :: PI, K, CONST

        PI = 2.D0*DASIN(1.D0)
        EC = LineEnergy(Z_ANODE, LINE(N))
        UZ = VTUBE/EC
        IF (EC.EQ.0) THEN
            PELLA_ANODECHAR = 0._16
            RETURN
        ENDIF
        TAU = ((UZ*LOG(UZ/(UZ-1)))-1)
        R = CALC_KR(Z_ANODE, EC, N, KAPPA)*TAU
        IF (R.EQ. 0) THEN
            PELLA_ANODECHAR = PELLA_ANODECONT(EC)
            RETURN
        ENDIF
        PELLA_ANODECHAR = PELLA_ANODECONT(EC)&
                    *R&
                    *(EC**2/12.396)
        IF (PELLA_ANODECHAR.EQ.0) PELLA_ANODECHAR = PELLA_ANODECONT(EC)
        RETURN
    END FUNCTION PELLA_ANODECHAR
END MODULE PELLA
