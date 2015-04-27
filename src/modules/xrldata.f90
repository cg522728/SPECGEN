MODULE XRLDATA
    USE :: xraylib
    USE :: TYPES
    USE :: CFGDATA
    implicit none
    INTEGER, DIMENSION(111)   :: SHELL
    INTEGER, DIMENSION(111)   :: TOSHELL
    INTEGER, DIMENSION(111)   :: LINE
    REAL(QP), DIMENSION(92)   :: NIST_ZA_RATIO
    REAL(WP), DIMENSION(92)   :: NIST_ION_POT
    REAL(QP), DIMENSION(92)   :: NIST_DENSITY
    INCLUDE 'shell.f90'
    INCLUDE 'line.f90'
CONTAINS
    FUNCTION CS_FLUOR_CHG(Z, N, E) RESULT(CS)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE X-RAY FLUORESCENCE PRODUCTION CROSS SECTION       #
        !#FOR AN ELEMENT                                                                 #
        !#################################################################################
        !#INPUT:                                                                         #
        !#      -Z      (INT)       ATOMIC NUMBER                                        #
        !#      -N      (INT)       LINE SELECTOR                                        #
        !#      -E      (DBLE)      ENERGY                                               #
        !#################################################################################
        IMPLICIT NONE
        INTEGER, INTENT(IN)     :: Z
        INTEGER, INTENT(IN)     :: N
        REAL(DP), INTENT(IN)    :: E
        REAL(WP)    :: CS

        REAL(DP)    :: E_EDGE

        E_EDGE = EDGE_ENERGY(Z, N)
        IF (E.LT.E_EDGE) THEN
            CS = 0_WP
            RETURN
        ENDIF
        CALL SetErrorMessages(0)
        !CS = CS_Fluorline_Kissel_Cascade(Z, LINE(N), E)
        CS = TRANSPROB(Z,N)*FLUOR_YIELD(Z,N)*CS_PHOTO_CHG(Z, E)
        CALL SetErrorMessages(1)
        RETURN
    END FUNCTION CS_FLUOR_CHG

    FUNCTION MAC(Z, E) RESULT(MU)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE MASS ABSORPTION COEFFICIENT FOR AN ELEMENT        #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                                            #
        !#      -E      (DBLE)  ENERGY                                                   #
        !#################################################################################
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Z
        REAL(DP), INTENT(IN) :: E
        REAL(WP) :: MU

        IF (E.LE.0) THEN
            MU = 0_WP
            RETURN
        ENDIF
        MU= CS_TOTAL(Z, E)
        RETURN
    END FUNCTION

    FUNCTION MAC_COMP(CP, E) RESULT(MU)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE MASS ABSORPTION COEFFICIENT FOR A COMPOUND        #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                                            #
        !#      -E      (DBLE)  ENERGY                                                   #
        !#################################################################################
        IMPLICIT NONE
        TYPE(CompoundData), POINTER, INTENT(IN) :: CP
        REAL(DP), INTENT(IN) :: E
        REAL(WP) :: MU

        INTEGER :: Z = 0
        REAL(WP) :: MAC_TMP = 0_WP

        INTEGER :: CNT

        DO CNT=1,CP%NELEMENTS
            Z = CP%ELEMENTS(CNT)
            MAC_TMP = MAC_TMP + MAC(Z, E)*CP%MASSFRACTIONS(CNT)
        END DO
        MU = MAC_TMP
        MAC_TMP = 0_WP
        RETURN
    END FUNCTION

    FUNCTION CS_FLUOR_CP_CHG(CP, N, E) RESULT(CS)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE X-RAY FLUORESCENCE PRODUCTION CROSS SECTION       #
        !#FOR A COMPOUND                                                                 #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -CP     (STRUC) XRAYLIB COMPOUND DATA                                    #
        !#      -N      SHELL SELECTOR                                                   #
        !#      -E      (DBLE)  ENERGY                                                   #
        !#################################################################################
        IMPLICIT NONE
        TYPE(CompoundData), POINTER, INTENT(IN) :: CP
        INTEGER, INTENT(IN) :: N
        REAL(DP), INTENT(IN) :: E
        REAL(WP) :: CS

        REAL(WP) :: CS_TMP = 0_WP

        INTEGER :: CNT

        DO CNT=1,CP%NELEMENTS
            CS_TMP = CS_TMP&
                        + CS_FLUOR_CHG(CP%ELEMENTS(CNT), N, E)&
                            *CP%MASSFRACTIONS(CNT)
        END DO
        CS = CS_TMP
        CS_TMP = 0_WP
        RETURN
    END FUNCTION

    FUNCTION TRANSPROB(Z, N)
        !#################################################################################
        !#THIS FUNCTION CALCULATES TRANSITION PROBABILITY AS DEFINED IN                  #
        !#  X-Ray Fluorescence Cross Sections for K and L X-Rays of the Elements         #
        !3  M.O. Krause, C.W. Nestor, Jr., C.J. Sparks, Jr. and E. Ricci                 #
        !#  ORNL-5399, pag 120, 1978                                                     #
        !#THIS IS DONE BY MULTIPLYING THE RELATIVE RADIATIVE RATE WITH                   #
        !#THE ATOMIC LEVEL WIDTH                                                         #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                                            #
        !#      -N      (INT)   LINE SELECTOR                                            #
        !#################################################################################
        IMPLICIT NONE
        REAL(QP)    :: TRANSPROB
        INTEGER, INTENT(IN)     :: N
        INTEGER, INTENT(IN)     :: Z
        REAL(WP) :: LW

        IF (Z.EQ.Z_ANODE) THEN
            LW = AtomicLevelWidth(Z, SHELL(N))
        ELSE
            LW = 1
        ENDIF
        TRANSPROB = RadRate(Z, LINE(N))*LW*1E3
        RETURN
    END FUNCTION TRANSPROB
    SUBROUTINE LOAD_NIST()
        !#################################################################################
        !#THIS SUBROUTINE LOADS THE FOLLOWING ELEMENTAL DATA FROM A FILE:                #
        !#      -Z/A-RATIO                                                               #
        !#      -MEAN IONISATION POTENTIAL IN eV                                         #
        !#      -DENSITY IN g/cm^3                                                       #
        !#SOME DENSITY VALUES ARE NOMINAL, VALUES FOR Z=85 AND Z=87 WERE SET TO 10       #
        !#TO COMPLETE THE SET.                                                           #
        !#THE DATA IS STORED IN THE FOLLOWING ARRAYS:                                    #
        !#  -Z/A-RATIO                       : NIST_ZA_RATIO                             #
        !#  -MEAN IONISATION POTENTIAL       : NIST_ION_POT                              #
        !#  -DENSITY                         : NIST_DENSITY                              #
        !#THIS DATASET WAS OBTAINED FROM:                                                #
        !#  http://http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html            #
        !#################################################################################
        IMPLICIT NONE
        INTEGER     :: Z
        CHARACTER(LEN=2)    :: ELEMENT_SYM
        CHARACTER(LEN=9)    :: ELEMENT_NAME
        REAL(QP)    :: ZA_RATIO
        REAL(WP)    :: I
        REAL(QP)    :: DENSITY

        NIST_ZA_RATIO = 0_QP
        NIST_ION_POT = 0_WP
        NIST_DENSITY = 0_QP

        OPEN(UNIT=1,FILE='nist.dat')
        DO
            READ (1,*,END=999) Z, ELEMENT_SYM, ELEMENT_NAME, ZA_RATIO, I, DENSITY
                NIST_ZA_RATIO(Z) = ZA_RATIO
                NIST_ION_POT(Z) = I*(1E-3)
                NIST_DENSITY(Z) = DENSITY
        END DO
999     RETURN
    END SUBROUTINE LOAD_NIST

    FUNCTION CS_PHOTO_CHG(Z, E) RESULT(CS)
        !#################################################################################
        !#THIS FUNCTION PROVIDES AN INTERFACE TO THE XRAYLIB FUNCTION                    #
        !#      EdgeEnergy(Z, N)                                                         #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                                            #
        !#      -N      (INT)   SHELL SELECTOR                                           #
        !#################################################################################
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Z
        REAL(DP), INTENT(IN) :: E
        REAL(WP) :: CS

        CS = CS_Photo(Z, E)
        RETURN
    END FUNCTION CS_PHOTO_CHG

    FUNCTION EDGE_ENERGY(Z, N) RESULT(E)
        !#################################################################################
        !#THIS FUNCTION PROVIDES AN INTERFACE TO THE XRAYLIB FUNCTION                    #
        !#      EdgeEnergy(Z, N)                                                         #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                                            #
        !#      -N      (INT)   SHELL SELECTOR                                           #
        !#################################################################################
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(DP) :: E

        E = EdgeEnergy(Z, SHELL(N))
        RETURN
    END FUNCTION EDGE_ENERGY

    FUNCTION FLUOR_YIELD(Z, N) RESULT(W)
        !#################################################################################
        !#THIS FUNCTION PROVIDES AN INTERFACE TO THE XRAYLIB FUNCTION                    #
        !#      FluorYield(Z, N)                                                         #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                                            #
        !#      -N      (INT)   SHELL SELECTOR                                           #
        !#################################################################################
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: W

        W = FluorYield(Z, SHELL(N))
        RETURN
    END FUNCTION FLUOR_YIELD

    FUNCTION LINE_ENERGY(Z, N) RESULT(E)
        !#################################################################################
        !#THIS FUNCTION PROVIDES AN INTERFACE TO THE XRAYLIB FUNCTION                    #
        !#      LineEnergy(Z, N)                                                         #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                                            #
        !#      -N      (INT)   LINE SELECTOR                                            #
        !#################################################################################
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: Z
        INTEGER, INTENT(IN) :: N
        REAL(DP) :: E

        E = LineEnergy(Z, LINE(N))
        RETURN
    END FUNCTION LINE_ENERGY

    FUNCTION DCS_R(E, A) RESULT(DCS)
        !#################################################################################
        !#THIS FUNCTION PROVIDES AN INTERFACE TO THE XRAYLIB FUNCTION                    #
        !#      DCS_Rayls(Z, E, A)                                                       #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                       [-]                  #
        !#      -E      (DBLE)  ENERGY                              [keV]                #
        !#      -A      (WP)    SCATTERING POLAR ANGLE              [rad]                #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        REAL(WP), INTENT(IN) :: A
        REAL(WP) :: DCS

        REAL(WP) :: DCS_TMP = 0_WP

        INTEGER :: CNT

        DO CNT = 1, CP_ST%NELEMENTS
            DCS_TMP = DCS_TMP + (DCS_Rayl(CP_ST%ELEMENTS(CNT), E, DBLE(A))&
                        *CP_ST%MASSFRACTIONS(CNT))
        END DO
        DCS = DCS_TMP
        DCS_TMP = 0_WP
        RETURN
    END FUNCTION DCS_R

    FUNCTION DCS_C(E, A) RESULT(DCS)
        !#################################################################################
        !#THIS FUNCTION PROVIDES AN INTERFACE TO THE XRAYLIB FUNCTION                    #
        !#      DCS_Compt(Z, E, A)                                                       #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z      (INT)   ATOMIC NUMBER                       [-]                  #
        !#      -E      (DBLE)  ENERGY                              [keV]                #
        !#      -A      (WP)    SCATTERING POLAR ANGLE              [rad]                #
        !#################################################################################
        REAL(DP), INTENT(IN) :: E
        REAL(WP), INTENT(IN) :: A
        REAL(WP) :: DCS

        REAL(WP) :: DCS_TMP = 0_WP

        INTEGER :: CNT

        DO CNT = 1, CP_ST%NELEMENTS
            DCS_TMP = DCS_TMP + (DCS_Compt(CP_ST%ELEMENTS(CNT), E, DBLE(A))&
                        *CP_ST%MASSFRACTIONS(CNT))
        END DO
        DCS = DCS_TMP
        DCS_TMP = 0_WP
        RETURN
    END FUNCTION DCS_C

    FUNCTION COMPTON_ENERGY(E, A) RESULT(EI)
        REAL(DP), INTENT(IN) :: E
        REAL(WP), INTENT(IN) :: A
        REAL(DP) :: EI

        EI = ComptonEnergy(E, DBLE(A))
        RETURN
    END FUNCTION COMPTON_ENERGY
END MODULE XRLDATA
