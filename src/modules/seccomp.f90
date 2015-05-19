MODULE SECCOMP
    USE :: xraylib
    USE :: TYPES
    USE :: XRLDATA
    USE :: CFGDATA
    USE :: ANODE
    USE :: MATHCHG
    implicit none
    !PRIVATE
    !PUBLIC ST_CHAR , ST_SCAT_R, ST_SCAT_C, ST_SCAT_CONT, DERIV_ST_SCAT
CONTAINS
    FUNCTION ST_CHAR(Z_ST, N) RESULT(N_CHAR)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTENSITY OF THE CHARACTERISTIC LINE OF           #
        !#OF THE SECONDARY TARGET AS DECRIBED IN :                                       #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z_ST   (INT)       ATOMIC NUMBER                   [-]                  #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        INTEGER, INTENT(IN) :: Z_ST
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: N_CHAR

        REAL(WP) :: I_X_AN_CHAR
        REAL(WP) :: I_X_AN_CONT
        REAL(DP) :: E_ST_EDGE
        REAL(DP) :: E_ST_LINE
        REAL(WP) :: CONST
        REAL(WP) :: PI

        PI = 2.D0*DASIN(1.D0)
        E_ST_EDGE = EDGE_ENERGY(Z_ST, N)
        E_ST_LINE = LINE_ENERGY(Z_ST, N)

        IF (E_ST_EDGE.EQ.0 .OR. E_ST_LINE.EQ.0) THEN
            N_CHAR = 0_WP
            RETURN
        ENDIF

        I_X_AN_CHAR = CHAR_X_AN_CHAR(Z_ST, N)
        I_X_AN_CONT = CHAR_X_AN_CONT(Z_ST, N)
        CONST = SA_ST_IN/(4*PI*SIN(A_ST_INCID))
        N_CHAR = CONST*(I_X_AN_CHAR + I_X_AN_CONT)*(ATOMICWEIGHT(Z_ST)/ST_MASS)
        RETURN
    END FUNCTION ST_CHAR

    FUNCTION CHAR_X_AN_CHAR(Z_ST, N) RESULT(I)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTENSITY OF THE CHARACTERISTIC LINE OF           #
        !#OF THE SECONDARY TARGET WHICH IS EXCITED BY THE TUBE CONTINUUM                 #
        !#AS DECRIBED IN :                                                               #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z_ST   (INT)       ATOMIC NUMBER                   [-]                  #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        INTEGER, INTENT(IN) :: Z_ST
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: I

        INTEGER :: Z_AN
        REAL(DP) :: E_ST_EDGE
        REAL(DP) :: E_ST_LINE
        REAL(DP) :: E_AN_LINE
        REAL(WP) :: TMP
        REAL(WP) :: CHI
        REAL(WP) :: CS
        REAL(WP) :: I_TMP
        REAL(WP) :: I_SUM

        INTEGER :: CNT

        Z_AN = Z_ANODE
        E_ST_EDGE = EDGE_ENERGY(Z_ST, N)
        E_ST_LINE = LINE_ENERGY(Z_ST, N)

        DO CNT = 1, SIZE(LINE)
            E_AN_LINE = LINE_ENERGY(Z_AN, CNT)
            IF (E_AN_LINE.EQ.0) CYCLE
            IF (E_AN_LINE.LT.E_ST_EDGE) CYCLE
            IF (E_AN_LINE.LT.EMIN) CYCLE
            I_TMP = ANODE_CHAR(Z_AN, CNT)
            CS = CS_FLUOR_CP_CHG(CP_ST, N, E_AN_LINE)
            CHI = CALC_CHI(E_AN_LINE, A_ST_INCID, E_ST_LINE, A_ST_TAKE_OFF)
            TMP = (I_TMP*CS)/CHI
            I_SUM = I_SUM + TMP
        END DO
        I = I_SUM
        RETURN
    END FUNCTION CHAR_X_AN_CHAR

    FUNCTION CHAR_X_AN_CONT(Z_ST, N) RESULT(I)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTENSITY OF THE CHARACTERISTIC LINE OF           #
        !#OF THE SECONDARY TARGET WHICH IS EXCITED BY THE CHARACTERISIC LINES OF         #
        !#THE TUBE AS DECRIBED IN :                                                      #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -Z_ST   (INT)       ATOMIC NUMBER                   [-]                  #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        INTEGER, INTENT(IN) :: Z_ST
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: I

        REAL(DP) :: E1
        REAL(DP) :: E2
        REAL(DP) :: E_ST_EDGE
        REAL(DP) :: E_ST_LINE

        E1 = EDGE_ENERGY(Z_ST, N)
        E2 = VTUBE

        I = INTEGRATE(DERIV_CHAR_X_AN_CONT, E1, E2, NSTEP, Z_ST, N)
        RETURN
    END FUNCTION CHAR_X_AN_CONT

    FUNCTION DERIV_CHAR_X_AN_CONT(E, Z_ST, N) RESULT(I)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE DERIVATIVE OF THE INTENSITY OF                    #
        !#THE CHARACTERISTIC LINE OF OF THE SECONDARY TARGET WHICH IS EXCITED            #
        !#BY THE CONTINUUM OF THE TUBE AS DECRIBED IN :                                  #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY IN CONTINUUM             [keV]                #
        !#      -Z_ST   (INT)       ATOMIC NUMBER                   [-]                  #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        REAL(DP), INTENT(IN) :: E
        INTEGER, INTENT(IN) :: Z_ST
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: I

        INTEGER :: Z_AN
        REAL(DP) :: E_ST_EDGE
        REAL(DP) :: E_ST_LINE
        REAL(WP) :: TMP
        REAL(WP) :: CHI

        Z_AN = Z_ANODE
        E_ST_EDGE = EDGE_ENERGY(Z_ST, N)
        E_ST_LINE = LINE_ENERGY(Z_ST, N)

        IF (E.LT.E_ST_EDGE .OR. E.GT.VTUBE) THEN
            I = 0_QP
            RETURN
        ENDIF
        TMP = ANODE_CONT(Z_AN, E)&
            *CS_FLUOR_CP_CHG( CP_ST, N, E)
        CHI = CALC_CHI( E, A_ST_INCID, E_ST_EDGE, A_ST_TAKE_OFF)
        I = (TMP/CHI)
        RETURN
    END FUNCTION DERIV_CHAR_X_AN_CONT

    FUNCTION CALC_CHI(E1, A1, E2, A2, Z_ST) RESULT(CHI)
        !#################################################################################
        !#THIS FUNCTION CALCULATES A FACTOR CHI USED IN OTHER FUNCTIONS                  #
        !#AS DECRIBED IN :                                                               #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E1     (DBLE)      ENERGY @ INCIDENCE              [keV]                #
        !#      -A1     (WP)        INCIDENCE ANGLE                 [rad]                #
        !#      -E2     (DBLE)      ENERGY @ EXIT                   [keV]                #
        !#      -A2     (WP)        TAKE OFF ANGLE                  [rad]                #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        REAL(DP), INTENT(IN) :: E1
        REAL(WP), INTENT(IN) :: A1
        REAL(DP), INTENT(IN) :: E2
        REAL(WP), INTENT(IN) :: A2
        INTEGER, INTENT(IN), OPTIONAL :: Z_ST
        REAL(WP) :: CHI

        LOGICAL :: Z_PRESENT
        REAL(WP) :: MAC1
        REAL(WP) :: MAC2

        Z_PRESENT = PRESENT(Z_ST)
        IF (Z_PRESENT) THEN
            MAC1 = MAC(Z_ST, E1)
            MAC2 = MAC(Z_ST, E2)
        ELSEIF (.NOT.Z_PRESENT) THEN
            MAC1 = MAC_COMP(CP_ST, E1)
            MAC2 = MAC_COMP(CP_ST, E2)
        ENDIF

        CHI = (MAC1/SIN(A1))&
            +(MAC2/SIN(A2))
        RETURN
    END FUNCTION CALC_CHI

    FUNCTION ST_SCAT_R(N) RESULT(N_RAYL)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTENSITY OF THE RAYLEIGH SCATTERED LINES         #
        !#OF THE TUBE AS DESCRIBED IN:                                                   #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: N_RAYL

        INTEGER :: Z_AN
        REAL(WP) :: CHI
        REAL(DP) :: E_AN_CHAR
        REAL(DP) :: E_AN_SCAT
        REAL(WP) :: I_CHAR
        REAL(WP) :: CONST
        REAL(WP) :: DCS
        REAL(WP) :: TMP

        Z_AN = Z_ANODE
        E_AN_CHAR = LINE_ENERGY(Z_AN, N)
        E_AN_SCAT = E_AN_CHAR
        CONST = SA_ST_IN/SIN(A_ST_INCID)
        I_CHAR = ANODE_CHAR(Z_AN, N)
        DCS = DCS_R(E_AN_CHAR, A_ST_POL)
        CHI = CALC_CHI(E_AN_CHAR, A_ST_INCID, E_AN_SCAT, A_ST_TAKE_OFF)

        N_RAYL = CONST*I_CHAR*DCS/CHI
        RETURN
    END FUNCTION ST_SCAT_R

    FUNCTION ST_SCAT_C(N) RESULT(N_RAYL)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTENSITY OF THE COMPTON SCATTERED LINES          #
        !#OF THE TUBE AS DESCRIBED IN:                                                   #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -N      (INT)       LINE SELECTOR                   [-]                  #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        INTEGER, INTENT(IN) :: N
        REAL(WP) :: N_RAYL

        INTEGER :: Z_AN
        REAL(WP) :: CHI
        REAL(DP) :: E_AN_CHAR
        REAL(DP) :: E_AN_SCAT
        REAL(WP) :: I_CHAR
        REAL(WP) :: CONST
        REAL(WP) :: DCS

        Z_AN = Z_ANODE
        E_AN_CHAR = LINE_ENERGY(Z_AN, N)
        E_AN_SCAT = COMPTON_ENERGY(E_AN_CHAR, A_ST_POL)
        CONST = SA_ST_IN/SIN(A_ST_INCID)
        I_CHAR = ANODE_CHAR(Z_AN, N)
        DCS = DCS_C(E_AN_CHAR, A_ST_POL)
        CHI = CALC_CHI(E_AN_CHAR, A_ST_INCID, E_AN_SCAT, A_ST_TAKE_OFF)

        N_RAYL = CONST*I_CHAR*DCS/CHI
        RETURN
    END FUNCTION ST_SCAT_C

    FUNCTION ST_SCAT_CONT(E) RESULT(N_SCAT)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE INTENSITY OF THE SCATTERED CONTINUUM              #
        !#OF THE TUBE AS DESCRIBED IN:                                                   #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        REAL(DP), INTENT(IN) :: E
        REAL(WP) :: N_SCAT

        REAL(DP) :: E1
        REAL(DP) :: E2
        REAL(WP) :: CONST

        CONST = SA_ST_IN/SIN(A_ST_INCID)
        E1 = E - (ESTEP/2)
        E2 = E + (ESTEP/2)

        IF (E1.GE.VTUBE) THEN
            N_SCAT = 0_WP
            RETURN
        ENDIF
        IF (E2.GE.VTUBE) THEN
            N_SCAT = 0_WP
            RETURN
        ENDIF
        N_SCAT = CONST*INTEGRATE(DERIV_ST_SCAT, E1, E2, INT(25))
        RETURN
    END FUNCTION ST_SCAT_CONT

    FUNCTION DERIV_ST_SCAT(E) RESULT(DN_SCAT)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE DERIVATIVE INTENSITY OF THE SCATTERED CONTINUUM   #
        !#OF THE TUBE AS DESCRIBED IN:                                                   #
        !#  ZARKADAS, Ch; KARYDAS, A. G.; PARADELLIS, Th. Theoretical study of           #
        !#  a secondary target XRF setup at different operational tube voltages.         #
        !#  X‐Ray Spectrometry, 2001, 30.2: 99-109.                                      #
        !#################################################################################
        !#INPUTS:                                                                        #
        !#      -E      (DBLE)      ENERGY                          [keV]                #
        !#################################################################################
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  WP = selected_real_kind(8)
        REAL(DP), INTENT(IN) :: E
        REAL(WP) :: DN_SCAT

        REAL(WP) :: CONST
        REAL(WP) :: N_SR
        REAL(WP) :: N_SC
        REAL(WP) :: C_R
        REAL(WP) :: C_C
        REAL(WP) :: I_R
        REAL(WP) :: I_C
        REAL(WP) :: DCS
        REAL(DP) :: EI

        I_R = DERIV_CONT(E, Z_ANODE)*DCS_R(E, A_ST_POL)
        C_R = CALC_CHI(E, A_ST_INCID, E, A_ST_TAKE_OFF)

        N_SR = I_R/C_R

        EI = COMPTON_ENERGY(E, A_ST_POL)

        I_C = DERIV_CONT(EI, Z_ANODE)*DCS_C(EI, A_ST_POL)
        C_C = CALC_CHI(EI, A_ST_INCID, E, A_ST_TAKE_OFF)

        N_SC = I_C/C_C

        DN_SCAT = (N_SR + N_SC)
        RETURN
    END FUNCTION DERIV_ST_SCAT
END MODULE SECCOMP
