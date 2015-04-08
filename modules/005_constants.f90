MODULE CONSTANTS
    USE     :: ISO_FORTRAN_ENV
    USE :: xraylib
    USE :: XRLDATA
    USE :: CFGDATA
    USE :: ANODE

CONTAINS
    FUNCTION CHECKKAPPA(DELTA)
        !#################################################################################
        !#THIS FUNCTION CALCULATES A VALUE FOR THE PROPORTIONALITY CONSTANT KAPPA        #
        !#BY CALCULATING CALC_KR() FOR EACH LINE AND INCREASING KAPPA BY 2^DELTA EACH    #
        !#LOOP. THIS LOOP CONTINUES UNTIL THE CHARACTERISTIC LINE INTENSITY IS GREATER   #
        !#THAN THE CONTINUUM.                                                            #
        !#################################################################################
        USE     :: ISO_FORTRAN_ENV
        IMPLICIT NONE
        INTEGER(INT64) :: N, BASE
        INTEGER :: CNT
        REAL(16) :: EI, EC
        REAL(16) :: TAU
        REAL(16) :: R, PI, DELTA
        REAL(16) :: KAPPA
        REAL(16) :: CHECKKAPPA
        REAL(16)    :: CON, K, NR

        DATA BASE/2/

        PI = 2.D0*DASIN(1.D0)
        K = 1.D0
        KAPPA = K
        DO N = 1, SIZE(LINE)
        WRITE (6,'(I4,1H/,I4,A1,$)',ADVANCE='NO') N, SIZE(LINE), CHAR(13)
        NR = N
        EC = ANINT(LineEnergy(Z_ANODE, LINE(N))/ESTEP)*ESTEP
        TAU = (((VTUBE/EC)*LOG((VTUBE/EC)/((VTUBE/EC)-1)))-1)
1       R = CALC_KR(Z_ANODE, EC, INT(N), K)*TAU
        IF (R.EQ. 0) CYCLE
        IF (R.LT. 1) THEN
            CON = 2_16**DELTA
            K = K + CON
            GO TO 1
        ENDIF
        IF (K.GE. KAPPA) KAPPA = K
        END DO
        CHECKKAPPA = KAPPA
        RETURN
    END FUNCTION CHECKKAPPA
END MODULE CONSTANTS
