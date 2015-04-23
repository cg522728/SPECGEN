MODULE XRLDATA
    USE :: xraylib
    USE :: TYPES
    implicit none
    INTEGER, DIMENSION(111)   :: SHELL
    INTEGER, DIMENSION(111)   :: LINE
    INCLUDE 'shell.f90'
    INCLUDE 'line.f90'
CONTAINS
    FUNCTION CS_FLUOR_CHG(Z, N, EI)
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        REAL(QP)    :: CS_FLUOR_CHG
        REAL(QP)    :: CS
        REAL(DP)    :: EI
        REAL(DP)    :: E_EDGE
        INTEGER     :: Z
        INTEGER     :: N

        E_EDGE = EdgeEnergy(Z, SHELL(N))
        IF (EI.LT.E_EDGE) THEN
            CS_FLUOR_CHG = 0_DP
            RETURN
        ENDIF
        CALL SetErrorMessages(0)
        CS = CS_Fluorline_Kissel_Cascade(Z, LINE(N), EI)
        CALL SetErrorMessages(1)
        CS_FLUOR_CHG = CS
        RETURN
    END FUNCTION CS_FLUOR_CHG
    FUNCTION FLUORYIELD_CHG(Z, SHELLIN)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE FLUORESCENCE YIELD                                #
        !#OF AN ATOMIC SHELL USING THE FORMULAS                                          #
        !#ACCORDING TO HUBBEL ET AL (1994) AS                                            #
        !#LISTED IN 'HANDBOOK OF X-RAY SPECTROMETRY' 2ND ED                              #
        !#(PRACTICAL SPECTROSCOPY SERIES VOLUME 29)                                      #
        !#USING THE FOLLOWING INPUTS:                                                    #
        !#  -Z          (INT)   :ATOMIC NUMBER                                           #
        !#  -SHELLIN    (INT)   :SHELL IDENTIFIER AS DEFINED IN XRAYLIB                  #
        !#################################################################################
        USE :: xraylib
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        REAL(QP), PARAMETER, DIMENSION(0:3) :: CK = (/0.0370, 0.03112, 5.44E-5, -1.25E-6/)
        REAL(QP), PARAMETER, DIMENSION(0:3) :: CL = (/0.17765, 0.00298937, 8.91297E-5, -2.67184E-7/)
        REAL(QP), PARAMETER, DIMENSION(2)   :: CM = (/1.29E-9, 13./)
        INTEGER :: Z
        INTEGER :: SHELLIN
        REAL(QP) :: FLUORYIELD_CHG
        REAL(QP) :: OMEGA1, OMEGA2
        INTEGER :: N

        DATA    OMEGA1/0/
        DATA    OMEGA2/0/

        IF (SHELLIN.EQ. K_SHELL) GO TO 101
        IF (SHELLIN.EQ. L3_SHELL) GO TO 102
        IF (SHELLIN.EQ. L2_SHELL) GO TO 102
        IF (SHELLIN.EQ. L1_SHELL) GO TO 102
        IF (SHELLIN.EQ. M5_SHELL) GO TO 103
        IF (SHELLIN.EQ. M4_SHELL) GO TO 103
        IF (SHELLIN.EQ. M3_SHELL) GO TO 103
        IF (SHELLIN.EQ. M2_SHELL) GO TO 103
        IF (SHELLIN.EQ. M1_SHELL) GO TO 103
        FLUORYIELD_CHG = 0
        RETURN

101     N = 0
        DO N=0,3
            OMEGA1 = OMEGA1 + CK(N)*Z**N
        END DO
        FLUORYIELD_CHG = OMEGA1**4*(1+OMEGA1**4)**(-1)
        RETURN
102     IF (Z.GE. 3 .AND. Z.LE.36) GO TO 2
        IF (Z.GE. 37 .AND. Z.LE.100) GO TO 3
2       FLUORYIELD_CHG = 1.939E-8*Z**(3.8874)
        RETURN
3       N = 0
        DO N=0,3
            OMEGA1 = OMEGA1 + CL(N)*Z**N
        END DO
        FLUORYIELD_CHG = OMEGA1**4*(1+OMEGA1**4)**(-1)
        RETURN
103     FLUORYIELD_CHG = CM(1)*(Z-CM(2))**4
        RETURN
    END FUNCTION
    FUNCTION MAC(Z, E)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE MASS ABSORPTION COEFFICIENT FOR AN ELEMENT        #
        !#BASED ON THE FOLLOWING INPUTS:                                                 #
        !#  -Z      (INT)   :ATOMIC NUMBER                                               #
        !#  -E      (DBLE)  :ENERGY                                                      #
        !#THE CALCULATIONS ARE PERFORMED USING THE TOTAL CROSS SECTIONS IN BARNS/ATOM    #
        !#AS DELIVERED BY THE XRAYLIBFUNCTION CSb_Total()                                #
        !#################################################################################
        USE :: xraylib
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        INTEGER :: Z
        REAL(DP) :: E
        REAL(QP), PARAMETER :: U = 1.6605402E-24
        REAL(QP) :: MAC

        MAC = CS_TOTAL(Z, E)
        RETURN
    END FUNCTION
    FUNCTION MAC_COMP(CP, E)
        !#################################################################################
        !#THIS FUNCTION CALCULATES THE MASS ABSORPTION COEFFICIENT FOR A COMPOUND        #
        !#BASED ON THE FOLLOWING INPUTS:                                                 #
        !#  -Z      (INT)   :ATOMIC NUMBER                                               #
        !#  -E      (DBLE)  :ENERGY                                                      #
        !#THE CALCULATIONS ARE PERFORMED USING THE TOTAL CROSS SECTIONS IN BARNS/ATOM    #
        !#AS DELIVERED BY THE XRAYLIBFUNCTION CSb_Total()                                #
        !#################################################################################
        USE :: xraylib
        USE :: CFGDATA
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        INTEGER :: Z, CNT
        REAL(DP) :: E
        REAL(QP), DIMENSION(:), ALLOCATABLE :: MAC_TMP
        REAL(QP) :: MAC_COMP
        TYPE(CompoundData), POINTER :: CP

        ALLOCATE(MAC_TMP(CP_ST%NELEMENTS))
        DO CNT=1,CP%NELEMENTS
            Z = CP%ELEMENTS(CNT)
            MAC_TMP(CNT) = MAC(Z, E)*CP%MASSFRACTIONS(CNT)
        END DO
        MAC_COMP = SUM(MAC_TMP)
        MAC_COMP = MAC_COMP
        RETURN
    END FUNCTION
        FUNCTION CS_FLUOR_CP_CHG(CP, N, EC)
        USE :: xraylib
        USE :: CFGDATA
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        INTEGER :: Z, CNT, N
        REAL(DP) :: EC, E_EDGE
        REAL(QP), PARAMETER :: U = 1.6605402E-24
        REAL(QP), DIMENSION(:), ALLOCATABLE :: TMP
        REAL(QP) :: CS_FLUOR_CP_CHG
        TYPE(CompoundData), POINTER :: CP

        ALLOCATE(TMP(CP%NELEMENTS))
        DO CNT=1,CP%NELEMENTS
            Z = CP%ELEMENTS(CNT)
            E_EDGE = EdgeEnergy(Z, SHELL(N))
            IF (EC.LT.E_EDGE) THEN
                TMP(CNT) = 0_QP
                CYCLE
            ENDIF
            CALL SetErrorMessages(0)
            TMP(CNT) = CS_Fluorline_Kissel_Cascade(CP%ELEMENTS(CNT), LINE(N), DBLE(EC))&
                        *CP%MASSFRACTIONS(CNT)
            CALL SetErrorMessages(1)
        END DO
        CS_FLUOR_CP_CHG = SUM(TMP)
        CS_FLUOR_CP_CHG = CS_FLUOR_CP_CHG
    END FUNCTION
    FUNCTION TRANSPROB(Z, N)
        IMPLICIT NONE
        !f2py INTEGER, PARAMETER ::  QP = selected_real_kind(8)
        !f2py INTEGER, PARAMETER ::  DP = selected_real_kind(8)
        REAL(QP)    :: TRANSPROB
        INTEGER     :: N
        INTEGER     :: Z
        INTEGER     :: CNT
        REAL(QP)    :: RIJ
        REAL(QP), DIMENSION(:), ALLOCATABLE    :: T
        ALLOCATE(T(SIZE(LINE)))
        T = 0._QP

        DO CNT = 1, SIZE(LINE)
            T(CNT) = RadRate(Z, LINE(CNT))
        END DO
        RIJ = RadRate(Z, N)/SUM(T)
        TRANSPROB = RIJ
        RETURN
    END FUNCTION TRANSPROB
END MODULE XRLDATA
