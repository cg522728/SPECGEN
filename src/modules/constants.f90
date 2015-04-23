MODULE CONSTANTS
    USE :: ISO_FORTRAN_ENV
    USE :: TYPES
    USE :: xraylib
    USE :: XRLDATA
    USE :: CFGDATA
    implicit none

CONTAINS
    FUNCTION DETEFF(EI)
        implicit none
        REAL(QP)   :: DETEFF
        REAL(DP), INTENT(IN)    :: EI
        REAL(QP)    :: K1
        REAL(QP)    :: K2
        K1 = MAC(Z_DET_WINDOW, EI)*D_DET_WINDOW*1E-4&
                +MAC(Z_DET_DL, EI)*D_DET_DL*1E-4&
                +MAC(Z_DET_GAP, EI)*D_DET_GAP*1E-4
        K2 = MAC(Z_DET_BODY, EI)*D_DET_BODY*1E-4
        DETEFF = EXP(-K1)*(1-EXP(-K2))
        RETURN
    END FUNCTION
END MODULE CONSTANTS
