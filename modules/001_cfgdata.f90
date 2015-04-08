MODULE CFGDATA
    USE :: xraylib
    IMPLICIT NONE
    PUBLIC
    REAL(16)                     :: VTUBE
    REAL(16)                     :: ITUBE
    REAL(16)                     :: ESTEP, EMIN
    REAL(16)                     :: CONC
    INTEGER                     :: NSTEP
    INTEGER                     :: Z_ANODE
    INTEGER                     :: Z_WINDOW
    INTEGER                     :: Z_FILTER
    REAL(16)                     :: D_WINDOW
    REAL(16)                     :: D_FILTER
    REAL(16)                     :: A_TAKE_OFF
    REAL(16)                     :: A_ST_POL
    REAL(16)                     :: A_ST_AZIM_IN
    REAL(16)                     :: A_ST_AZIM_OUT
    REAL(16)                     :: A_SAM_IN
    REAL(16)                     :: A_SAM_OUT
    REAL(16)                     :: SA_ANODE_OUT
    REAL(16)                     :: SA_ST_IN
    REAL(16)                     :: SA_ST_OUT
    CHARACTER(LEN=16)           :: STR_SECTARGET
    CHARACTER(LEN=16)           :: STR_SAMPLE
    CHARACTER(LEN=16)           :: STR_FILTER
    CHARACTER(LEN=3)            :: STR_TYPE
    TYPE(CompoundData), POINTER :: CP_ST
    TYPE(CompoundData), POINTER :: CP_SAM

    CHARACTER(LEN=255)          :: STRING

CONTAINS

    !======================================================================================================
    SUBROUTINE INITCFG()
        !#################################################################################
        !#THIS SUBROUTINE CHECKS IF THE CONFIG FILES ARE                                 #
        !#PRESENT. IF THEY ARE NOT PRESENT IT CREATES                                    #
        !#AND FILLS THEM WITH THE STANDARD DATA FOR                                      #
        !#THE PANALYTICAL EPSILON 5                                                      #
        !#################################################################################
        USE :: xraylib
        IMPLICIT NONE
        LOGICAL :: DIR_E

        INQUIRE(FILE='tube.dat', EXIST=DIR_E)
        IF (.NOT.DIR_E) THEN
            OPEN(UNIT=1,FILE='tube.dat')
            WRITE (1, 200) 'PAN-Gd', 64, 26., 1., 4, 300.
            CLOSE(1)
        ENDIF
        INQUIRE(FILE='sect.dat', EXIST=DIR_E)
        IF (.NOT.DIR_E) THEN
            OPEN(UNIT=1,FILE='sect.dat')
            WRITE (1, 202) 'Al', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Ti', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Fe', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Ge', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Zr', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Mo', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Ag', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'CsI', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Ce2O3', 'SCT', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'B4C', 'POL', 90., 45., 45., 1., 1.
            WRITE (1, 202) 'Al2O3', 'POL', 90., 45., 45., 1., 1.
            CLOSE(1)
        ENDIF
        INQUIRE(FILE='filt.dat', EXIST=DIR_E)
        IF (.NOT.DIR_E) THEN
            OPEN(UNIT=1,FILE='filt.dat')
            WRITE (1, 203) 'Mo', 250.
            WRITE (1, 203) 'Zr', 125.
            WRITE (1, 203) 'Cu', 100.
            WRITE (1, 203) 'Al', 500.
            WRITE (1, 203) 'Al', 250.
            CLOSE(1)
        ENDIF
200     FORMAT(A16, I3, 2F6.2, I3, F6.2)
202     FORMAT(A16, 2X, A3, 5F6.2)
203     FORMAT(A16, F6.2)
    END SUBROUTINE INITCFG
    !======================================================================================================
    SUBROUTINE READCFG()
        !#################################################################################
        !#THIS SUBROUTINE READS THE BASIC DATA FROM THE COMMANDLINE                      #
        !#ARGUMENTS. THE PARAMETERS IT READS ARE (IN ORDER):                             #
        !#   -TUBE POTENTIAL IN kV                                                       #
        !#   -TUBE CURRENT IN mA                                                         #
        !#   -THE SECONDARY TARGET USED IN SYMB.                                         #
        !#   -THE FILTER USED IN TERMS OF ATOMIC NUMBER                                  #
        !#   -THE THICKNESS OF THE FILTER IN mm                                          #
        !#   -THE COMPOSITION OF THE SAMPLE IN SYMB.                                     #
        !#   -THE CONCENTRATION OF THE SAMPLE IN ug/mm^2                                 #
        !#   -THE STEPSIZE FOR THE SPECTRUM AND INTEGRATIONS                             #
        !#   -THE MINIMAL ENERGY OF THE SPECTRUM                                         #
        !#THE EXPRESSIONS IN SYMB. ARE EVALUATED WITH THE COMPOUND                       #
        !#PARSER, WHICH IS PART OF XRAYLIB.                                              #
        !#################################################################################
        USE :: xraylib
        IMPLICIT NONE
        CHARACTER(LEN=16)   :: ARG1, ARG2, ARG3, ARG4, ARG5, ARG6, ARG7
        CHARACTER(LEN=16)   :: ARG8, ARG9
        REAL(16) :: PI
        PI = 2.D0*DASIN(1.D0)

        !READING THE ARGUMENTS
        !ARGUMENT 1 IS THE SELECTOR!
        CALL GET_COMMAND_ARGUMENT(2, ARG1)! VTUBE
        CALL GET_COMMAND_ARGUMENT(3, ARG2)! ITUBE
        CALL GET_COMMAND_ARGUMENT(4, ARG3)! STR_SECTARGER
        CALL GET_COMMAND_ARGUMENT(5, ARG4)! Z_FILTER
        CALL GET_COMMAND_ARGUMENT(6, ARG5)! D_FILTER
        CALL GET_COMMAND_ARGUMENT(7, ARG6)! STR_SAMPLE
        CALL GET_COMMAND_ARGUMENT(8, ARG7)! CONC
        CALL GET_COMMAND_ARGUMENT(9, ARG8)! ESTEP
        CALL GET_COMMAND_ARGUMENT(10, ARG9)! EMIN

        !REMOVING SPACES FROM THE READ STRINGS
        ARG1 = ADJUSTL(ARG1)
        ARG2 = ADJUSTL(ARG2)
        ARG3 = ADJUSTL(ARG3)
        ARG4 = ADJUSTL(ARG4)
        ARG5 = ADJUSTL(ARG5)
        ARG6 = ADJUSTL(ARG6)
        ARG7 = ADJUSTL(ARG7)
        ARG8 = ADJUSTL(ARG8)
        ARG9 = ADJUSTL(ARG9)

        !CONVERTING STRING TO NUMERIC WHERE NECCESSARY
        READ(ARG1,'(F6.2)') VTUBE
        READ(ARG2,'(F6.2)') ITUBE
        READ(ARG4,'(I3)') Z_FILTER
        READ(ARG5,'(F6.2)') D_FILTER
        READ(ARG7,'(F6.2)') CONC
        READ(ARG8,'(F6.2)') ESTEP
        READ(ARG9,'(F6.2)') EMIN
        STR_SECTARGET = ARG3
        STR_SAMPLE = ARG6

        !SETTING DEFAULT VALUES FOR PANALYTICAL EPSILON 5
        !TUBE AND GEOMETRY
        Z_ANODE = 64
        A_TAKE_OFF = 26
        SA_ANODE_OUT = 1
        Z_WINDOW = 4
        D_WINDOW = 300.
        A_ST_POL = 45
        A_ST_AZIM_IN = 45
        A_ST_AZIM_OUT = 45
        SA_ST_IN = 1
        SA_ST_OUT = 1

        !CONVERTING ANGLES IN DEGREES TO ANGLES IN RADIAN
        A_TAKE_OFF = DEG2RAD(A_TAKE_OFF)
        A_ST_POL = DEG2RAD(A_ST_POL)
        A_ST_AZIM_IN = DEG2RAD(A_ST_AZIM_IN)
        A_ST_AZIM_OUT = DEG2RAD(A_ST_AZIM_OUT)
        !CALCULATING NUMBER OF STEPS IN CONTINUUM AND INTEGRATION
        NSTEP = INT((VTUBE-EMIN)*ESTEP**(-1))
        !PARSING COMPOUND STRINGS WITH COMPOUNDPARSER FROM XRAYLIB
        CP_ST => COMPOUNDPARSER(ADJUSTL(STR_SECTARGET))
        CP_SAM => COMPOUNDPARSER(ADJUSTL(STR_SAMPLE))
    END SUBROUTINE READCFG
    !======================================================================================================
    SUBROUTINE LOADCFG()
        !#################################################################################
        !#THIS SUBROUTINE READS THE STORED CONFIGURATIONS                                #
        !#FOR TUBE, SECONDARY TARGET, AND FILTER.                                        #
        !#IT OFFERS A MENU FORM WHICH THE USER CAN CHOOSE                                #
        !#EITHER A STORED CONFIGURATION OF TO CREATE A NEW                               #
        !#CONFIGURATION                                                                  #
        !#################################################################################
        USE :: xraylib
        IMPLICIT NONE
        CHARACTER(LEN=16) :: NAME
        INTEGER :: CHOICE, NUM, CNT
        !CHECKING THE NUMBER OF CONFIGURATIONS STORED IN
        !'tube.dat'
        OPEN(UNIT=1, FILE='tube.dat',ACCESS='SEQUENTIAL')
        NUM = 0
        DO
            READ (1,*,END=1) STRING
            NUM = NUM + 1
        END DO

        !DISPLAY MENU OF TUBE-CONFIGURATIONS AND ASKING FOR CHOICE
1       WRITE (*,*) "TUBE CONFIG IN FILE"
        REWIND(1)
        DO CNT = 1, NUM
            READ (1,200,END=2) STRING, Z_ANODE, A_TAKE_OFF, SA_ANODE_OUT,&
                Z_WINDOW, D_WINDOW
            WRITE (*,201) CNT, STRING
        END DO
2       WRITE (*,201) 0, 'NEW CONFIG'
        WRITE (*,'(A,$)') 'CHOICE='
        READ (*,*) CHOICE
        IF (CHOICE.NE. 0) GO TO 3
        CLOSE(1)
        !DESIGN AND WRITE NEW CONFIGURATION TO FILE
        !OPENING CONFIG FILE 'tube.dat' FOR APPEND
        OPEN(UNIT=1, FILE='tube.dat', ACCESS='APPEND')
        WRITE (*,'(A,$)') "NAME = "
        READ (*,*) NAME
        WRITE (*,'(A,$)') "ANODE (SYMB) = "
        READ (*,*) STRING
        Z_ANODE = SymbolToAtomicNumber(STRING)  !CONVERT SYMB. TO ATOMIC NUMBER
        WRITE (*,'(A,$)') "TAKE OFF ANGLE (IN DEG) = "
        READ (*,*) A_TAKE_OFF
        WRITE (*,'(A,$)') "SOLID ANGLE (IN DEG) = "
        READ (*,*) SA_ANODE_OUT
        WRITE (*,'(A,$)') "Z_WINDOW = "
        READ (*,*) STRING
        Z_WINDOW = SYMBOLTOATOMICNUMBER(STRING)
        WRITE (*,'(A,$)') "D_WINDOW (IN mm) = "
        READ (*,*) D_WINDOW
        !WRITE CONFIGURATION TO FILE 'tube.dat'
        WRITE (1,200) NAME, Z_ANODE, A_TAKE_OFF, SA_ANODE_OUT,&
            Z_WINDOW, D_WINDOW
        GO TO 4
        !IF CHOICE IS NOT 0, THEN REWIND FILE AND READ UNTIL CORRECT
        !CONFIGURATION HAS BEEN READ
3       REWIND(1)
        DO CNT = 1, CHOICE
            READ (1,200,END=3) NAME, Z_ANODE, A_TAKE_OFF, SA_ANODE_OUT,&
                Z_WINDOW, D_WINDOW
        END DO
        CLOSE(1)
        !CHECKING THE NUMBER OF CONFIGURATIONS STORED IN
        !'sect.dat'
4       OPEN(UNIT=2, FILE='sect.dat',ACCESS='SEQUENTIAL')
        NUM = 0
        DO
            READ (2,*,END=5) STRING
            NUM = NUM + 1
        END DO

5       WRITE (*,*) "SECONDARY TARGET CONFIG IN FILE"
        REWIND(2)
        !DISPLAY MENU OF ST-CONFIGURATIONS AND ASKING FOR CHOICE
        DO CNT = 1, NUM
            READ (2,202,END=6) STR_SECTARGET, STR_TYPE, A_ST_POL, A_ST_AZIM_IN, A_ST_AZIM_OUT,&
            SA_ST_IN, SA_ST_OUT
            WRITE (*,201) CNT, STR_SECTARGET
        END DO
6       WRITE (*,201) 0, 'NEW CONFIG'
        WRITE (*,'(A,$)') 'CHOICE='
        READ (*,*) CHOICE
        IF (CHOICE.NE. 0) GO TO 7
        CLOSE(2)
        !DESIGN AND WRITE NEW CONFIGURATION TO FILE
        !OPENING CONFIG FILE 'sect.dat' FOR APPEND
        OPEN(UNIT=2, FILE='sect.dat', ACCESS='APPEND')
        WRITE (*,'(A,$)') "SECONDARY TARGET (SYMB) = "
        READ (*,*) STR_SECTARGET
        WRITE (*,'(A,$)') "TYPE (POL/SCT)) = "
        READ (*,*) STR_TYPE
        WRITE (*,'(A,$)') "POLAR ANGLE (IN DEG) = "
        READ (*,*) A_ST_POL
        WRITE (*,'(A,$)') "AZIMUTHAL ANGLE ENTRY (IN DEG) = "
        READ (*,*) A_ST_AZIM_IN
        WRITE (*,'(A,$)') "AZIMUTHAL ANGLE EXIT (IN DEG) = "
        READ (*,*) A_ST_AZIM_OUT
        WRITE (*,'(A,$)') "SOLID ANGLE ENTRY (IN SR) = "
        READ (*,*) SA_ST_IN
        WRITE (*,'(A,$)') "SOLID ANGLE EXIT (IN SR) = "
        READ (*,*) SA_ST_OUT
        !WRITE CONFIGURATION TO FILE 'sect.dat'
        WRITE (2,202) STR_SECTARGET, STR_TYPE, A_ST_POL, A_ST_AZIM_IN, A_ST_AZIM_OUT,&
            SA_ST_IN, SA_ST_OUT
        GO TO 8
7       REWIND(2)
        !IF CHOICE IS NOT 0, THEN REWIND FILE AND READ UNTIL CORRECT
        !CONFIGURATION HAS BEEN READ
        DO CNT = 1, CHOICE
            READ (2,202,END=7) STR_SECTARGET, STR_TYPE, A_ST_POL, A_ST_AZIM_IN, A_ST_AZIM_OUT, SA_ST_IN,&
                SA_ST_OUT
        END DO
8       CLOSE(2)
        !CHECKING THE NUMBER OF CONFIGURATIONS STORED IN
        !'filt.dat'
        OPEN(UNIT=3, FILE='filt.dat',ACCESS='SEQUENTIAL')
        NUM = 0
        DO
            READ (3,*,END=10) STRING
            NUM = NUM + 1
        END DO

10      WRITE (*,*) "FILTER CONFIG IN FILE"
        REWIND(3)
        !DISPLAY MENU OF ST-CONFIGURATIONS AND ASKING FOR CHOICE
        DO CNT = 1, NUM
            READ (3,203,END=11) STR_FILTER, D_FILTER
            WRITE (*,201) CNT, STR_FILTER
        END DO
11      WRITE (*,201) 0, 'NEW CONFIG'
        WRITE (*,201) 1, 'NO FILTER'
        WRITE (*,'(A,$)') 'CHOICE='
        READ (*,*) CHOICE
        !IF CHOICE IS 1, THEN NO FILTER IS REQUIRED
        !TO REMOVE THE FILTER FROM ALL EQUATIONS,
        !THE THICKNESS IS SET TO 0 mm.
        !THIS LEADS TO ALL FILTER CONTRIBUTIONS
        !BEING EQUAL TO 1
        !EXP(Mu*0)=1
        IF (CHOICE.EQ. 1) THEN
            Z_FILTER = 4
            D_FILTER = 0
            GO TO 13
        ENDIF
        IF (CHOICE.NE. 0) GO TO 12
        CLOSE(3)
        !DESIGN AND WRITE NEW CONFIGURATION TO FILE
        !OPENING CONFIG FILE 'filt.dat' FOR APPEND
        OPEN(UNIT=3, FILE='filt.dat', ACCESS='APPEND')
        WRITE (*,'(A,$)') "FILTER (SYMB) = "
        READ (*,*) STR_FILTER
        Z_FILTER = SYMBOLTOATOMICNUMBER(STR_FILTER)
        WRITE (*,'(A,$)') "D_FILTER = "
        READ (*,*) D_FILTER
        !WRITE CONFIGURATION TO FILE 'sect.dat'
        WRITE (3,203) STR_FILTER, D_FILTER
        GO TO 8
12      REWIND(2)
        !IF CHOICE IS NOT 0, THEN REWIND FILE AND READ UNTIL CORRECT
        !CONFIGURATION HAS BEEN READ
        DO CNT = 1, CHOICE
            READ (3,203,END=12) STR_FILTER, D_FILTER
        END DO
13      CLOSE(3)
        WRITE (*,'(A,$)') "SAMPLE (SYMB) = "
        READ (*,*) STR_SAMPLE
        WRITE (*,'(A,$)') "CONC (IN ug/mm^2) = "
        READ (*,*) CONC
        WRITE (*,'(A,$)') "V_TUBE (IN kV) = "
        READ (*,*) VTUBE
        WRITE (*,'(A,$)') "I_TUBE (IN mA) = "
        READ (*,*) ITUBE
        WRITE (*,'(A,$)') "E_STEP (IN keV) = "
        READ (*,*) ESTEP
        WRITE (*,'(A,$)') "E_MIN (IN keV) = "
        READ (*,*) EMIN
        CLOSE(1)
200     FORMAT(A16, I3, 2F6.2, I3, F6.2)
201     FORMAT(T10, 1H[, I2, 1H], 2X, A16)
202     FORMAT(A16, 2X, A3, 5F6.2)
203     FORMAT(A16, F6.2)

        !CONVERTING ANGLES IN DEGREES TO ANGLES IN RADIAN
        A_TAKE_OFF = DEG2RAD(A_TAKE_OFF)
        A_ST_POL = DEG2RAD(A_ST_POL)
        A_ST_AZIM_IN = DEG2RAD(A_ST_AZIM_IN)
        A_ST_AZIM_OUT = DEG2RAD(A_ST_AZIM_OUT)

        !CALCULATING NUMBER OF STEPS IN CONTINUUM AND INTEGRATION
        NSTEP = INT((VTUBE-EMIN)*ESTEP**(-1))

        !PARSING COMPOUND STRINGS WITH COMPOUNDPARSER FROM XRAYLIB
        STR_SECTARGET = ADJUSTL(STR_SECTARGET)
        CP_ST => CompoundParser(STR_SECTARGET)
        RETURN
    END SUBROUTINE LOADCFG
    !======================================================================================================
    FUNCTION DEG2RAD(ANG)
        !#################################################################################
        !#THIS FUNCTION CONVERTS ANGLES IN DEGREES                                       #
        !#TO ANGLES IN RADIAN                                                            #
        !#################################################################################
        IMPLICIT NONE
        REAL(16) :: ANG
        REAL(16) :: DEG2RAD
        REAL(16) :: PI
        PI = 2.D0*DASIN(1.D0)
        DEG2RAD = ANG*(PI/180)
    END FUNCTION DEG2RAD
    !======================================================================================================
    SUBROUTINE DEBUGCFG()
        !#################################################################################
        !#WRITING ALL POSSIBLE PARAMETERS TO 'DEBUG.LOG'                                 #
        !#################################################################################
        WRITE (999,*) Z_ANODE
        WRITE (999,*) A_TAKE_OFF
        WRITE (999,*)  SA_ANODE_OUT
        WRITE (999,*) Z_WINDOW
        WRITE (999,*)  D_WINDOW
        WRITE (999,*) STR_SECTARGET
        WRITE (999,*) A_ST_POL
        WRITE (999,*) A_ST_AZIM_IN
        WRITE (999,*) A_ST_AZIM_OUT
        WRITE (999,*) SA_ST_IN
        WRITE (999,*) SA_ST_OUT
        WRITE (999,*) STR_FILTER
        WRITE (999,*) D_FILTER
        WRITE (999,*) VTUBE
        WRITE (999,*) ITUBE
        WRITE (999,*) ESTEP
        WRITE (999,*) NSTEP
        WRITE (999,*) EMIN
    END SUBROUTINE DEBUGCFG
END MODULE CFGDATA
