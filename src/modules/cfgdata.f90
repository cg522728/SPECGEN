MODULE CFGDATA
    USE :: xraylib
    USE :: TYPES
    IMPLICIT NONE
    PUBLIC
    REAL(DP)                     :: VTUBE
    REAL(QP)                     :: ITUBE
    REAL(DP)                     :: ESTEP, EMIN
    REAL(QP)                     :: CONC
    INTEGER                     :: NSTEP
    INTEGER                     :: Z_ANODE
    INTEGER                     :: Z_WINDOW
    INTEGER                     :: Z_FILTER
    REAL(QP)                     :: D_WINDOW
    REAL(QP)                     :: D_FILTER
    REAL(QP)                     :: A_INCID
    REAL(QP)                     :: A_TAKE_OFF
    REAL(QP)                     :: A_ST_POL
    REAL(QP)                     :: A_ST_AZIM_IN
    REAL(QP)                     :: A_ST_AZIM_OUT
    REAL(QP)                     :: A_SAM_IN
    REAL(QP)                     :: A_SAM_OUT
    REAL(QP)                     :: SA_ANODE_OUT
    REAL(QP)                     :: SA_ST_IN
    REAL(QP)                     :: SA_ST_OUT
    REAL(QP)                    :: D_DET_WINDOW
    INTEGER                    :: Z_DET_WINDOW
    REAL(QP)                    :: D_DET_GAP
    INTEGER                    :: Z_DET_GAP
    INTEGER                    :: Z_DET_DL
    REAL(QP)                    :: D_DET_DL
    REAL(QP)                    :: D_DET_BODY
    INTEGER                    :: Z_DET_BODY
    CHARACTER(LEN=16)           :: STR_SECTARGET
    CHARACTER(LEN=16)           :: STR_SAMPLE
    CHARACTER(LEN=16)           :: STR_FILTER
    CHARACTER(LEN=16)           :: STR_ANODE
    CHARACTER(LEN=3)            :: STR_TYPE
    TYPE(CompoundData), POINTER :: CP_ST
    TYPE(CompoundData), POINTER :: CP_SAM
    TYPE(CompoundData), POINTER :: CP_AN

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
        INQUIRE(FILE='det.dat', EXIST=DIR_E)
        IF (.NOT.DIR_E) THEN
            OPEN(UNIT=1,FILE='det.dat')
            WRITE (1, 204) 'PAN-32', 4, 8., 4, 0., 32, 50.E-3, 32, 5.E3
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
        INQUIRE(FILE='nist.dat', EXIST=DIR_E)
        IF (.NOT.DIR_E) THEN
            OPEN(UNIT=1,FILE='nist.dat')
            WRITE (1,*) '1   H   Hydrogen    0.99212     19.2    8.375E-05'
            WRITE (1,*) '2   He  Helium      0.49968     41.8    1.663E-04'
            WRITE (1,*) '3   Li  Lithium     0.43221     40.0    5.340E-01'
            WRITE (1,*) '4   Be  Beryllium   0.44384     63.7    1.848E+00'
            WRITE (1,*) '5   B   Boron       0.46245     76.0    2.370E+00'
            WRITE (1,*) '6   C   Carbon      0.49954     78.0    1.700E+00'
            WRITE (1,*) '7   N   Nitrogen    0.49976     82.0    1.165E-03'
            WRITE (1,*) '8   O   Oxygen      0.50002     95.0    1.332E-03'
            WRITE (1,*) '9   F   Fluorine    0.47372     115.0   1.580E-03'
            WRITE (1,*) '10  Ne  Neon        0.49555     137.0   8.385E-04'
            WRITE (1,*) '11  Na  Sodium      0.47847     149.0   9.710E-01'
            WRITE (1,*) '12  Mg  Magnesium   0.49373     156.0   1.740E+00'
            WRITE (1,*) '13  Al  Aluminum    0.48181     166.0   2.699E+00'
            WRITE (1,*) '14  Si  Silicon     0.49848     173.0   2.330E+00'
            WRITE (1,*) '15  P   Phosphorus  0.48428     173.0   2.200E+00'
            WRITE (1,*) '16  S   Sulfur      0.49897     180.0   2.000E+00'
            WRITE (1,*) '17  Cl  Chlorine    0.47951     174.0   2.995E-03'
            WRITE (1,*) '18  Ar  Argon       0.45059     188.0   1.662E-03'
            WRITE (1,*) '19  K   Potassium   0.48595     190.0   8.620E-01'
            WRITE (1,*) '20  Ca  Calcium     0.49903     191.0   1.550E+00'
            WRITE (1,*) '21  Sc  Scandium    0.46712     216.0   2.989E+00'
            WRITE (1,*) '22  Ti  Titanium    0.45948     233.0   4.540E+00'
            WRITE (1,*) '23  V   Vanadium    0.45150     245.0   6.110E+00'
            WRITE (1,*) '24  Cr  Chromium    0.46157     257.0   7.180E+00'
            WRITE (1,*) '25  Mn  Manganese   0.45506     272.0   7.440E+00'
            WRITE (1,*) '26  Fe  Iron        0.46556     286.0   7.874E+00'
            WRITE (1,*) '27  Co  Cobalt      0.45815     297.0   8.900E+00'
            WRITE (1,*) '28  Ni  Nickel      0.47708     311.0   8.902E+00'
            WRITE (1,*) '29  Cu  Copper      0.45636     322.0   8.960E+00'
            WRITE (1,*) '30  Zn  Zinc        0.45879     330.0   7.133E+00'
            WRITE (1,*) '31  Ga  Gallium     0.44462     334.0   5.904E+00'
            WRITE (1,*) '32  Ge  Germanium   0.44071     350.0   5.323E+00'
            WRITE (1,*) '33  As  Arsenic     0.44046     347.0   5.730E+00'
            WRITE (1,*) '34  Se  Selenium    0.43060     348.0   4.500E+00'
            WRITE (1,*) '35  Br  Bromine     0.43803     343.0   7.072E-03'
            WRITE (1,*) '36  Kr  Krypton     0.42959     352.0   3.478E-03'
            WRITE (1,*) '37  Rb  Rubidium    0.43291     363.0   1.532E+00'
            WRITE (1,*) '38  Sr  Strontium   0.43369     366.0   2.540E+00'
            WRITE (1,*) '39  Y   Yttrium     0.43867     379.0   4.469E+00'
            WRITE (1,*) '40  Zr  Zirconium   0.43848     393.0   6.506E+00'
            WRITE (1,*) '41  Nb  Niobium     0.44130     417.0   8.570E+00'
            WRITE (1,*) '42  Mo  Molybdenum  0.43777     424.0   1.022E+01'
            WRITE (1,*) '43  Tc  Technetium  0.43919     428.0   1.150E+01'
            WRITE (1,*) '44  Ru  Ruthenium   0.43534     441.0   1.241E+01'
            WRITE (1,*) '45  Rh  Rhodium     0.43729     449.0   1.241E+01'
            WRITE (1,*) '46  Pd  Palladium   0.43225     470.0   1.202E+01'
            WRITE (1,*) '47  Ag  Silver      0.43572     470.0   1.050E+01'
            WRITE (1,*) '48  Cd  Cadmium     0.42700     469.0   8.650E+00'
            WRITE (1,*) '49  In  Indium      0.42676     488.0   7.310E+00'
            WRITE (1,*) '50  Sn  Tin         0.42120     488.0   7.310E+00'
            WRITE (1,*) '51  Sb  Antimony    0.41889     487.0   6.691E+00'
            WRITE (1,*) '52  Te  Tellurium   0.40752     485.0   6.240E+00'
            WRITE (1,*) '53  I   Iodine      0.41764     491.0   4.930E+00'
            WRITE (1,*) '54  Xe  Xenon       0.41130     482.0   5.485E-03'
            WRITE (1,*) '55  Cs  Cesium      0.41383     488.0   1.873E+00'
            WRITE (1,*) '56  Ba  Barium      0.40779     491.0   3.500E+00'
            WRITE (1,*) '57  La  Lanthanum   0.41035     501.0   6.154E+00'
            WRITE (1,*) '58  Ce  Cerium      0.41395     523.0   6.657E+00'
            WRITE (1,*) '59  Pr  Praseodymium    0.41871     535.0   6.710E+00'
            WRITE (1,*) '60  Nd  Neodymium   0.41597     546.0   6.900E+00'
            WRITE (1,*) '61  Pm  Promethium  0.42094     560.0   7.220E+00'
            WRITE (1,*) '62  Sm  Samarium    0.41234     574.0   7.460E+00'
            WRITE (1,*) '63  Eu  Europium    0.41457     580.0   5.243E+00'
            WRITE (1,*) '64  Gd  Gadolinium  0.40699     591.0   7.900E+00'
            WRITE (1,*) '65  Tb  Terbium     0.40900     614.0   8.229E+00'
            WRITE (1,*) '66  Dy  Dysprosium  0.40615     628.0   8.550E+00'
            WRITE (1,*) '67  Ho  Holmium     0.40623     650.0   8.795E+00'
            WRITE (1,*) '68  Er  Erbium      0.40655     658.0   9.066E+00'
            WRITE (1,*) '69  Tm  Thulium     0.40844     674.0   9.321E+00'
            WRITE (1,*) '70  Yb  Ytterbium   0.40453     684.0   6.730E+00'
            WRITE (1,*) '71  Lu  Lutetium    0.40579     694.0   9.840E+00'
            WRITE (1,*) '72  Hf  Hafnium     0.40338     705.0   1.331E+01'
            WRITE (1,*) '73  Ta  Tantalum    0.40343     718.0   1.665E+01'
            WRITE (1,*) '74  W   Tungsten    0.40250     727.0   1.930E+01'
            WRITE (1,*) '75  Re  Rhenium     0.40278     736.0   2.102E+01'
            WRITE (1,*) '76  Os  Osmium      0.39958     746.0   2.257E+01'
            WRITE (1,*) '77  Ir  Iridium     0.40058     757.0   2.242E+01'
            WRITE (1,*) '78  Pt  Platinum    0.39984     790.0   2.145E+01'
            WRITE (1,*) '79  Au  Gold        0.40108     790.0   1.932E+01'
            WRITE (1,*) '80  Hg  Mercury     0.39882     800.0   1.355E+01'
            WRITE (1,*) '81  Tl  Thallium    0.39631     810.0   1.172E+01'
            WRITE (1,*) '82  Pb  Lead        0.39575     823.0   1.135E+01'
            WRITE (1,*) '83  Bi  Bismuth     0.39717     823.0   9.747E+00'
            WRITE (1,*) '84  Po  Polonium    0.40195     830.0   9.320E+00'
            WRITE (1,*) '85  At  Astatine    0.40479     825.0   1.000E+01'
            WRITE (1,*) '86  Rn  Radon       0.38736     794.0   9.066E-03'
            WRITE (1,*) '87  Fr  Francium    0.39010     827.0   1.000E+01'
            WRITE (1,*) '88  Ra  Radium      0.38934     826.0   5.000E+00'
            WRITE (1,*) '89  Ac  Actinium    0.39202     841.0   1.007E+01'
            WRITE (1,*) '90  Th  Thorium     0.38787     847.0   1.172E+01'
            WRITE (1,*) '91  Pa  Protactinium    0.39388     878.0   1.537E+01'
            WRITE (1,*) '92  U   Uranium     0.38651     890.0   1.895E+01'
            CLOSE(1)
        ENDIF
200     FORMAT(A16, I3, 2ES32.20E3, I3, ES32.20E3)
202     FORMAT(A16, 2X, A3, 5ES32.20E3)
203     FORMAT(A16, ES32.20E3)
204     FORMAT(A16, 4(I3, ES32.20E3))
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
        INTEGER     :: Z_SAMPLE
        REAL(16) :: PI
        PI = 2.D0*DASIN(1.D0)

        !READING THE ARGUMENTS
        !ARGUMENT 1 IS THE SELECTOR!
        CALL GET_COMMAND_ARGUMENT(2, ARG1)! VTUBE = 100
        CALL GET_COMMAND_ARGUMENT(3, ARG2)! ITUBE = 6
        CALL GET_COMMAND_ARGUMENT(4, ARG3)! STR_SECTARGER = 'CsI'
        CALL GET_COMMAND_ARGUMENT(5, ARG4)! Z_FILTER = 4
        CALL GET_COMMAND_ARGUMENT(6, ARG5)! D_FILTER = 0.
        CALL GET_COMMAND_ARGUMENT(7, ARG6)! STR_SAMPLE = 'Fe'
        CALL GET_COMMAND_ARGUMENT(8, ARG7)! CONC = 50.
        CALL GET_COMMAND_ARGUMENT(9, ARG8)! ESTEP = 0.01
        CALL GET_COMMAND_ARGUMENT(10, ARG9)! EMIN = 0

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
        STR_FILTER = AtomicNumberToSymbol(Z_FILTER)
        IF (EMIN.EQ.0) EMIN = EMIN + 1E-10

        !SETTING DEFAULT VALUES FOR PANALYTICAL EPSILON 5
        !TUBE AND GEOMETRY
        Z_ANODE = 64
        A_INCID = 90
        A_TAKE_OFF = 26
        SA_ANODE_OUT = 1
        Z_WINDOW = 4
        D_WINDOW = 300.
        A_ST_POL = 90
        A_ST_AZIM_IN = 45
        A_ST_AZIM_OUT = 45
        SA_ST_IN = 1
        SA_ST_OUT = 1

        D_DET_WINDOW = 8.
        Z_DET_WINDOW = 4
        D_DET_GAP = 0.
        Z_DET_GAP = 4
        Z_DET_DL = 32
        D_DET_DL = 50.E-3
        D_DET_BODY = 5.E3
        Z_DET_BODY = 32

        !CONVERTING ANGLES IN DEGREES TO ANGLES IN RADIAN
        A_INCID = DEG2RAD(A_INCID)
        A_TAKE_OFF = DEG2RAD(A_TAKE_OFF)
        A_ST_POL = DEG2RAD(A_ST_POL)
        A_ST_AZIM_IN = DEG2RAD(A_ST_AZIM_IN)
        A_ST_AZIM_OUT = DEG2RAD(A_ST_AZIM_OUT)
        !CALCULATING NUMBER OF STEPS IN CONTINUUM AND INTEGRATION
        NSTEP = INT((VTUBE-EMIN)*ESTEP**(-1))
        !PARSING COMPOUND STRINGS WITH COMPOUNDPARSER FROM XRAYLIB
        CP_ST => COMPOUNDPARSER(ADJUSTL(STR_SECTARGET))
        STR_ANODE = AtomicNumberToSymbol(Z_ANODE)
        CP_AN => COMPOUNDPARSER(STR_ANODE)
        READ (ARG6, '(I3)', ERR=10) Z_SAMPLE
        STR_SAMPLE = AtomicNumberToSymbol(Z_SAMPLE)
        CP_SAM => COMPOUNDPARSER(ADJUSTL(STR_SAMPLE))
        RETURN
10      READ (ARG6, '(A16)') STR_SAMPLE
        CP_SAM => COMPOUNDPARSER(ADJUSTL(STR_SAMPLE))
        RETURN
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
        !'det.dat'
        OPEN(UNIT=1, FILE='det.dat',ACCESS='SEQUENTIAL')
        NUM = 0
        DO
            READ (1,*,END=4) STRING
            NUM = NUM + 1
        END DO

        !DISPLAY MENU OF DETECTOR-CONFIGURATIONS AND ASKING FOR CHOICE
4       WRITE (*,*) "DETECTOR CONFIG IN FILE"
        REWIND(1)
        DO CNT = 1, NUM
            READ (1,204,END=5) STRING, Z_DET_WINDOW, D_DET_WINDOW, Z_DET_GAP, D_DET_GAP,&
 Z_DET_DL, D_DET_DL, Z_DET_BODY, D_DET_BODY
            WRITE (*,201) CNT, STRING
        END DO
5       WRITE (*,201) 0, 'NEW CONFIG'
        WRITE (*,'(A,$)') 'CHOICE='
        READ (*,*) CHOICE
        IF (CHOICE.NE. 0) GO TO 6
        CLOSE(1)
        !DESIGN AND WRITE NEW CONFIGURATION TO FILE
        !OPENING CONFIG FILE 'tube.dat' FOR APPEND
        OPEN(UNIT=1, FILE='det.dat', ACCESS='APPEND')

        WRITE (*,'(A,$)') "NAME = "
        READ (*,*) NAME
        WRITE (*,'(A,$)') "Z_DET_WINDOW = "
        READ (*,*) Z_DET_WINDOW
        WRITE (*,'(A,$)') "D_DET_WINDOW = "
        READ (*,*) D_DET_WINDOW
        WRITE (*,'(A,$)') "Z_DET_GAP = "
        READ (*,*) Z_DET_GAP
        WRITE (*,'(A,$)') "D_DET_GAP = "
        READ (*,*) D_DET_GAP
        WRITE (*,'(A,$)') "Z_DET_DL = "
        READ (*,*) Z_DET_DL
        WRITE (*,'(A,$)') "D_DET_DL = "
        READ (*,*) D_DET_DL
        WRITE (*,'(A,$)') "Z_DET_BODY = "
        READ (*,*) Z_DET_BODY
        WRITE (*,'(A,$)') "D_DET_BODY = "
        READ (*,*) D_DET_BODY

        !WRITE CONFIGURATION TO FILE 'tube.dat'
        WRITE (1,204) NAME, Z_DET_WINDOW, D_DET_WINDOW, Z_DET_GAP, D_DET_GAP,&
 Z_DET_DL, D_DET_DL, Z_DET_BODY, D_DET_BODY
        GO TO 7
        !IF CHOICE IS NOT 0, THEN REWIND FILE AND READ UNTIL CORRECT
        !CONFIGURATION HAS BEEN READ
6       REWIND(1)
        DO CNT = 1, CHOICE
            READ (1,204,END=6) NAME, Z_DET_WINDOW, D_DET_WINDOW, Z_DET_GAP, D_DET_GAP,&
 Z_DET_DL, D_DET_DL, Z_DET_BODY, D_DET_BODY
        END DO
        CLOSE(1)


        !CHECKING THE NUMBER OF CONFIGURATIONS STORED IN
        !'sect.dat'
7       OPEN(UNIT=2, FILE='sect.dat',ACCESS='SEQUENTIAL')
        NUM = 0
        DO
            READ (2,*,END=8) STRING
            NUM = NUM + 1
        END DO

8       WRITE (*,*) "SECONDARY TARGET CONFIG IN FILE"
        REWIND(2)
        !DISPLAY MENU OF ST-CONFIGURATIONS AND ASKING FOR CHOICE
        DO CNT = 1, NUM
            READ (2,202,END=9) STR_SECTARGET, STR_TYPE, A_ST_POL, A_ST_AZIM_IN, A_ST_AZIM_OUT,&
            SA_ST_IN, SA_ST_OUT
            WRITE (*,201) CNT, STR_SECTARGET
        END DO
9       WRITE (*,201) 0, 'NEW CONFIG'
        WRITE (*,'(A,$)') 'CHOICE='
        READ (*,*) CHOICE
        IF (CHOICE.NE. 0) GO TO 10
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
10      REWIND(2)
        !IF CHOICE IS NOT 0, THEN REWIND FILE AND READ UNTIL CORRECT
        !CONFIGURATION HAS BEEN READ
        DO CNT = 1, CHOICE
            READ (2,202,END=10) STR_SECTARGET, STR_TYPE, A_ST_POL, A_ST_AZIM_IN, A_ST_AZIM_OUT, SA_ST_IN,&
                SA_ST_OUT
        END DO
        CLOSE(2)
        !CHECKING THE NUMBER OF CONFIGURATIONS STORED IN
        !'filt.dat'
        OPEN(UNIT=3, FILE='filt.dat',ACCESS='SEQUENTIAL')
        NUM = 0
        DO
            READ (3,*,END=13) STRING
            NUM = NUM + 1
        END DO

13      WRITE (*,*) "FILTER CONFIG IN FILE"
        REWIND(3)
        !DISPLAY MENU OF ST-CONFIGURATIONS AND ASKING FOR CHOICE
        DO CNT = 2, NUM
            READ (3,203,END=14) STR_FILTER, D_FILTER
            WRITE (*,201) CNT, STR_FILTER
        END DO
14      WRITE (*,201) 0, 'NEW CONFIG'
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
            GO TO 16
        ENDIF
        IF (CHOICE.NE. 0) GO TO 15
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
15      REWIND(2)
        !IF CHOICE IS NOT 0, THEN REWIND FILE AND READ UNTIL CORRECT
        !CONFIGURATION HAS BEEN READ
        DO CNT = 2, CHOICE
            READ (3,203,END=15) STR_FILTER, D_FILTER
        END DO
16      CLOSE(3)
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
200     FORMAT(A16, I3, 2ES32.20E3, I3, ES32.20E3)
201     FORMAT(T10, 1H[, I2, 1H], 2X, A16)
202     FORMAT(A16, 2X, A3, 5ES32.20E3)
203     FORMAT(A16, ES32.20E3)
204     FORMAT(A16, 4(I3, ES32.20E3))

        !CONVERTING ANGLES IN DEGREES TO ANGLES IN RADIAN
        A_INCID = 90.
        A_INCID = DEG2RAD(A_INCID)
        A_TAKE_OFF = DEG2RAD(A_TAKE_OFF)
        A_ST_POL = DEG2RAD(A_ST_POL)
        A_ST_AZIM_IN = DEG2RAD(A_ST_AZIM_IN)
        A_ST_AZIM_OUT = DEG2RAD(A_ST_AZIM_OUT)

        !CALCULATING NUMBER OF STEPS IN CONTINUUM AND INTEGRATION
        NSTEP = INT((VTUBE-EMIN)*ESTEP**(-1))

        !PARSING COMPOUND STRINGS WITH COMPOUNDPARSER FROM XRAYLIB
        STR_ANODE = AtomicNumberToSymbol(Z_ANODE)
        CP_AN => COMPOUNDPARSER(STR_ANODE)
        STR_SAMPLE = ADJUSTL(STR_SAMPLE)
        CP_SAM => CompoundParser(STR_SAMPLE)
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
        REAL(QP) :: ANG
        REAL(QP) :: DEG2RAD
        REAL(QP) :: PI
        PI = 2.D0*DASIN(1.D0)
        DEG2RAD = PI*(ANG/180)
    END FUNCTION DEG2RAD
    !======================================================================================================
    SUBROUTINE DEBUGCFG()
        !#################################################################################
        !#WRITING ALL POSSIBLE PARAMETERS TO 'DEBUG.LOG'                                 #
        !#################################################################################
        WRITE (999,*) Z_ANODE
        WRITE (999,*) A_TAKE_OFF
        WRITE (999,*) A_INCID
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
        WRITE (999,*) STR_SAMPLE
        WRITE (999,*) CONC
        WRITE (999,*) VTUBE
        WRITE (999,*) ITUBE
        WRITE (999,*) ESTEP
        WRITE (999,*) NSTEP
        WRITE (999,*) EMIN
    END SUBROUTINE DEBUGCFG
    SUBROUTINE DISPCFG()
        !#################################################################################
        !#DISPLAY ALL POSSIBLE PARAMETERS                                                #
        !#################################################################################
        WRITE (6,'(53("#"))')
        WRITE (6,'(1H#,A16,3H = , I32,1H#)') 'Z_ANODE', Z_ANODE
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'A_INCID', A_INCID
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'A_TAKE_OFF', A_TAKE_OFF
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)')  'SA_ANODE_OUT', SA_ANODE_OUT
        WRITE (6,'(1H#,A16,3H = , I32,1H#)') 'Z_WINDOW', Z_WINDOW
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)')  'D_WINDOW', D_WINDOW
        WRITE (6,'(1H#,A16,3H = , A32,1H#)') 'STR_SECTARGET', TRIM(STR_SECTARGET)
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'A_ST_POL', A_ST_POL
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'A_ST_AZIM_IN', A_ST_AZIM_IN
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'A_ST_AZIM_OUT', A_ST_AZIM_OUT
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'SA_ST_IN', SA_ST_IN
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'SA_ST_OUT', SA_ST_OUT
        WRITE (6,'(1H#,A16,3H = , A32,1H#)') 'STR_FILTER', TRIM(STR_FILTER)
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'D_FILTER', D_FILTER
        WRITE (6,'(1H#,A16,3H = , A32,1H#)') 'STR_SAMPLE', TRIM(STR_SAMPLE)
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'CONC', CONC
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'VTUBE', VTUBE
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'ITUBE', ITUBE
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'ESTEP', ESTEP
        WRITE (6,'(1H#,A16,3H = , I32,1H#)') 'NSTEP', NSTEP
        WRITE (6,'(1H#,A16,3H = , F32.4,1H#)') 'EMIN', EMIN
        WRITE (6,'(53("#"))')
    END SUBROUTINE DISPCFG
END MODULE CFGDATA
