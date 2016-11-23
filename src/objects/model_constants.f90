module model_constants

    implicit none

    ! Model source code version number
    character(len=*),public,parameter :: kVERSION_STRING = "0.3"

    ! ------------------------------------------------
    ! String constants
    ! ------------------------------------------------
    integer,public,parameter :: MAXSTRINGLENGTH  = 1024     ! maximum random string length
    integer,public,parameter :: SM_STRING_LENGTH = 255      ! well defined small strings (unused at present)
    integer,public,parameter :: MAXFILELENGTH    = 1024     ! maximum file name length
    integer,public,parameter :: MAXVARLENGTH     = 1024     ! maximum variable name length
    integer,public,parameter :: MAXDIMLENGTH     = 1024     ! maximum variable name length

    ! ------------------------------------------------
    ! Size constants
    ! ------------------------------------------------
    integer,public,parameter :: MAX_NUMBER_TIMES = 48       ! maximum number of time steps to integrate from an input source
    integer,public,parameter :: MAX_NUMBER_VARS  = 255      ! maximum number of permitted variables to process
    integer,public,parameter :: MAX_NUMBER_FILES = 100000   ! maximum number of permitted input files
                                                            ! 100000 = 1 file/day for ~274 years
    integer,public,parameter :: N_ATM_QM_SEGMENTS = 300     ! Number of quantiles to use for quantile mapping atmospheric data

    ! ------------------------------------------------
    ! Input Type Constants
    ! ------------------------------------------------
    integer,public,parameter :: kGCM_TYPE         = 1
    integer,public,parameter :: kREANALYSIS_TYPE  = 2
    integer,public,parameter :: kGEFS_TYPE        = 3
    integer,public,parameter :: kOBS_TYPE         = 4


    ! ------------------------------------------------
    ! Data Transformation Type Constants
    ! ------------------------------------------------
    integer,public,parameter :: kNO_TRANSFORM     = 0
    integer,public,parameter :: kQUANTILE_MAPPING = 1
    integer,public,parameter :: kLOG_TRANSFORM    = 2
    integer,public,parameter :: kCUBE_ROOT        = 3
    integer,public,parameter :: kFIFTH_ROOT       = 4

    ! ------------------------------------------------
    ! Geospatial Interpolation Type Constants
    ! ------------------------------------------------
    integer,public,parameter :: kNEAREST          = 1
    integer,public,parameter :: kBILINEAR         = 2

    ! ------------------------------------------------
    ! Normalization Method
    ! ------------------------------------------------
    integer,public,parameter :: kPREDICTIONDATA   = 1
    integer,public,parameter :: kTRAININGDATA     = 2

    ! ------------------------------------------------
    ! Other Constants
    ! ------------------------------------------------
    integer,public,parameter :: kFILL_VALUE       = -9999
    character(len=*),public,parameter :: kDEFAULT_OPTIONS_FILENAME = "downscale_options.txt"



end module model_constants
