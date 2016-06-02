module model_constants

    implicit none
    
    ! ------------------------------------------------
    ! Model constants (string lengths)
    ! ------------------------------------------------
    integer,public,parameter :: MAXSTRINGLENGTH  = 1024     ! maximum random string length
    integer,public,parameter :: SM_STRING_LENGTH = 255      ! well defined small strings (unused at present)
    integer,public,parameter :: MAXFILELENGTH    = 1024     ! maximum file name length
    integer,public,parameter :: MAXVARLENGTH     = 1024     ! maximum variable name length
    integer,public,parameter :: MAXDIMLENGTH     = 1024     ! maximum variable name length
    integer,public,parameter :: MAX_NUMBER_VARS  = 255      ! maximum number of permitted variables to process 
    integer,public,parameter :: MAX_NUMBER_FILES = 100000   ! maximum number of permitted input files 
                                                            ! 100000 = 1 file/day for ~274 years

    ! ------------------------------------------------
    ! Input Type Constants
    ! ------------------------------------------------
    integer,public,parameter :: kGCM_TYPE         = 1
    integer,public,parameter :: kREANALYSIS_TYPE  = 2
    integer,public,parameter :: kGEFS_TYPE        = 3
    integer,public,parameter :: kOBS_TYPE         = 4
    
    ! ------------------------------------------------
    ! Physical Constants
    ! ------------------------------------------------
    ! real,public, parameter :: LH_vaporization=2260000.0 ! J/kg
    ! ! could be calculated as 2.5E6 + (-2112.0)*temp_degC ?
    ! real,public, parameter :: Rd  = 287.058   ! J/(kg K) specific gas constant for dry air
    ! real,public, parameter :: Rw  = 461.5     ! J/(kg K) specific gas constant for moist air
    ! real,public, parameter :: cp  = 1012.0    ! J/kg/K   specific heat capacity of moist STP air? 
    ! real,public, parameter :: gravity= 9.81   ! m/s^2    gravity
    ! real,public, parameter :: pi  = 3.1415927 ! pi
    ! real,public, parameter :: stefan_boltzmann = 5.67e-8 ! the Stefan-Boltzmann constant
    ! real,public, parameter :: karman = 0.41   ! the von Karman constant


end module model_constants
