!>------------------------------------------------
!! Contains type definitions for a variety of model data strucutres
!! Also defines model constants (e.g. gravity, and MAXFILELENGTH)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module data_structures
    use time
    implicit none

! ------------------------------------------------
! Model constants (string lengths)
! ------------------------------------------------
    integer,parameter::MAXSTRINGLENGTH  = 1024  ! maximum random string length
    integer,parameter::MAXFILELENGTH    = 1024  ! maximum file name length
    integer,parameter::MAXVARLENGTH     = 1024  ! maximum variable name length
    integer,parameter::MAX_NUMBER_FILES = 50000 ! maximum number of permitted input files (probably a bit extreme)

! ------------------------------------------------
! Physical Constants
! ------------------------------------------------
    real, parameter :: LH_vaporization=2260000.0 ! J/kg
    ! could be calculated as 2.5E6 + (-2112.0)*temp_degC ?
    real, parameter :: Rd  = 287.058   ! J/(kg K) specific gas constant for dry air
    real, parameter :: Rw  = 461.5     ! J/(kg K) specific gas constant for moist air
    real, parameter :: cp  = 1012.0    ! J/kg/K   specific heat capacity of moist STP air? 
    real, parameter :: gravity= 9.81   ! m/s^2    gravity
    real, parameter :: pi  = 3.1415927 ! pi
    real, parameter :: stefan_boltzmann = 5.67e-8 ! the Stefan-Boltzmann constant
    real, parameter :: karman = 0.41   ! the von Karman constant
    
! ------------------------------------------------
!   various data structures for use in geographic interpolation routines
! ------------------------------------------------
    ! contains the location of a specific grid point
    type position
        integer::x,y
    end type position
    ! contains location of surrounding 4 grid cells
    type fourpos
        integer::x(4),y(4)
    end type fourpos
    
    ! a geographic look up table for spatial interpolation, from x,y with weight w
    type geo_look_up_table
        ! x,y index positions, [n by m by 4] where there are 4 surrounding low-res points 
        ! for every high resolution point grid point to interpolate to
        integer,allocatable,dimension(:,:,:)::x,y
        ! weights to use for each of the 4 surrounding gridpoints.  Sum(over axis 3) must be 1.0
        real,allocatable,dimension(:,:,:)::w
    end type geo_look_up_table

    ! ------------------------------------------------
    ! A look up table for vertical interpolation. from z with weight w
    ! ------------------------------------------------
    type vert_look_up_table
        ! z index positions for all x,y,z points (x 2 for above and below z levels)
        integer,allocatable,dimension(:,:,:,:)::z
        ! weights to use for each of the two surrounding points.  Sum (over axis 1) must be 1.0
        real,allocatable,dimension(:,:,:,:)::w
    end type vert_look_up_table

    ! ------------------------------------------------
    ! generic interpolable type so geo interpolation routines will work on winds, domain, or boundary conditions. 
    ! ------------------------------------------------
    type interpolable_type
        ! all interpolables must have position (lat, lon, z)
        real, allocatable, dimension(:,:) :: lat,lon
        real, allocatable, dimension(:,:,:) :: z
        
        ! these are the look up tables that describe how to interpolate vertically (vert_lut) and horizontally (geolut)
        type(vert_look_up_table)::vert_lut
        type(geo_look_up_table)::geolut
            
        ! used to keep track of whether or not a particular error has been printed yet for this structure
        logical :: dx_errors_printed=.False.
        logical :: dy_errors_printed=.False.
    end type interpolable_type

    type variable_type
        character(len=MAXVARLENGTH)         :: name ! name of the variable
        real, allocatable, dimension(:,:,:) :: data ! raw data
        ! need to store transformation information (e.g. quantile mapping?) here?
        ! or maybe other information (e.g. continuous? positive? zero-filled ala precip)
    end type variable_type

    ! ------------------------------------------------
    ! type to contain external wind fields, only real addition is nfiles... maybe this could be folded in elsewhere?
    ! ------------------------------------------------
    type, extends(interpolable_type) :: atm_type
        type(variable_type), allocatable, dimension(:) :: variables
        integer :: n_variables
    end type atm_type

    ! ------------------------------------------------
    ! types for the options for each sub-component (since these are identical for now, could we just use one...)
    ! ------------------------------------------------
    type gcm_options
        character (len=MAXFILELENGTH), allocatable, dimension(:)   :: file_names
        type(Time_type),               allocatable, dimension(:,:) :: file_start, file_end
        character (len=MAXVARLENGTH),  allocatable, dimension(:)   :: var_names
        character (len=MAXVARLENGTH) :: lat_name, lon_name, time_name, name
    end type gcm_options

    type obs_options
        character (len=MAXFILELENGTH), allocatable, dimension(:)   :: file_names
        type(Time_type),               allocatable, dimension(:,:) :: file_start, file_end
        character (len=MAXVARLENGTH),  allocatable, dimension(:)   :: var_names
        character (len=MAXVARLENGTH) :: lat_name, lon_name, time_name, name
    end type obs_options
    
    type training_options
        character (len=MAXFILELENGTH), allocatable, dimension(:)   :: file_names
        type(Time_type),               allocatable, dimension(:,:) :: file_start, file_end
        character (len=MAXVARLENGTH),  allocatable, dimension(:)   :: var_names
        character (len=MAXVARLENGTH) :: lat_name, lon_name, time_name, name
    end type training_options


    ! ------------------------------------------------
    ! store all model options
    ! ------------------------------------------------
    type options_type
        character (len=MAXVARLENGTH) :: version,comment

        ! file names
        character (len=MAXFILELENGTH) :: output_file
        
        ! options for each sub-component
        type(gcm_options)       :: gcm
        type(obs_options)       :: obs
        type(training_options)  :: training
        
        ! date/time parameters
        type(Time_type) :: training_start, training_stop  ! define the period over which the model should be trained
        type(Time_type) :: first_time,     last_time      ! define the period over which the model should be applied
        
        integer :: first_point, last_point ! start and end positions to run the model for
        
        integer :: warning_level        ! level of warnings to issue when checking options settings 0-10.  
                                        ! 0  = Don't print anything
                                        ! 1  = print serious warnings
        ! (DEFAULT if debug=True)       ! 2  = print all warnings
                                        ! 3-4 ... nothing specified equivalent to 2
        ! (DEFAULT if debug=False)      ! 5  = Stop for options that are likely to break the model (print all warnings) 
                                        ! 6-8... nothing specified equivalent to 5
                                        ! 9  = stop on serious warnings only
                                        ! 10 = stop on all warnings
    end type options_type
end module data_structures  
