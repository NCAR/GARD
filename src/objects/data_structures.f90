!>------------------------------------------------
!! Contains type definitions for a variety of model data strucutres
!! Also defines model constants (e.g. gravity, and MAXFILELENGTH)
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module data_structures
    use model_constants
    use time
    implicit none
    
    ! private
    ! public, config, atm,  obs,  results
    
    
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
    
    ! ------------------------------------------------
    ! a geographic look up table for spatial interpolation, from x,y with weight w
    ! ------------------------------------------------
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
    
    ! ------------------------------------------------
    ! type to store quantile mapping data
    ! ------------------------------------------------
    type qm_correction_type
        real, allocatable, dimension(:) :: start_idx, end_idx
        real, allocatable, dimension(:) :: slope, offset
    end type qm_correction_type
    
    ! ------------------------------------------------
    ! type to contain a single variable (atm or obs)
    ! ------------------------------------------------
    type variable_type
        character(len=MAXVARLENGTH)         :: name      ! name of the variable
        real, allocatable, dimension(:,:,:) :: data      ! raw data
        integer                             :: data_type ! Type of data.  e.g. precip, temperature, or other. 
    end type variable_type
    
    ! ------------------------------------------------
    ! adds a quantile mapping capability to atmospheric input
    ! ------------------------------------------------
    type, extends(variable_type) :: atm_variable_type
        type(qm_correction_type), allocatable, dimension(:,:,:) :: qm ! per gridpoint (per month? or DOY?) QM
        real, allocatable, dimension(:,:) :: mean, stddev ! per gridpoint mean and standard deviation (for normalization)
    end type atm_variable_type

    ! ------------------------------------------------
    ! type to contain atmospheric fields
    ! ------------------------------------------------
    type, extends(interpolable_type) :: atm
        type(atm_variable_type), allocatable, dimension(:) :: variables
        type(Time_type), allocatable, dimension(:) :: times
        integer :: n_variables, n_times
        character (len=MAXSTRINGLENGTH) :: name
    end type atm
    
    ! ------------------------------------------------
    ! type to contain observation data
    ! ------------------------------------------------
    type, extends(interpolable_type) :: obs
        type(variable_type), allocatable, dimension(:) :: variables
        integer :: n_variables
    end type obs
    

    ! ------------------------------------------------
    ! types for the options for each sub-component (since these are identical for now, could we just use one...)
    ! ------------------------------------------------
    type input_config
        character (len=MAXFILELENGTH), allocatable, dimension(:,:) :: file_names
        character (len=MAXVARLENGTH),  allocatable, dimension(:)   :: var_names
        integer :: n_variables
        
        type(Time_type),               allocatable, dimension(:,:) :: file_start, file_end
        character (len=MAXSTRINGLENGTH) :: calendar
        integer :: calendar_start_year
        
        integer :: data_type ! Type of data.  e.g. gcm, reanalysis, forecast
        integer :: time_file
        character (len=MAXVARLENGTH) :: lat_name, lon_name, time_name
        character (len=MAXSTRINGLENGTH) :: name
        logical :: debug
    end type input_config
    
    type, extends(input_config) :: prediction_config
    end type prediction_config
    
    type, extends(input_config) :: obs_config
    end type obs_config
    
    type, extends(input_config) :: training_config
    end type training_config
    
    ! ------------------------------------------------
    ! store all model options
    ! ------------------------------------------------
    type config
        character (len=MAXSTRINGLENGTH) :: version, comment

        character (len=MAXFILELENGTH) :: options_filename
        character (len=MAXFILELENGTH) :: name
        ! file names
        character (len=MAXFILELENGTH) :: output_file
        
        ! options for each sub-component
        type(training_config)      :: training
        type(obs_config)           :: obs
        type(prediction_config)    :: prediction
        
        ! date/time parameters
        type(Time_type) :: training_start, training_stop  ! define the period over which the model should be trained
        type(Time_type) :: first_time,     last_time      ! define the period over which the model should be applied
        
        integer :: first_point, last_point ! start and end positions to run the model for
        
        logical :: debug
        integer :: warning_level        ! level of warnings to issue when checking options settings 0-10.  
                                        ! 0  = Don't print anything
                                        ! 1  = print serious warnings
        ! (DEFAULT if debug=True)       ! 2  = print all warnings
                                        ! 3-4 ... nothing specified equivalent to 2
        ! (DEFAULT if debug=False)      ! 5  = Stop for options that are likely to break the model (print all warnings) 
                                        ! 6-8... nothing specified equivalent to 5
                                        ! 9  = stop on serious warnings only
                                        ! 10 = stop on all warnings
    end type config
end module data_structures  
