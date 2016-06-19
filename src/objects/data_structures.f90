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
        character(len=MAXVARLENGTH), allocatable, dimension(:) :: attributes_names
        character(len=MAXVARLENGTH), allocatable, dimension(:) :: attributes_values
    end type variable_type
    
    ! ------------------------------------------------
    ! adds a quantile mapping capability to atmospheric input
    ! ------------------------------------------------
    type, extends(variable_type) :: atm_variable_type
        type(qm_correction_type), allocatable, dimension(:,:,:) :: qm ! per gridpoint (per month? or DOY?) QM
        real, allocatable, dimension(:,:) :: mean, stddev ! per gridpoint mean and standard deviation (for normalization)
    end type atm_variable_type
    
    ! ------------------------------------------------
    ! adds mean and stddev statistics
    ! ------------------------------------------------
    type, extends(variable_type) :: obs_variable_type
        real, allocatable, dimension(:,:) :: mean, stddev ! per gridpoint mean and standard deviation (for normalization?)
        integer :: transformation                         ! type of transformation applied to data (e.g. sqrt, log, ???)
        real :: logistic_threshold                        ! threshold to use in logistic regression
    end type obs_variable_type

    type, extends(variable_type) :: output_variable_type
        real, allocatable, dimension(:,:)       :: mean, stddev ! per gridpoint mean and standard deviation (for normalization?)
        real, allocatable, dimension(:,:,:)     :: errors       ! store pre grid point expected errors from the downscaling code
        real, allocatable, dimension(:,:,:,:)   :: coefficients ! store pre grid point regression coefficients
        real, allocatable, dimension(:,:,:)     :: obs
        real, allocatable, dimension(:,:,:,:)   :: predictors
        real, allocatable, dimension(:,:,:,:)   :: training
        real, allocatable, dimension(:,:,:)     :: logistic
        real :: logistic_threshold
    end type output_variable_type
    
    type, extends(interpolable_type) :: base_data_type
        type(Time_type), allocatable, dimension(:) :: times
        integer :: n_variables, n_times
        character (len=MAXSTRINGLENGTH) :: name
        integer :: first_time, last_time
        integer :: training_start, training_stop
        integer :: transform_start, transform_stop
    end type base_data_type
    
    ! ------------------------------------------------
    ! type to contain atmospheric fields
    ! ------------------------------------------------
    type, extends(base_data_type) :: atm
        type(atm_variable_type), allocatable, dimension(:) :: variables
    end type atm
    
    ! ------------------------------------------------
    ! type to contain observation data
    ! ------------------------------------------------
    type, extends(base_data_type) :: obs
        type(obs_variable_type), allocatable, dimension(:) :: variables
        logical, dimension(:,:), allocatable :: mask
    end type obs
    
    type, extends(base_data_type) :: results
        type(output_variable_type), allocatable, dimension(:) :: variables
    end type results

    ! ------------------------------------------------
    ! types for the options for each sub-component (since these are identical for now, could we just use one...)
    ! ------------------------------------------------
    type input_config
        character (len=MAXSTRINGLENGTH) :: name
        
        character (len=MAXFILELENGTH), allocatable, dimension(:,:) :: file_names
        character (len=MAXVARLENGTH),  allocatable, dimension(:)   :: var_names
        integer :: n_variables, nfiles
        
        ! file prefix to read data from instead of a long list of file_names
        ! filenames are defined as <preloaded>_<variable_name>.nc will be read
        ! e.g. predictor_ua.nc
        character (len=MAXFILELENGTH) :: preloaded
        
        ! will try to read these from the input data files units attribute, but can be specified here if not
        character (len=MAXSTRINGLENGTH) :: calendar
        integer :: calendar_start_year          ! year to start the time data calendar (e.g. 1900-01-01)
        double precision :: time_gain      = 1  ! to convert file "time" data into days (from e.g. seconds) calculated from units
        
        integer :: data_type ! Type of data.  e.g. gcm, reanalysis, forecast
        
        integer :: time_file ! specify the variable number to use when selecting files to read time data from
        character (len=MAXVARLENGTH)    :: lat_name, lon_name, time_name
        
        integer :: selected_time   = -1         ! to just use a single time from each file (for e.g. GEFS forecast)
        
        logical :: debug
    end type input_config
    
    type, extends(input_config) :: atm_config
        integer, dimension(:), allocatable :: selected_level ! to just use a specific vertical level from each file (e.g. a pressure level)
                                                             ! the array is to provide one level for each variable
        integer :: interpolation_method
        
        integer, dimension(:), allocatable :: time_indices   ! specific time indices to average over (e.g. multiple hours in a daily file)
        ! tranformation to apply to each atmophseric variable
        integer, dimension(:), allocatable :: transformations
        integer, dimension(:), allocatable :: input_Xforms
    end type atm_config
    
    type, extends(atm_config) :: prediction_config
    end type prediction_config
    
    type, extends(atm_config) :: training_config
    end type training_config
    
    type, extends(input_config) :: obs_config
        real :: logistic_threshold
        real :: mask_value
        integer :: mask_variable
    end type obs_config
    
    
    ! ------------------------------------------------
    ! store all model options
    ! ------------------------------------------------
    type config
        character (len=MAXSTRINGLENGTH) :: version, comment

        character (len=MAXFILELENGTH) :: options_filename
        character (len=MAXFILELENGTH) :: training_file
        character (len=MAXFILELENGTH) :: prediction_file
        character (len=MAXFILELENGTH) :: observation_file
        
        character (len=MAXFILELENGTH) :: name
        ! file names
        character (len=MAXFILELENGTH) :: output_file
        
        logical :: pure_analog
        logical :: analog_regression
        logical :: pure_regression
        
        ! if not equal to kFILL_VALUE then it will be used to generate a probability of exceedance
        real    :: logistic_threshold
        
        ! options for each sub-component
        type(training_config)      :: training
        type(obs_config)           :: obs
        type(prediction_config)    :: prediction
        
        ! date/time parameters
        type(Time_type) :: training_start, training_stop    ! define the period over which the model should be trained
        type(Time_type) :: first_time,     last_time        ! define the period over which the model should be applied
        type(Time_type) :: transform_start, transform_stop  ! define the period over which any transformations should be developed
                                                            ! e.g. to Quantile map GCM data into training atm data space
                                                            
        integer :: first_point, last_point ! start and end positions to run the model for(?)
        integer :: n_analogs
        
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
