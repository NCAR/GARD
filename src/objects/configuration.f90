!>------------------------------------------------
!! Reads the configuration input file
!! Also looks for commandline arguments to for an config filename
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module config_mod
    
    use data_structures
    use model_constants
    
    implicit none
    private
    
    character(len=MAXFILELENGTH), parameter :: DEFAULT_OPTIONS_FILENAME = "downscale_options.txt"
    character(len=MAXFILELENGTH), parameter :: VERSION_STRING = "0.1"
    
    public :: read_config
contains
    
    !>------------------------------------------------
    !! Finds the options filename to use and reads the configuration
    !!
    !!------------------------------------------------
    function read_config() result(options)
        implicit none
        type(config) :: options
        character(len=MAXFILELENGTH) :: filename
        
        options%version = VERSION_STRING
        
        ! Look for commandline argument number 1
        CALL get_command_argument(1, filename)
        if (len_trim(filename) > 0) then
            options%options_filename = filename
        else
            ! if no commandline options were given, use the default filename
            options%options_filename = DEFAULT_OPTIONS_FILENAME
        endif
        
        options%name = options%options_filename
        options%debug = .True.
        
        options%training%name = "Training data"
        
        ! this is the time to make predictions over
        call options%first_time%init("gregorian")
        call options%first_time%set("2000-01-01 00:00:00")
        call options%last_time%init("gregorian")
        call options%last_time%set("2001-01-01 00:00:00")
        
        ! this it the time period to use for calibration of the regression variables
        call options%training_start%init("gregorian")
        call options%training_start%set("1990-01-01 00:00:00")
        call options%training_stop%init("gregorian")
        call options%training_stop%set("1991-01-01 00:00:00")

        call options%transform_start%init("gregorian")
        call options%transform_start%set("1990-01-01 00:00:00")
        call options%transform_stop%init("gregorian")
        call options%transform_stop%set("1991-01-01 00:00:00")
        
        
        options%training    = read_training_options(options)
        options%prediction  = read_prediction_options(options)
        options%obs         = read_obs_options(options)
        
    end function read_config


    !>------------------------------------------------
    !! Read the training configuration
    !!
    !!------------------------------------------------
    function read_training_options(options) result(training_options)
        implicit none
        type(config), intent(in) :: options
        type(training_config) :: training_options
        integer :: nfiles, nvars

        nfiles = 1
        nvars = 3
        
        allocate(training_options%file_names(nfiles,nvars))
        allocate(training_options%var_names(nvars))
        
        training_options%name           = "GEFS test"
        training_options%n_variables    = nvars
        training_options%file_names(1,1)= "/d5/gefs/all/2010/20101127/spfh_2m_2010112700_mean.nc"
        training_options%var_names(1)   = "SPFH_2maboveground"
        training_options%file_names(1,2)= "/d5/gefs/all/2010/20101127/ugrd_80m_2010112700_mean.nc"
        training_options%var_names(2)   = "UGRD_80maboveground"
        training_options%file_names(1,3)= "/d5/gefs/all/2010/20101127/vgrd_80m_2010112700_mean.nc"
        training_options%var_names(3)   = "VGRD_80maboveground"
        training_options%lat_name       = "latitude"
        training_options%lon_name       = "longitude"
        training_options%data_type      = kGEFS_TYPE
        training_options%debug          = options%debug
        training_options%calendar       = "standard"
        training_options%time_gain      = 1/86400.0D0
        training_options%calendar_start_year = 1970
        training_options%time_file      = 1 
        training_options%time_name      = "time"
        training_options%selected_time  = 1
        
    end function read_training_options



    !>------------------------------------------------
    !! Read the prediction configuration
    !!
    !!------------------------------------------------
    function read_prediction_options(options) result(prediction_options)
        implicit none
        type(config), intent(in) :: options
        type(prediction_config) :: prediction_options
        integer :: nfiles, nvars

        nfiles = 1
        nvars = 3
        
        allocate(prediction_options%file_names(nfiles,nvars))
        allocate(prediction_options%var_names(nvars))
        
        prediction_options%name           = "GCM test"
        prediction_options%n_variables    = nvars
        prediction_options%file_names(1,1)= "/d4/gutmann/cmip/daily/ccsm/subset/hus_day_CCSM4_historical_r6i1p1_19750101-19791231.nc"
        prediction_options%var_names(1)   = "hus"
        prediction_options%file_names(1,2)= "/d4/gutmann/cmip/daily/ccsm/subset/ua_day_CCSM4_historical_r6i1p1_19750101-19791231.nc"
        prediction_options%var_names(2)   = "ua"
        prediction_options%file_names(1,3)= "/d4/gutmann/cmip/daily/ccsm/subset/va_day_CCSM4_historical_r6i1p1_19750101-19791231.nc"
        prediction_options%var_names(3)   = "va"
        prediction_options%lat_name       = "lat"
        prediction_options%lon_name       = "lon"
        prediction_options%data_type      = kGCM_TYPE
        prediction_options%debug          = options%debug
        prediction_options%calendar       = "noleap"
        prediction_options%calendar_start_year = 1850
        prediction_options%time_file      = 1 
        prediction_options%time_name      = "time"
        
    end function read_prediction_options

    !>------------------------------------------------
    !! Read the training configuration
    !!
    !!------------------------------------------------
    function read_obs_options(options) result(obs_options)
        implicit none
        type(config), intent(in) :: options
        type(obs_config) :: obs_options
        integer :: nfiles, nvars

        nfiles = 2
        nvars = 2
        
        allocate(obs_options%file_names(nfiles,nvars))
        allocate(obs_options%var_names(nvars))
        
        obs_options%name           = "OBS test"
        obs_options%n_variables    = nvars
        obs_options%file_names(1,1)= "/d2/gutmann/usbr/stat_data/DAILY/obs/maurer.125/pr/nldas_met_update.obs.daily.pr.1979.nc"
        obs_options%file_names(2,1)= "/d2/gutmann/usbr/stat_data/DAILY/obs/maurer.125/pr/nldas_met_update.obs.daily.pr.1980.nc"
        obs_options%var_names(1)   = "pr"
        obs_options%file_names(1,2)= "/d2/gutmann/usbr/stat_data/DAILY/obs/maurer.125/tasmax/nldas_met_update.obs.daily.tasmax.1979.nc"
        obs_options%file_names(2,2)= "/d2/gutmann/usbr/stat_data/DAILY/obs/maurer.125/tasmax/nldas_met_update.obs.daily.tasmax.1980.nc"
        obs_options%var_names(2)   = "tasmax"
        obs_options%lat_name       = "lat"
        obs_options%lon_name       = "lon"
        obs_options%data_type      = kOBS_TYPE
        obs_options%debug          = options%debug
        obs_options%calendar       = "standard"
        obs_options%calendar_start_year = 1940
        obs_options%time_file      = 1 
        obs_options%time_name      = "time"
        
    end function read_obs_options
    
end module config_mod
