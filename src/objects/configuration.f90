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
        
        options%training = read_training_options(options)
        
    end function read_config


    function read_training_options(options) result(training_options)
        implicit none
        type(config), intent(in) :: options
        type(training_config) :: training_options
        integer :: nfiles, nvars

        nfiles = 1
        nvars = 1
        
        allocate(training_options%file_names(nfiles,nvars))
        allocate(training_options%var_names(nvars))
        
        training_options%name           = "GCM test"
        training_options%n_variables    = nvars
        training_options%file_names(1,1)= "ua_short.nc"
        training_options%var_names(1)   = "ua"
        training_options%data_type      = kGCM_TYPE
        training_options%debug          = options%debug
        
        ! type(Time_type),               allocatable, dimension(:,:) :: file_start, file_end
        ! character (len=MAXVARLENGTH) :: lat_name, lon_name, time_name
    end function read_training_options
    
end module config_mod
