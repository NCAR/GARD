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
    use string,      only : str
    use io_routines, only : io_newunit
    
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
        
        call read_base_options(options)
        
        options%training    = read_training_options( options%training_file,      options%debug)
        options%prediction  = read_prediction_options( options%prediction_file,  options%debug)
        options%obs         = read_obs_options( options%observation_file,        options%debug)
        
    end function read_config

    subroutine read_base_options(options)
        implicit none
        type(config), intent(inout) :: options
        
        options%version = VERSION_STRING
        options%options_filename = get_options_file()
        
        options%name = options%options_filename
        options%debug = .True.
        
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
        
        options%training_file    = options%options_filename
        options%prediction_file  = options%options_filename
        options%observation_file = options%options_filename
        
    end subroutine read_base_options


    !>------------------------------------------------
    !! Read the training configuration
    !!
    !!------------------------------------------------
    function read_training_options(filename, debug) result(training_options)
        implicit none
        character(len=*), intent(in) :: filename
        logical, intent(in)          :: debug
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
        training_options%debug          = debug
        training_options%calendar       = "standard"
        training_options%time_gain      = 1/86400.0D0
        training_options%calendar_start_year = 1970
        training_options%time_file      = 1 
        training_options%time_name      = "time"
        training_options%selected_time  = 1
        
    end function read_training_options


    function get_options_file() result(options_file)
        implicit none
        character(len=MAXFILELENGTH) ::options_file
        integer :: error
        logical :: file_exists
    
        if (command_argument_count()>0) then
            call get_command_argument(1,options_file, status=error)
            if (error>0) then
                options_file = DEFAULT_OPTIONS_FILENAME
            elseif (error==-1) then
                write(*,*) "Options filename = ", trim(options_file), " ...<cutoff>"
                write(*,*) "Maximum filename length = ", MAXFILELENGTH
                stop "ERROR: options filename too long"
            endif
        else
            options_file = DEFAULT_OPTIONS_FILENAME
        endif
        INQUIRE(file=trim(options_file), exist=file_exists)
        if (.not.file_exists) then
            write(*,*) "Using options file = ", trim(options_file)
            stop "Options file does not exist. "
        endif
    end function

    !>------------------------------------------------
    !! Read the prediction configuration
    !!
    !!------------------------------------------------
    function read_prediction_options(filename, debug) result(prediction_options)
        implicit none
        character(len=*), intent(in) :: filename
        logical, intent(in)          :: debug
        type(prediction_config) :: prediction_options
        
        integer :: name_unit, i

        ! namelist variables to be read
        integer :: nfiles, nvars, calendar_start_year
        character(len=MAXSTRINGLENGTH)  :: name, data_type, calendar
        character(len=MAXVARLENGTH)     :: lat_name, lon_name, time_name
        character(len=MAXFILELENGTH), dimension(:), allocatable :: file_list
        character(len=MAXVARLENGTH),  dimension(:), allocatable :: var_names

        ! setup the namelist
        namelist /prediction_parameters/ nfiles, nvars, name, data_type,   &
                                      lat_name, lon_name, time_name,    &
                                      file_list, var_names
        !defaults :
        nfiles      = -1
        nvars       = -1
        name        = ""
        data_type   = ""
        lat_name    = ""
        lon_name    = ""
        time_name   = ""
        file_list   = ""
        var_names   = ""
        calendar    = ""
        calendar_start_year = 1900
        
        ! read namelists
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=prediction_parameters)
        close(name_unit)
        
        ! allocate necessary arrays
        allocate(prediction_options%file_names(nfiles,nvars))
        allocate(prediction_options%var_names(nvars))
        
        ! finally, store the data into the config structure
        prediction_options%name           = name
        prediction_options%n_variables    = nvars
        ! prediction_options%file_names(1,1)= "/d4/gutmann/cmip/daily/ccsm/subset/hus_day_CCSM4_historical_r6i1p1_19750101-19791231.nc"
        do i=1,nvars
            prediction_options%var_names(i)    = var_names(i)
            prediction_options%file_names(:,i) = read_files_list(file_list(i))
        end do
        prediction_options%lat_name       = lat_name
        prediction_options%lon_name       = lon_name
        prediction_options%time_name      = time_name
        prediction_options%calendar       = calendar
        prediction_options%calendar_start_year = calendar_start_year
        prediction_options%time_file      = 1 
        prediction_options%data_type      = read_data_type(data_type)
        prediction_options%debug          = debug
        
    end function read_prediction_options

    !>------------------------------------------------
    !! Read the training configuration
    !!
    !!------------------------------------------------
    function read_obs_options(filename, debug) result(obs_options)
        implicit none
        character(len=*), intent(in) :: filename
        logical, intent(in)          :: debug
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
        obs_options%debug          = debug
        obs_options%calendar       = "standard"
        obs_options%calendar_start_year = 1940
        obs_options%time_file      = 1 
        obs_options%time_name      = "time"
        
    end function read_obs_options
    
    function read_data_type(type_name) result(data_type)
        implicit none
        character(len=*) :: type_name
        integer :: data_type
        
        select case(trim(type_name))
        case("GEFS")
            data_type = kGEFS_TYPE
        case("GCM")
            data_type = kGCM_TYPE
        case("obs")
            data_type = kOBS_TYPE
        case default
            print*, "ERROR: unknown data type: "//trim(type_name)
            print*, "Must be one of: GEFS, GCM, obs"
            stop
        end select
        
    end function read_data_type

    function read_files_list(filename) result(file_list)
        implicit none
        character(len=*) :: filename
        
        character(len=MAXFILELENGTH), dimension(MAX_NUMBER_FILES) :: forcing_files
        character(len=MAXFILELENGTH), dimension(:), allocatable   :: file_list
        integer :: nfiles
        integer :: file_unit
        integer :: i, error
        character(len=MAXFILELENGTH) :: temporary_file
        
        open(unit=io_newunit(file_unit), file=filename)
        i=0
        error=0
        do while (error==0)
            read(file_unit, *, iostat=error) temporary_file
            if (error==0) then
                if ((i+1) > MAX_NUMBER_FILES) then
                    write(*,*) "ERROR reading: "//trim(filename)
                    stop "Too many files to read"
                endif
                i=i+1
                forcing_files(i) = temporary_file
            endif
        enddo
        close(file_unit)
        
        nfiles = i
        ! print out a summary
        write(*,*) "Boundary conditions files to be used:"
        if (nfiles>10) then
            write(*,*) "  nfiles=", trim(str(nfiles)), ", too many to print."
            write(*,*) "  First file:", trim(forcing_files(1))
            write(*,*) "  Last file: ", trim(forcing_files(nfiles))
        else
            do i=1,nfiles
                write(*,*) "    ",trim(forcing_files(i))
            enddo
        endif
        
        allocate(file_list(nfiles))
        file_list(1:nfiles) = forcing_files(1:nfiles)

    end function read_files_list
    
end module config_mod
