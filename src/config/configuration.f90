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
    
    character(len=MAXFILELENGTH), parameter :: kDEFAULT_OPTIONS_FILENAME = "downscale_options.txt"
    character(len=MAXFILELENGTH), parameter :: kVERSION_STRING = "0.1"
    
    logical :: module_debug
    public :: read_config
    public :: read_files_list, read_data_type, get_options_file ! only need to be public for test_config
contains
    
    !>------------------------------------------------
    !! Finds the options filename to use and reads the configuration
    !!
    !!------------------------------------------------
    function read_config() result(options)
        implicit none
        type(config) :: options
        
        call read_base_options(options)
        
        options%prediction  = read_prediction_options(  options%prediction_file,    options%debug)
        options%training    = read_training_options(    options%training_file,      options%debug)
        options%obs         = read_obs_options(         options%observation_file,   options%debug)
    end function read_config


    subroutine read_base_options(options)
        implicit none
        type(config), intent(inout)     :: options
        
        integer :: name_unit, n_analogs, n_log_analogs
        character(len=MAXSTRINGLENGTH)  :: name, start_date, end_date, start_train, end_train
        character(len=MAXSTRINGLENGTH)  :: start_transform, end_transform
        character(len=MAXFILELENGTH)    :: training_file, prediction_file, observation_file, output_file
        logical :: pure_analog, analog_regression, pure_regression, debug
        logical :: sample_analog, logistic_from_analog_exceedance
        real    :: logistic_threshold, analog_threshold
        
        ! setup the namelist
        namelist /parameters/   name, debug,                                        &
                                training_file, prediction_file, observation_file,   &
                                output_file,                                        &
                                start_date, end_date, start_train, end_train,       &
                                start_transform, end_transform,                     &
                                n_analogs, n_log_analogs, logistic_threshold,       &
                                pure_analog, analog_regression, pure_regression,    &
                                sample_analog, logistic_from_analog_exceedance,     &
                                analog_threshold

        options%version = kVERSION_STRING
        options%options_filename = get_options_file()
        
        ! defaults
        training_file    = options%options_filename
        prediction_file  = options%options_filename
        observation_file = options%options_filename
        start_date       = ""
        end_date         = ""
        start_train      = ""
        end_train        = ""
        start_transform  = ""
        end_transform    = ""
        output_file      = "downscaled_output.nc"
        n_analogs        = -1
        n_log_analogs    = -1
        analog_threshold = -1
        pure_analog      = .False.
        analog_regression= .True.
        pure_regression  = .False.
        debug            = .True.
        logistic_threshold= kFILL_VALUE
        sample_analog    = .False.
        logistic_from_analog_exceedance = .False.
        
        options%name = options%options_filename
        
        ! read namelists
        open(io_newunit(name_unit), file=trim(options%options_filename))
        read(name_unit,nml=parameters)
        close(name_unit)
        
        ! this is the time to make predictions over
        call options%first_time%init("gregorian")
        call options%first_time%set(start_date)
        call options%last_time%init("gregorian")
        call options%last_time%set(end_date)
        
        ! this it the time period to use for calibration of the regression variables
        call options%training_start%init("gregorian")
        call options%training_start%set(start_train)
        call options%training_stop%init("gregorian")
        call options%training_stop%set(end_train)

        ! this is the time period to use when calculating e.g. quantile mapping transformations
        call options%transform_start%init("gregorian")
        call options%transform_start%set(start_transform)
        call options%transform_stop%init("gregorian")
        call options%transform_stop%set(end_transform)
        
        options%training_file       = training_file
        options%prediction_file     = prediction_file
        options%observation_file    = observation_file
        options%output_file         = output_file
        options%n_analogs           = n_analogs
        options%n_log_analogs       = n_log_analogs
        options%analog_threshold    = analog_threshold
        options%pure_analog         = pure_analog
        options%analog_regression   = analog_regression
        options%pure_regression     = pure_regression
        
        options%logistic_threshold  = logistic_threshold
        options%sample_analog       = sample_analog
        options%logistic_from_analog_exceedance  = logistic_from_analog_exceedance


        options%debug = debug
        module_debug = options%debug

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
        
        integer :: name_unit, i

        ! namelist variables to be read
        integer :: nfiles, nvars, calendar_start_year, selected_time, interpolation_method
        integer, dimension(MAX_NUMBER_TIMES) :: time_indices
        integer, dimension(MAX_NUMBER_VARS)  :: selected_level
        character(len=MAXSTRINGLENGTH)       :: name, data_type, calendar
        character(len=MAXVARLENGTH)          :: lat_name, lon_name, time_name
        character(len=MAXFILELENGTH)         :: preloaded
        character(len=MAXFILELENGTH), dimension(MAX_NUMBER_VARS) :: file_list
        character(len=MAXVARLENGTH),  dimension(MAX_NUMBER_VARS) :: var_names
        integer, dimension(MAX_NUMBER_VARS) :: input_transformations

        ! setup the namelist
        namelist /training_parameters/ nfiles, nvars, name, data_type,    &
                                         lat_name, lon_name, time_name,   &
                                         file_list, var_names,            &
                                         calendar, calendar_start_year,   &
                                         selected_time, time_indices,     &
                                         interpolation_method, preloaded, &
                                         selected_level, input_transformations
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
        selected_time = -1
        time_indices = -1
        interpolation_method = kNEAREST
        input_transformations = kNO_TRANSFORM
        preloaded   = ""
        selected_level = -1
        
        ! read namelists
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=training_parameters)
        close(name_unit)
        
        
        if (nfiles <= 0) stop "Number of training files (nfiles) is not valid. "
        if (nvars  <= 0) stop "Number of training variables (nvars) is not valid. "
        
        ! allocate necessary arrays
        allocate(training_options%file_names(nfiles,nvars))
        allocate(training_options%var_names(nvars))
        allocate(training_options%selected_level(nvars))
        allocate(training_options%input_Xforms(nvars))
        
        ! finally, store the data into the config structure
        training_options%name           = name
        training_options%n_variables    = nvars
        training_options%nfiles         = nfiles
        do i=1,nvars
            training_options%var_names(i)    = var_names(i)
            nfiles = read_files_list(file_list(i), training_options%file_names(:,i))
            if (nfiles /= training_options%nfiles) then
                write(*,*) "WARNING read a different number of training input files than specified"
                write(*,*) nfiles, training_options%nfiles, trim(var_names(i))
            endif
        end do
        training_options%lat_name       = lat_name
        training_options%lon_name       = lon_name
        training_options%time_name      = time_name
        training_options%calendar       = calendar
        training_options%calendar_start_year = calendar_start_year
        training_options%time_file      = 1
        training_options%selected_time  = selected_time
        call copy_array_i(time_indices,   training_options%time_indices)
        training_options%selected_level = selected_level(1:nvars)
        training_options%data_type      = read_data_type(data_type)
        training_options%interpolation_method = interpolation_method
        training_options%input_Xforms   = input_transformations(1:nvars)
        training_options%debug          = debug
        training_options%preloaded      = preloaded
        
        where(training_options%selected_level == -1) training_options%selected_level = 1
        call check_training_options(training_options)
    end function read_training_options
    
    subroutine copy_array_i(input, output, invalid)
        implicit none
        integer, dimension(:), intent(in)   :: input
        integer, dimension(:), allocatable, intent(inout):: output
        integer, optional :: invalid
        integer :: i,n, invalid_test
        
        if (present(invalid)) then
            invalid_test = invalid
        else
            invalid_test = -1
        endif
        
        !  first find the number of valid entries in the input
        n = 0
        do i = 1, size(input)
            if (input(i) /= invalid_test) n = n+1
        end do
        n = max(1, n)
        
        if (allocated(output)) then
            if (size(output)/=n) then
                deallocate(output)
                allocate(output(n))
            endif
        else
            allocate(output(n))
        endif
        
        output(:n)  = input(:n)
    end subroutine copy_array_i
    
    !>------------------------------------------------
    !!Verify the options read from the training options namelist
    !!
    !!------------------------------------------------
    subroutine check_training_options(opt)
        implicit none
        type(training_config) :: opt
        integer :: i, j
        
        if (trim(opt%lat_name) == "")  stop "Error : lat_name not supplied in training options. "
        if (trim(opt%lon_name) == "")  stop "Error : lon_name not supplied in training options. "
        if (trim(opt%calendar) == "")  stop "Error : Calendar not supplied in training options. "
        if (trim(opt%time_name) == "") stop "Error : time_name not supplied in training options. "
        
        
        do i = 1, opt%n_variables
            if (trim(opt%var_names(i)) == "") stop "Invalid or not enough variable names specified in training options. "
            do j = 1, opt%nfiles
                if (trim(opt%file_names(j,i)) == "") stop "Invalid or not enough file names specified in training options. "
            enddo
        end do

    end subroutine check_training_options

    
    !>------------------------------------------------
    !! Read the prediction configuration
    !!
    !!------------------------------------------------
    function read_prediction_options(filename, debug) result(prediction_options)
        implicit none
        character(len=*), intent(in) :: filename
        logical, intent(in)          :: debug
        type(prediction_config) :: prediction_options
        
        integer :: name_unit, i, j

        ! namelist variables to be read
        integer :: nfiles, nvars, calendar_start_year, selected_time, interpolation_method
        integer, dimension(MAX_NUMBER_TIMES) :: time_indices
        integer, dimension(MAX_NUMBER_VARS) :: selected_level
        character(len=MAXSTRINGLENGTH)  :: name, data_type, calendar
        character(len=MAXVARLENGTH)     :: lat_name, lon_name, time_name
        character(len=MAXFILELENGTH)    :: preloaded
        character(len=MAXFILELENGTH), dimension(MAX_NUMBER_VARS) :: file_list
        character(len=MAXVARLENGTH),  dimension(MAX_NUMBER_VARS) :: var_names
        integer, dimension(MAX_NUMBER_VARS) :: transformations, input_transformations

        ! setup the namelist
        namelist /prediction_parameters/ nfiles, nvars, name, data_type,    &
                                         lat_name, lon_name, time_name,     &
                                         file_list, var_names,              &
                                         calendar, calendar_start_year,     &
                                         selected_time,                     &
                                         input_transformations, transformations,&
                                         interpolation_method, preloaded,   &
                                         selected_level, time_indices

        !defaults :
        nfiles              = -1
        nvars               = -1
        name                = ""
        data_type           = ""
        lat_name            = ""
        lon_name            = ""
        time_name           = ""
        file_list           = ""
        var_names           = ""
        calendar            = ""
        calendar_start_year = 1900
        selected_time       = -1
        transformations     = kQUANTILE_MAPPING
        input_transformations = kNO_TRANSFORM
        interpolation_method = kNEAREST
        preloaded           = ""
        time_indices        = -1
        selected_level      = -1
        
        ! read namelists
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=prediction_parameters)
        close(name_unit)
        
        if (nfiles <= 0) stop "Number of prediction files (nfiles) is not valid. "
        if (nvars  <= 0) stop "Number of prediction variables (nvars) is not valid. "
        
        ! allocate necessary arrays
        allocate(prediction_options%file_names(nfiles,nvars))
        allocate(prediction_options%var_names(nvars))
        allocate(prediction_options%selected_level(nvars))
        allocate(prediction_options%transformations(nvars))
        allocate(prediction_options%input_Xforms(nvars))
        
        ! finally, store the data into the config structure
        prediction_options%name           = name
        prediction_options%n_variables    = nvars
        prediction_options%nfiles         = nfiles
        do i=1,nvars
            prediction_options%var_names(i)    = var_names(i)
            nfiles = read_files_list(file_list(i), prediction_options%file_names(:,i))
            if (nfiles/=prediction_options%nfiles) then
                write(*,*) nfiles, prediction_options%nfiles
                stop "Error reading the correct number of prediction input files"
            endif
        end do
        prediction_options%lat_name       = lat_name
        prediction_options%lon_name       = lon_name
        prediction_options%time_name      = time_name
        prediction_options%calendar       = calendar
        prediction_options%calendar_start_year = calendar_start_year
        prediction_options%selected_time  = selected_time
        prediction_options%selected_level = selected_level(1:nvars)
        call copy_array_i(time_indices, prediction_options%time_indices)
        prediction_options%time_file      = 1 
        prediction_options%data_type      = read_data_type(data_type)
        prediction_options%transformations= transformations
        prediction_options%input_Xforms   = input_transformations(1:nvars)
        prediction_options%interpolation_method = interpolation_method
        prediction_options%debug          = debug
        prediction_options%preloaded      = preloaded
        
        where(prediction_options%selected_level == -1) prediction_options%selected_level = 1
        call check_prediction_options(prediction_options)
        
    end function read_prediction_options
    
    !>------------------------------------------------
    !!Verify the options read from the prediction options namelist
    !!
    !!------------------------------------------------
    subroutine check_prediction_options(opt)
        implicit none
        type(prediction_config) :: opt
        integer :: i,j
        
        if (trim(opt%lat_name)  == "") stop "Error : lat_name not supplied in prediction options. "
        if (trim(opt%lon_name)  == "") stop "Error : lon_name not supplied in prediction options. "
        if (trim(opt%calendar)  == "") stop "Error : Calendar not supplied in prediction options. "
        if (trim(opt%time_name) == "") stop "Error : time_name not supplied in prediction options. "
        
        do i = 1, opt%n_variables
            if (trim(opt%var_names(i)) == "") stop "Invalid or not enough variable names specified in prediction options. "
            do j = 1, opt%nfiles
                if (trim(opt%file_names(j,i)) == "") stop "Invalid or not enough file names specified in prediction options. "
            enddo
        end do

    end subroutine check_prediction_options

    !>------------------------------------------------
    !! Read the obs configuration
    !!
    !!------------------------------------------------
    function read_obs_options(filename, debug) result(obs_options)
        implicit none
        character(len=*), intent(in) :: filename
        logical, intent(in)          :: debug
        type(obs_config) :: obs_options
        
        integer :: name_unit, i

        ! namelist variables to be read
        integer :: nfiles, nvars, calendar_start_year
        character(len=MAXSTRINGLENGTH)  :: name, data_type, calendar
        character(len=MAXVARLENGTH)     :: lat_name, lon_name, time_name
        character(len=MAXFILELENGTH)    :: preloaded
        character(len=MAXFILELENGTH), dimension(MAX_NUMBER_VARS) :: file_list
        character(len=MAXVARLENGTH),  dimension(MAX_NUMBER_VARS) :: var_names
        integer,                      dimension(MAX_NUMBER_VARS) :: input_transformations
        real :: mask_value, logistic_threshold
        integer :: mask_variable

        ! setup the namelist
        namelist /obs_parameters/ nfiles, nvars, name, data_type,           &
                                         lat_name, lon_name, time_name,     &
                                         file_list, var_names,              &
                                         calendar, calendar_start_year,     &
                                         mask_value, mask_variable,         &
                                         preloaded, logistic_threshold,     &
                                         input_transformations
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
        mask_value  = 1e20
        mask_variable = 1
        preloaded   = ""
        input_transformations = 0
        logistic_threshold = kFILL_VALUE
        
        ! read namelists
        open(io_newunit(name_unit), file=filename)
        read(name_unit,nml=obs_parameters)
        close(name_unit)
        
        if (nfiles <= 0) stop "Number of obs files (nfiles) is not valid. "
        if (nvars  <= 0) stop "Number of obs variables (nvars) is not valid. "
        
        ! allocate necessary arrays
        allocate(obs_options%file_names(nfiles,nvars))
        allocate(obs_options%var_names(nvars))
        allocate(obs_options%input_Xforms(nvars))
        
        ! finally, store the data into the config structure
        obs_options%name           = name
        obs_options%n_variables    = nvars
        obs_options%nfiles         = nfiles
        do i=1,nvars
            obs_options%var_names(i)    = var_names(i)
            nfiles = read_files_list(file_list(i), obs_options%file_names(:,i))
            if (nfiles/=obs_options%nfiles) then
                write(*,*) nfiles, obs_options%nfiles
                stop "Error reading the correct number of obs input files"
            endif
        end do
        obs_options%lat_name       = lat_name
        obs_options%lon_name       = lon_name
        obs_options%time_name      = time_name
        obs_options%calendar       = calendar
        obs_options%calendar_start_year = calendar_start_year
        obs_options%time_file      = 1 
        obs_options%data_type      = read_data_type(data_type)
        obs_options%mask_value     = mask_value
        obs_options%mask_variable  = mask_variable
        obs_options%debug          = debug
        obs_options%preloaded      = preloaded
        obs_options%logistic_threshold = logistic_threshold
        obs_options%input_Xforms   = input_transformations(1:nvars)
        
        call check_obs_options(obs_options)
    end function read_obs_options
    
    !>------------------------------------------------
    !!Verify the options read from the obs options namelist
    !!
    !!------------------------------------------------
    subroutine check_obs_options(opt)
        implicit none
        type(obs_config) :: opt
        integer :: i,j
        
        if (trim(opt%lat_name) == "")  stop "Error : lat_name not supplied in obs options. "
        if (trim(opt%lon_name) == "")  stop "Error : lon_name not supplied in obs options. "
        if (trim(opt%calendar) == "")  stop "Error : Calendar not supplied in obs options. "
        if (trim(opt%time_name) == "") stop "Error : time_name not supplied in obs options. "
        
        do i = 1, opt%n_variables
            if (trim(opt%var_names(i)) == "") stop "Invalid or not enough variable names specified in obs options. "
            do j = 1, opt%nfiles
                if (trim(opt%file_names(j,i)) == "") stop "Invalid or not enough file names specified in obs options. "
            enddo
        end do

    end subroutine check_obs_options

    
    !>-----------------------------------------------
    !! Convert a data type string into its corresponding constant integer expression
    !!
    !!  @param  type_name   [in]    character string describing a data type
    !!  @retval data_type   [out]   integer constant corresponding to the input string
    !!
    !!------------------------------------------------
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
            write(*,*) "ERROR: unknown data type: "//trim(type_name)
            write(*,*) "Must be one of: GEFS, GCM, obs"
            stop
        end select
        
    end function read_data_type

    !>-----------------------------------------------
    !! Read a long list of files (one per line) from an input file (filename)
    !!
    !!  @param  filename    [in]    character string containing the name of the file to read
    !!  @param  file_list   [out]   allocatable array of filenames read from the input file
    !!  @retval nfiles      [out]   number of files read in
    !!
    !!------------------------------------------------
    function read_files_list(filename, file_list) result(nfiles)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=MAXFILELENGTH), dimension(:), intent(inout) :: file_list
        
        integer :: nfiles
        
        character(len=MAXFILELENGTH), dimension(MAX_NUMBER_FILES)            :: forcing_files
        integer :: file_unit
        integer :: i, error
        character(len=MAXFILELENGTH) :: temporary_file
        
        if (module_debug) print*, "Reading: ",trim(filename)
        open(unit=io_newunit(file_unit), file=filename)
        i=0
        error=0
        temporary_file = ""
        do while (error==0)
            read(file_unit, *, iostat=error) temporary_file
            if (error==0) then
                i=i+1
                
                if (i > MAX_NUMBER_FILES) then
                    write(*,*) "ERROR reading: "//trim(filename)
                    stop "Too many files to read"
                endif
                
                forcing_files(i) = temporary_file
            endif
            
        enddo
        
        close(file_unit)
        
        nfiles = i
        
        file_list(1:nfiles) = forcing_files(1:nfiles)

    end function read_files_list

    !>------------------------------------------
    !! Get the first commandline option if one exists, otherwise return a default filename
    !!
    !!  @retval options_file    [out]   Filename to read the configuration options from
    !!
    !!------------------------------------------
    function get_options_file() result(options_file)
        implicit none
        character(len=MAXFILELENGTH) ::options_file
        integer :: error
        logical :: file_exists
    
        ! if a commandline argument was given
        if (command_argument_count()>0) then
            ! read the commandline argument
            call get_command_argument(1,options_file, status=error)
            
            ! if there was an error, return the default filename
            if (error>0) then
                options_file = kDEFAULT_OPTIONS_FILENAME
            ! if the error was -1, then the filename supplied was just too long
            elseif (error==-1) then
                write(*,*) "Options filename = ", trim(options_file), " ...<cutoff>"
                write(*,*) "Maximum filename length = ", MAXFILELENGTH
                stop "ERROR: options filename too long"
            endif
        else
            ! if not commandline options were given, assume a default filename
            options_file = kDEFAULT_OPTIONS_FILENAME
        endif
        
        ! check to see if the expected filename even exists on disk
        INQUIRE(file=trim(options_file), exist=file_exists)
        ! if it does not exist, print an error and stop
        if (.not.file_exists) then
            write(*,*) "Using options file = ", trim(options_file)
            stop "Options file does not exist. "
        endif
        write(*,*) "Using options file = ", trim(options_file)
    end function get_options_file
    
end module config_mod
