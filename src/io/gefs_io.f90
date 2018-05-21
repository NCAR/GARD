!>------------------------------------------------
!! Handle all IO for GEFS data
!!
!! Loops through GEFS variables reading data and computing basic statistics
!! Reads time and lat/lon variables from the GEFS files as well.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
submodule(gefs_mod) gefs_io_implementation

    use model_constants
    use basic_stats_mod,only: time_mean, time_stddev, time_minval
    use string,         only: str
    use io_routines,    only: io_read, io_getdims, io_maxDims, file_exists
    use time_io,        only: read_times
    use geo,            only: standardize_coordinates

    implicit none
    logical :: debug = .False.

contains

    !>------------------------------------------------
    !! Initialize the GEFS module
    !!
    !!------------------------------------------------
    module subroutine init_GEFS_io(options)
        implicit none
        type(config), intent(in) :: options

        debug = options%debug

    end subroutine init_GEFS_io

    !>------------------------------------------------
    !! Read in the GEFS data
    !!
    !! Uses variable names and filenames defined in the options data structure
    !!
    !! Loops through all input variables, then reads lat, lon, and time variables
    !!
    !!------------------------------------------------
    module function read_GEFS(options) result(GEFS_data)
        implicit none
        class(atm_config), intent(in) :: options
        type(atm) :: GEFS_data

        integer :: var_idx, ntimesteps
        integer :: nx,ny

        ! allocate space to store all of the variables to be read
        allocate( GEFS_data%variables( options%n_variables ))
        GEFS_data%name = options%name

        ! loop over variables reading them in
        do var_idx = 1, options%n_variables
            associate(var => GEFS_data%variables(var_idx))

                var = read_GEFS_variable(options%var_names(var_idx),        &
                                         options%file_names(:, var_idx),    &
                                         options%selected_time,             &
                                         options%time_indices,              &
                                         options%time_weights,              &
                                         options%agg_method(var_idx),       &
                                         options%preloaded)

                if (var_idx==1) then
                    ntimesteps = size(var%data, 1)
                else
                    if (ntimesteps /= size(var%data, 1)) then
                        write(*,*) "GEFS 1st variable ntime:", trim(str( ntimesteps )),  &
                                   "current ntime:", trim(str( size( var%data, 1) ))
                        write(*,*) " For variable:",trim(options%var_names( var_idx ))
                        stop "Error: Number of timesteps are not constant between variables"
                    endif
                endif
                ! then compute grid mean and stddev statistics for re-normalization
                ! could also add a standard normal transformation?
                call compute_grid_stats(var)
            end associate
        enddo

        allocate(GEFS_data%times(ntimesteps))
        if (debug) write(*,*) "Reading Time data"
        call read_times(options, GEFS_data%times, options%timezone_offset)
        if (debug) then
            write(*,*) "Times cover the period:"
            write(*,*) "   "//trim(GEFS_data%times(1)%as_string())
            write(*,*) "   "//trim(GEFS_data%times(ntimesteps)%as_string())
            write(*,*) ""
            write(*,*) "Reading Lat / Lon coordinates"
        endif
        call io_read(options%file_names(1, 1), options%lat_name, GEFS_data%lat)
        call io_read(options%file_names(1, 1), options%lon_name, GEFS_data%lon)

        call standardize_coordinates(GEFS_data)

    end function read_GEFS

    subroutine compute_grid_stats(var)
        implicit none
        type(atm_variable_type), intent(inout) :: var

        integer :: nx, ny

        nx = size(var%data,2)
        ny = size(var%data,3)

        if (allocated(var%mean))    deallocate(var%mean)
        if (allocated(var%stddev))  deallocate(var%stddev)
        if (allocated(var%min_val))  deallocate(var%min_val)
        allocate(var%mean(nx,ny))
        allocate(var%stddev(nx,ny))
        allocate(var%min_val(nx,ny))

        where(var%data>1e10) var%data=0

        call time_mean( var%data, var%mean )
        call time_stddev( var%data, var%stddev, mean_in=var%mean )
        call time_minval( var%data, var%min_val)

    end subroutine compute_grid_stats

    function read_GEFS_variable(varname, filenames, timestep, time_indices, time_weights, agg_method, preload) result(output)
        implicit none
        character(len=MAXVARLENGTH),    intent(in)              :: varname
        character(len=MAXFILELENGTH),   intent(in), dimension(:):: filenames
        integer,                        intent(in)              :: timestep
        integer,                        intent(in), dimension(:):: time_indices
        real,                           intent(in), dimension(:):: time_weights
        integer,                        intent(in)              :: agg_method
        character(len=MAXFILELENGTH),   intent(in), optional    :: preload

        type(atm_variable_type) :: output

        integer, dimension(io_maxDims) :: dims
        logical :: using_all_times

        output%name = varname

        if (present(preload)) then
            if (trim(preload) /= "" ) then
                if (file_exists(trim(preload)//trim(varname)//".nc")) then
                    if (debug) write(*,*) "Reading preloaded data: ", trim(preload)//trim(varname)//".nc"
                    call io_read(trim(preload)//trim(varname)//".nc", "data", output%data)
                    return
                endif
            endif
        endif

        ! if not subsets in by timestep or time_indices to average were specified then we are using all times in the file
        using_all_times = ((timestep == -1).and.(time_indices(1) == -1))
        dims = get_dims(varname, filenames, using_all_times)

        ! note, we reverse the order of the dimensions here to speed up later computations which will occur per grid cell over time
        allocate(output%data(dims(1), dims(2), dims(3)))

        call load_data(varname, filenames, output%data, timestep, time_indices, time_weights, agg_method)

    end function read_GEFS_variable

    !! requires all filenames to have the same number of time steps... no good for monthly files...
    function get_dims(varname, filenames, using_all_times) result(dims)
        implicit none
        character(len=MAXVARLENGTH),    intent(in)               :: varname
        character(len=MAXFILELENGTH),   intent(in), dimension(:) :: filenames
        logical,                        intent(in)               :: using_all_times

        integer :: file_idx, ntimesteps, nx, ny
        integer, dimension(io_maxDims) :: dims

        ntimesteps = 0

        ! make a passe through all the files to first find total the dimensions of the output data
        do file_idx = 1,size(filenames)

            call io_getdims(filenames(file_idx), varname, dims)
            ! dims(1) = the number of dimensions
            if (dims(1)<3) then
                ! if there are only two dimensions in the variable, assume there is only one time slice per file
                ! this should probably be checked somewhere
                ntimesteps = ntimesteps + 1
            else
                ! the last dimension is assumed to be time (+1 because dims(1)=ndims takes a slot)
                if (using_all_times) then
                    ntimesteps = ntimesteps + dims( dims(1)+1 )
                else
                    ntimesteps = ntimesteps + 1
                endif
            endif

            if (file_idx == 1) then
                nx = dims(2)
                ny = dims(3)
            else
                if ((nx/=dims(2)).or.(ny/=dims(3))) then
                    write(*,*) "GEFS 1st nx:", trim(str(nx)), "current nx:", trim(str(dims(2)))
                    write(*,*) " For file:",trim(filenames(file_idx)), " variable: ", trim(varname)
                    stop "Error: Grid dimensions are not constant over time"
                endif
            endif
        end do

        if (debug) then
            write(*,*) trim(str(ntimesteps)), " timesteps found in ",trim(str(size(filenames)))," files"
            write(*,*) "  The first file is: ", trim(filenames(1))
            write(*,*) "  For variable: ", trim(varname)
        endif

        ! only deals with 2D input data
        dims(1) = ntimesteps
        dims(2) = nx
        dims(3) = ny
    end function get_dims

    !>------------------------------------------------------------
    !!  Read in data from a given filename and variable into a pre-allocated output array
    !!
    !!  Handles 2d, 3d, and 4d file data
    !!    2D data are converted to a 1 x nx x ny
    !!    3D data have the third dimension (time) moved to the front
    !!    4D data have the last dimension (time) moved to the front,
    !!      and the 3rd dim subset to the first element for now
    !!
    !!------------------------------------------------------------
    subroutine load_data(varname, filenames, output, timestep, time_indices, time_weights, agg_method)
        implicit none
        character(len=MAXVARLENGTH),  intent(in)                    :: varname
        character(len=MAXFILELENGTH), intent(in),   dimension(:)    :: filenames
        real,                         intent(inout),dimension(:,:,:):: output
        integer,                        intent(in)                  :: timestep
        integer,                        intent(in), dimension(:)    :: time_indices
        real,                           intent(in), dimension(:)    :: time_weights
        integer,                        intent(in)                  :: agg_method

        ! array index counters
        integer :: file_idx, curstep, i
        integer, dimension(io_maxDims) :: dims

        ! array temporaries to store input data in
        real, allocatable, dimension(:,:) :: data_2d
        real, allocatable, dimension(:,:,:) :: data_3d
        real, allocatable, dimension(:,:,:,:) :: data_4d

        ! current time step in resultant data
        curstep = 1
        ! loop through input files
        do file_idx = 1,size(filenames,1)

            ! find the dimensions of the current file
            call io_getdims(filenames(file_idx), varname, dims)

            !---------------------------
            ! read in 3D + time data, subset 3rd dim
            !---------------------------
            if (dims(1)==4) then
                call io_read(filenames(file_idx), varname, data_4d)
                ! assign all timesteps to output data array
                if (timestep==-1) then
                    if (time_indices(1)==-1) then
                        do i=1,dims(5)
                            output(curstep,:,:) = data_4d(:,:,1,i)
                            curstep = curstep + 1
                        end do
                    else
                        output(curstep,:,:) = data_4d(:,:,1,time_indices(1))
                        if (agg_method == kAGG_TYPE_AVG) then
                            output(curstep,:,:) = output(curstep,:,:) * time_weights(1)
                        endif
                        do i=2,size(time_indices)
                            if (agg_method == kAGG_TYPE_AVG) then
                                output(curstep,:,:) = output(curstep,:,:) + data_4d(:,:,1,time_indices(i)) * time_weights(i)
                            else if (agg_method == kAGG_TYPE_SUM) then
                                output(curstep,:,:) = output(curstep,:,:) + data_4d(:,:,1,time_indices(i))
                            else if (agg_method == kAGG_TYPE_MIN) then
                                ! output(curstep,:,:) = array_minimum(output(curstep,:,:), data_4d(:,:,1,time_indices(i)))
                                stop "Not implemented: minimum aggregation type"
                            else if (agg_method == kAGG_TYPE_MAX) then
                                ! output(curstep,:,:) = array_maximum(output(curstep,:,:), data_4d(:,:,1,time_indices(i)))
                                stop "Not implemented: maximum aggregation type"
                            else
                                stop "Not implemented: unknown aggregation type"
                            endif
                        enddo
                        if (agg_method == kAGG_TYPE_AVG) then
                            output(curstep,:,:) = output(curstep,:,:) / sum(time_weights)
                        endif
                        curstep = curstep + 1
                    endif
                else
                    output(curstep,:,:) = data_4d(:,:,1,timestep)
                    curstep = curstep + 1
                endif
                ! deallocate input array because it is allocated in the read routine... bad form?
                deallocate(data_4d)
            !---------------------------
            ! read in 2D + time data move time to first dim
            !---------------------------
            else if (dims(1)==3) then
                call io_read(filenames(file_idx), varname, data_3d)
                ! assign all timesteps to output data array
                if (timestep==-1) then
                    if (time_indices(1)==-1) then
                        do i=1,dims(4)
                            output(curstep,:,:) = data_3d(:,:,i)
                            curstep = curstep + 1
                        end do
                    else
                        output(curstep,:,:) = data_3d(:,:,time_indices(1))
                        if (agg_method == kAGG_TYPE_AVG) then
                            output(curstep,:,:) = output(curstep,:,:) * time_weights(1)
                        endif
                        do i=2,size(time_indices)
                            if (agg_method == kAGG_TYPE_AVG) then
                                output(curstep,:,:) = output(curstep,:,:) + (data_3d(:,:,time_indices(i)) * time_weights(i))
                            else if (agg_method == kAGG_TYPE_SUM) then
                                output(curstep,:,:) = output(curstep,:,:) + data_3d(:,:,time_indices(i))
                            else if (agg_method == kAGG_TYPE_MIN) then
                                ! output(curstep,:,:) = array_minimum(output(curstep,:,:), data_4d(:,:,1,time_indices(i)))
                                stop "Not implemented: minimum aggregation type"
                            else if (agg_method == kAGG_TYPE_MAX) then
                                ! output(curstep,:,:) = array_maximum(output(curstep,:,:), data_4d(:,:,1,time_indices(i)))
                                stop "Not implemented: maximum aggregation type"
                            else
                                stop "Not implemented: unknown aggregation type"
                            endif
                        enddo
                        if (agg_method == kAGG_TYPE_AVG) then
                            output(curstep,:,:) = output(curstep,:,:) / sum(time_weights)
                        endif
                        curstep = curstep + 1
                    endif
                else
                    output(curstep,:,:) = data_3d(:,:,timestep)
                    curstep = curstep + 1
                endif
                ! deallocate input array because it is allocated in the read routine... bad form?
                deallocate(data_3d)
            !---------------------------
            ! read in 2D data, add a new (1 element) time dim
            !---------------------------
            else if (dims(1)==2) then
                ! read in 2d (no time) data
                call io_read(filenames(file_idx), varname, data_2d)
                output(curstep,:,:) = data_2d(:,:)
                ! just one time step to add to the current step
                curstep = curstep + 1
                deallocate(data_2d)
            endif
        end do

        if ( (curstep-1) /= size(output,1) ) then
            write(*,*) "Error: data read in (", trim(str(curstep-1)), &
                        ") did not match calculated time steps:", trim(str(size(output,1)))
            stop "Data loading error: Time dimension did not match"
        end if

    end subroutine load_data

end submodule gefs_io_implementation
