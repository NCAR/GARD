!>------------------------------------------------
!! Handle all IO for GCM data
!! 
!! Loops through GCM variables reading data and computing basic statistics
!! Reads time and lat/lon variables from the GCM files as well. 
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module gcm_mod

    use data_structures
    use model_constants
    use basic_stats_mod,only: time_mean, time_stddev
    use string,         only: str
    use io_routines,    only: io_read, io_getdims, io_maxDims
    use time_io,        only: read_times
    use geo,            only: standardize_coordinates
    
    implicit none
    logical :: debug
    
contains
    
    !>------------------------------------------------
    !! Initialize the GCM module
    !!
    !!------------------------------------------------
    subroutine init_gcm_io(options)
        implicit none
        type(config), intent(in) :: options
        
        debug = options%debug
        
    end subroutine init_gcm_io
    
    !>------------------------------------------------
    !! Read in the GCM data
    !!
    !! Uses variable names and filenames defined in the options data structure
    !!
    !! Loops through all input variables, then reads lat, lon, and time variables
    !!
    !!------------------------------------------------
    function read_gcm(options) result(gcm_data)
        implicit none
        class(input_config), intent(in) :: options
        type(atm) :: gcm_data
        
        integer :: var_idx, ntimesteps
        integer :: nx,ny
        
        ! allocate space to store all of the variables to be read
        allocate( gcm_data%variables( options%n_variables ))
        gcm_data%name = options%name
        
        ! loop over variables reading them in
        do var_idx = 1, options%n_variables
            associate(var => gcm_data%variables(var_idx))
                
                var = read_gcm_variable( options%var_names(var_idx),        &
                                         options%file_names(:, var_idx))
                
                if (var_idx==1) then
                    ntimesteps = size(var%data, 1)
                else
                    if (ntimesteps /= size(var%data, 1)) then
                        write(*,*) "GCM 1st variable ntime:", trim(str( ntimesteps )),  &
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
        
        allocate(gcm_data%times(ntimesteps))
        call read_times(options, gcm_data%times)
        
        call io_read(options%file_names(1, 1), options%lat_name, gcm_data%lat)
        call io_read(options%file_names(1, 1), options%lon_name, gcm_data%lon)
        call standardize_coordinates(gcm_data)
    end function read_gcm
    
    subroutine compute_grid_stats(var)
        implicit none
        type(atm_variable_type), intent(inout) :: var
        
        integer :: nx, ny
        
        nx = size(var%data,2)
        ny = size(var%data,3)
        
        if (allocated(var%mean))    deallocate(var%mean)
        if (allocated(var%stddev))  deallocate(var%stddev)
        allocate(var%mean(nx,ny))
        allocate(var%stddev(nx,ny))
        
        where(var%data>1e10) var%data=0
        
        call time_mean( var%data, var%mean )
        call time_stddev( var%data, var%stddev, mean_in=var%mean )
        
    end subroutine compute_grid_stats
    
    function read_gcm_variable(varname, filenames) result(output)
        implicit none
        character(len=MAXVARLENGTH),    intent(in)              :: varname
        character(len=MAXFILELENGTH),   intent(in), dimension(:):: filenames
        
        type(atm_variable_type) :: output
        
        integer, dimension(io_maxDims) :: dims
        
        output%name = varname
        
        dims = get_dims(varname, filenames)
        
        ! note, we reverse the order of the dimensions here to speed up later computations which will occur per grid cell over time
        allocate(output%data(dims(1), dims(2), dims(3)))
        
        call load_data(varname, filenames, output%data)
        
    end function read_gcm_variable
    
    !! requires all filenames to have the same number of time steps... no good for monthly files...
    function get_dims(varname, filenames) result(dims)
        implicit none
        character(len=MAXVARLENGTH),    intent(in)               :: varname
        character(len=MAXFILELENGTH),   intent(in), dimension(:) :: filenames
        
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
                ! the last dimension is assumed to be time (+1 because ndims takes a slot)
                ntimesteps = ntimesteps + dims( dims(1)+1 )
            endif
            
            if (file_idx == 1) then
                nx = dims(2)
                ny = dims(3)
            else
                if ((nx/=dims(2)).or.(ny/=dims(3))) then
                    write(*,*) "GCM 1st nx:", trim(str(nx)), "current nx:", trim(str(dims(2)))
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
    subroutine load_data(varname, filenames, output)
        implicit none
        character(len=MAXVARLENGTH),  intent(in)                    :: varname
        character(len=MAXFILELENGTH), intent(in),   dimension(:)    :: filenames
        real,                         intent(inout),dimension(:,:,:):: output
        
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
                do i=1,dims(5)
                    output(curstep,:,:) = data_4d(:,:,1,i)
                    curstep = curstep + 1
                end do
                ! deallocate input array because it is allocated in the read routine... bad form? 
                deallocate(data_4d)
            !---------------------------
            ! read in 2D + time data move time to first dim
            !---------------------------
            else if (dims(1)==3) then
                call io_read(filenames(file_idx), varname, data_3d)
                ! assign all timesteps to output data array
                do i=1,dims(4)
                    output(curstep,:,:) = data_3d(:,:,i)
                    curstep = curstep + 1
                end do
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
    
end module gcm_mod
