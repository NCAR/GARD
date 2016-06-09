module downscaling_mod

    use data_structures
    use model_constants
    use quantile_mapping,   only : develop_qm, apply_qm
    use time_util,          only : setup_time_indices
    implicit none
    
contains
    function downscale(training_atm, training_obs, predictors, options) result(output)
            implicit none
            type(atm),    intent(inout) :: training_atm, predictors
            type(obs),    intent(inout) :: training_obs
            type(config), intent(in) :: options
            
            type(results) :: output
            type(qm_correction_type) :: qm
            real, dimension(:,:), allocatable :: train_data, pred_data
            integer :: nx, ny, ntimes, ntrain, nobs, noutput
            integer :: n_obs_variables, n_atm_variables
            integer :: i, j, l, x, y, v
            real :: w
            integer :: p_start, p_end, t_tr_start, t_tr_stop, o_tr_start, o_tr_stop
            
            ! should this be done outside of "downscale" ? 
            call setup_timing(training_atm, training_obs, predictors, options)
            ! simple local variables to make code more legible later
            p_start    = predictors%first_time
            p_end      = predictors%last_time
            t_tr_start = training_atm%training_start
            t_tr_stop  = training_atm%training_stop
            o_tr_start = training_obs%training_start
            o_tr_stop  = training_obs%training_stop


            
            ntimes = size(predictors%variables(1)%data,1)
            ntrain = size(training_atm%variables(1)%data,1)
            nobs   = size(training_obs%variables(1)%data,1)
            nx     = size(training_obs%variables(1)%data,2)
            ny     = size(training_obs%variables(1)%data,3)
            
            noutput = predictors%last_time - predictors%first_time + 1
            
            n_atm_variables = size(training_atm%variables)
            n_obs_variables = size(training_obs%variables)

            allocate( train_data( ntrain, n_atm_variables+1 ) )
            allocate( pred_data(  ntimes, n_atm_variables+1 ) )
            allocate(output%variables(n_obs_variables))
            do i = 1,n_obs_variables
                allocate(output%variables(i)%data(noutput, nx, ny))
                allocate(output%variables(i)%errors(noutput, nx, ny))
            end do
            
            ! constant coefficient for regressions... might be better to keep this in the point downscaling section? 
            pred_data(:,1) = 1
            train_data(:,1) = 1
            do j=1,ny
                do i=1,nx
                    
                    ! if (training_obs%mask(i,j)) then
                        
                        do v=1,n_atm_variables
                            
                            pred_data(:,v+1)  = read_point( predictors%variables(v)%data,   i,j, predictors%geoLUT)
                            train_data(:,v+1) = read_point( training_atm%variables(v)%data, i,j, training_atm%geoLUT)
                            
                            ! perform (e.g.) quantile mapping
                            call transform_data(options%prediction%transformations(v),                                         &
                                                pred_data(:,v+1),  predictors%transform_start,   predictors%transform_stop,    &
                                                train_data(:,v+1), training_atm%transform_start, training_atm%transform_stop)

                        enddo
                        
                        do v=1,n_obs_variables
                            associate(                                       &
                                var        => output%variables(v)%data,      &
                                errors     => output%variables(v)%errors,    &
                                observed   => training_obs%variables(v)%data(o_tr_start:o_tr_stop, i, j) &
                                )

                                var(:,i,j) = downscale_point( pred_data( p_start   : p_end,    :),  &
                                                              train_data(t_tr_start: t_tr_stop,:),  & 
                                                              observed, errors(:,i,j), options )
                            end associate
                            
                            if (options%debug) then
                                print*, ""
                                print*, "-----------------------------"
                                print*, "   Downscaling point: ",i,j
                                print*, "   For variable     : ", v
                            endif
                            
                            call transform_data(kQUANTILE_MAPPING,                                  &
                                                output%variables(v)%data(:,i,j),       1, noutput,  &
                                                training_obs%variables(v)%data(:,i,j), 1, nobs)
                        enddo

                        
                    ! endif
                enddo
            enddo

            !     call develop_qm(var(p_start:p_end,6,1), &
            !                     observed(o_tr_start:o_tr_stop, 135,100), qm, n_segments = 300)
            ! 
            !     call apply_qm(var(p_start:p_end,6,1), var(p_start:p_end,7,1), qm)
            ! 
            !     
            !     end associate
            ! end do
            ! 
            
    end function downscale
    
    function read_point(input_data, i, j, geolut) result(output)
        implicit none
        real, dimension(:,:,:), intent(in) :: input_data
        integer,                intent(in) :: i, j
        type(geo_look_up_table),intent(in) :: geolut
        
        real, dimension(:), allocatable :: output
        integer :: k,x,y
        real    :: w
        
        ! create the output to be the same length as first dimension of the input 
        allocate( output( size(input_data,1) ) )
        
        ! interpolation is performed using the bilinear weights computing in the geolut
        do k = 1,4
            ! %x contains the x coordinate each elements
            x = geolut%x(k,i,j)
            ! %y contains the y coordinate each elements
            y = geolut%y(k,i,j)
            ! %w contains the weight each elements
            w = geolut%w(k,i,j)
            output = output + input_data(:, x, y) * w
        end do

        
    end function read_point
    
    subroutine transform_data(transform_type, pred_data,  &
                            p_xf_start, p_xf_stop,      &
                            train_data,                 &
                            t_xf_start, t_xf_stop)
        implicit none
        integer,            intent(in)    :: transform_type
        real, dimension(:), intent(inout) :: pred_data
        integer,            intent(in)    :: p_xf_start, p_xf_stop
        real, dimension(:), intent(in), optional :: train_data
        integer,            intent(in), optional :: t_xf_start, t_xf_stop
        
        ! local variables needed by the quantile mapping transform
        type(qm_correction_type) :: qm
        real, dimension(:), allocatable :: temporary
        ! transform the relevant data (e.g. quantile map or other)

        select case (transform_type)
        case (kQUANTILE_MAPPING)
            
            allocate( temporary( size(pred_data)) )
            
            call develop_qm(pred_data( p_xf_start:p_xf_stop), &
                            train_data(t_xf_start:t_xf_stop), &
                            qm, n_segments = N_ATM_QM_SEGMENTS)

            call apply_qm(pred_data, temporary, qm)
            
            pred_data = temporary
            
        case (kLOG_TRANSFORM)
            pred_data = log(pred_data)
        case (kCUBE_ROOT)
            ! this may get optimized as a cube root by most compilers.
            pred_data = pred_data ** 1/3.0
        end select

        
    end subroutine transform_data
    
    function downscale_point(predictor, atm, obs, errors, options) result(output)
        implicit none
        real,    dimension(:,:), intent(in)  :: predictor, atm
        real,    dimension(:),   intent(in)  :: obs, errors
        type(config),            intent(in)  :: options
        
        real,    dimension(:), allocatable :: output
        integer, dimension(:), allocatable :: analogs
        
        integer :: i,n, n_analogs, selected_analog
        real    :: rand
        
        n_analogs = options%n_analogs
        n = size(predictor)
        allocate(output(n))
        allocate(analogs(n_analogs))
        
        do i = 1, n
            analogs = find_analogs(predictor(i,:), atm, n_analogs)
            call random_number(rand)
            selected_analog = floor(rand * n_analogs)+1
            print*, analogs(selected_analog)
            output(i) = obs( analogs(selected_analog) )
        end do
        
    end function downscale_point
    
    function find_analogs(match, input, n) result(analogs)
        implicit none
        real,    intent(in), dimension(:)   :: match
        real,    intent(in), dimension(:,:) :: input
        integer, intent(in)     :: n
        integer, dimension(1:n)   :: analogs
        
        real,    dimension(:), allocatable :: distances
        logical, dimension(:), allocatable :: mask
        integer, dimension(1) :: min_location
        integer :: i, n_inputs, nvars
        
        ! if (match==0) then
        !     analogs = pick_n_random_zeros(input, n)
        ! else
            
            n_inputs = size(input,1)
            nvars    = size(input,2)
            allocate(mask(n_inputs))
            allocate(distances(n_inputs))
            distances = 0
            mask = .True.
            
            do i=1,nvars
                distances = distances + abs(input(:,i) - match(i))
            enddo
            
            do i = 1, n
                min_location = minloc(distances,mask=mask)
                mask(min_location) = .false.
                
                analogs(i) = min_location(1)
            end do
        ! endif
        
    end function find_analogs
    
    function pick_n_random_zeros(input, n) result(analogs)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in)     :: n
        integer, dimension(1:n)   :: analogs
        integer, dimension(:), allocatable :: zero_locations
        
        integer :: i, n_inputs, n_zeros, selected_analog
        real :: rand

        n_inputs = size(input)
        allocate(zero_locations(n_inputs))
        n_zeros = 0
        do i=1,n_inputs
            if (input(i)==0) then
                n_zeros = n_zeros + 1
                zero_locations(n_zeros) = i
            end if
        end do
        
        do i=1,n
            call random_number(rand)
            selected_analog = floor(rand * n_zeros)+1
            analogs(i) = zero_locations(selected_analog)
        end do

    end function pick_n_random_zeros
    
    subroutine setup_timing(training_atm, training_obs, predictors, options)
        implicit none
        type(atm),    intent(inout) :: training_atm, predictors
        type(obs),    intent(inout) :: training_obs
        type(config), intent(in) :: options
        
        print*, "-------------------"
        print*, "Training ATM"
        call setup_time_indices(training_atm, options)
        print*, "-------------------"
        print*, "Training Obs"
        call setup_time_indices(training_obs, options)
        print*, "-------------------"
        print*, "Predictors"
        call setup_time_indices(predictors, options)

    end subroutine setup_timing

end module downscaling_mod

! results type information for reference:
! type(obs_variable_type), allocatable, dimension(:) :: variables
    ! character(len=MAXVARLENGTH)         :: name      ! name of the variable
    ! real, allocatable, dimension(:,:,:) :: data      ! raw data
    ! integer                             :: data_type ! Type of data.  e.g. precip, temperature, or other. 
    ! real, allocatable, dimension(:,:) :: mean, stddev ! per gridpoint mean and standard deviation (for normalization?)
    ! integer :: transformation                         ! type of transformation applied to data (e.g. sqrt, log, ???)
! type(Time_type), allocatable, dimension(:) :: times
! integer :: n_variables, n_times
! character (len=MAXSTRINGLENGTH) :: name
