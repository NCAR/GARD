module downscaling_mod

    use data_structures
    use regression_mod,     only : compute_regression
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
            ! prediction period index variables
            integer :: p_start, p_end
            ! training index variables 
            integer :: t_tr_start, t_tr_stop, o_tr_start, o_tr_stop, tr_size
            
            ! should this be done outside of "downscale" ? 
            call setup_timing(training_atm, training_obs, predictors, options)
            ! simple local variables to make code more legible later
            p_start    = predictors%first_time
            p_end      = predictors%last_time
            t_tr_start = training_atm%training_start
            t_tr_stop  = training_atm%training_stop
            o_tr_start = training_obs%training_start
            o_tr_stop  = training_obs%training_stop
            tr_size    = o_tr_stop - o_tr_start + 1


            
            ntimes = size(predictors%variables(1)%data,1)
            ntrain = size(training_atm%variables(1)%data,1)
            nobs   = size(training_obs%variables(1)%data,1)
            nx     = size(training_obs%variables(1)%data,2)
            ny     = size(training_obs%variables(1)%data,3)
            
            noutput = p_end - p_start + 1
            
            n_atm_variables = size(training_atm%variables)
            n_obs_variables = size(training_obs%variables)
            
            write(*,*) "N atmospheric variables: ", n_atm_variables
            write(*,*) "N observed variables: ", n_obs_variables
            
            nx = 32
            ny = 132
            ! ny = 180

            allocate(output%variables(n_obs_variables))
            do i = 1,n_obs_variables
                allocate(output%variables(i)%data       (noutput, nx, ny))
                allocate(output%variables(i)%errors     (noutput, nx, ny))
                allocate(output%variables(i)%obs        (tr_size, nx, ny))
                allocate(output%variables(i)%training   (tr_size, nx, ny, n_atm_variables+1))
                allocate(output%variables(i)%predictors (noutput, nx, ny, n_atm_variables+1))
                
                output%variables(i)%data        = 0
                output%variables(i)%predictors  = 0
                output%variables(i)%training    = 0
                output%variables(i)%obs         = 0
            end do
            
            !$omp parallel default(shared)                  &
            !$omp private(train_data, pred_data, i, j, v)   &
            !$omp firstprivate(nx,ny,n_atm_variables, n_obs_variables, noutput, ntrain, nobs, ntimes) &
            !$omp firstprivate(p_start, p_end, t_tr_start, t_tr_stop, o_tr_start, o_tr_stop)
            allocate( train_data( ntrain, n_atm_variables+1 ) )
            allocate( pred_data(  ntimes, n_atm_variables+1 ) )
            train_data=0
            pred_data=0
            ! constant coefficient for regressions... might be better to keep this in the point downscaling section? 
            pred_data(:,1) = 1
            train_data(:,1) = 1
            
            do v=1,n_atm_variables
                ! does not need to be normalized if it is transformed to match training_atm anyway
                if (options%prediction%transformations(v) /= kQUANTILE_MAPPING) then
                    !$omp do
                    do j=1, size(predictors%variables(v)%data,3)
                        call normalize(predictors%variables(v), j)
                    enddo
                    !$omp end do
                endif
                !$omp do
                do j=1, size(training_atm%variables(v)%data,3)
                    call normalize(training_atm%variables(v), j)
                enddo
                !$omp end do
            enddo
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Normalization must be finished before running anything else
            !$omp barrier
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            !$omp do schedule(static, 1)
            do j=100,ny
                do i=1,nx
                    
                    if (training_obs%mask(i,j)) then
                        
                        do v=1,n_atm_variables
                            
                            call read_point( predictors%variables(v)%data,   pred_data(:,v+1),  i,j, predictors%geoLUT,   options%prediction%interpolation_method)
                            call read_point( training_atm%variables(v)%data, train_data(:,v+1), i,j, training_atm%geoLUT, options%training%interpolation_method)
                            
                            ! perform (e.g.) quantile mapping
                            call transform_data(options%prediction%transformations(v),                                         &
                                                pred_data(:,v+1),  predictors%transform_start,   predictors%transform_stop,    &
                                                train_data(:,v+1), training_atm%transform_start, training_atm%transform_stop)
                            
                            ! save these data for output debugging / algorithm development while we are at it. 
                            output%variables(1)%predictors(:,i,j,v+1) = pred_data( p_start   : p_end,    v+1)
                            output%variables(1)%training  (:,i,j,v+1) = train_data(t_tr_start: t_tr_stop,v+1)
                            
                        enddo

                        do v=1,n_obs_variables
                            
                            if (options%debug) then
                                !$omp critical
                                write(*,*) ""
                                write(*,*) "-----------------------------"
                                write(*,*) "   Downscaling point: ",i,j
                                write(*,*) "   For variable     : ", v
                                !$omp end critical
                            endif
                            
                            ! removed associate statment because it was causing problems with parallelization with ifort...
                            output%variables(v)%data(:,i,j) = downscale_point(                                          &
                                                            pred_data    (   p_start : p_end,     :),        &
                                                            train_data   (t_tr_start : t_tr_stop, :),        & 
                                                            training_obs%variables(v)%data(o_tr_start : o_tr_stop, i, j),  &
                                                            output%variables(v)%errors(:,i,j), options )
                                                          
                            output%variables(v)%obs       (:,i,j)   = training_obs%variables(v)%data(o_tr_start:o_tr_stop, i, j)
                            
                            ! should be possible to run this here, but seems to run into a shared memory problem
                            ! call transform_data(kQUANTILE_MAPPING,                                  &
                            !                     output%variables(v)%data(:,i,j),       1, noutput,  &
                            !                     training_obs%variables(v)%data(:,i,j), 1, nobs)
                        enddo
                    else
                        do v=1,n_obs_variables
                            output%variables(v)%data(:,i,j) = kFILL_VALUE
                        enddo
                        if (options%debug) then
                            !$omp critical
                            write(*,*) ""
                            write(*,*) "-----------------------------"
                            write(*,*) "   Masked point: ",i,j
                            write(*,*) "-----------------------------"
                            !$omp end critical
                        endif
                    endif
                enddo
            enddo
            !$omp end do
            !$omp end parallel

    end function downscale
    
    subroutine read_point(input_data, output, i, j, geolut, method)
        implicit none
        real, dimension(:,:,:), intent(in) :: input_data
        integer,                intent(in) :: i, j
        type(geo_look_up_table),intent(in) :: geolut
        integer,                intent(in) :: method
        
        real, dimension(:),     intent(inout) :: output
        integer :: k,x,y
        real    :: w
        
        ! create the output to be the same length as first dimension of the input 
        ! allocate( output( size(input_data,1) ) )
        
        ! interpolation is performed using the bilinear weights computing in the geolut
        select case (method)
        case (kBILINEAR)
            output = 0
            do k = 1,4
                ! %x contains the x coordinate each elements
                x = geolut%x(k,i,j)
                ! %y contains the y coordinate each elements
                y = geolut%y(k,i,j)
                ! %w contains the weight each elements
                w = geolut%w(k,i,j)
                output = output + input_data(:, x, y) * w
            end do
        case (kNEAREST)
            ! %x contains the nearest x coordinate
            x = geolut%x(1,i,j)
            ! %y contains the nearest y coordinate
            y = geolut%y(1,i,j)
            output = input_data(:, x, y)
        end select

        
    end subroutine read_point
    
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
        
        real,    dimension(:),   allocatable :: output
        real,    dimension(:),   allocatable :: obs_analogs
        real,    dimension(:,:), allocatable :: regression_data
        integer, dimension(:),   allocatable :: analogs
        real(8), dimension(:),   allocatable :: coefficients
        
        integer :: i, j, n, nvars
        integer :: a, n_analogs, selected_analog
        real    :: rand
        
        n_analogs = options%n_analogs
        nvars = size(atm,2)
        n = size(predictor,1)
        
        allocate(coefficients(nvars))
        allocate(output(n))
        allocate(analogs(n_analogs))
        
        if (options%analog_regression) then
            allocate(regression_data(n_analogs, nvars))
            allocate(obs_analogs(n_analogs))
        endif
        
        if (options%pure_regression) then
            output(1) = compute_regression(predictor(1,:), atm, obs, coefficients)
        endif
        
        do i = 1, n
            if (.not.options%pure_regression) then
                analogs = find_analogs(predictor(i,:), atm, n_analogs)
            endif
            
            if (options%pure_analog) then
                
                call random_number(rand)
                selected_analog = floor(rand * n_analogs)+1

                output(i) = obs( analogs(selected_analog) )
                
            elseif (options%analog_regression) then
                do a=1,n_analogs
                    obs_analogs(a) = obs(analogs(a))
                    regression_data(a,:) = atm(analogs(a),:)
                enddo
                output(i) = compute_regression(predictor(i,:), regression_data, obs_analogs, coefficients)
                
            elseif (options%pure_regression) then
                output(i) = 0
                do j=1,nvars
                    output(i) = output(i) + predictor(i,j) * coefficients(j)
                enddo
            endif
            
        end do
        
    end function downscale_point
    
    ! normalize an atmospheric variable by subtracting the mean and dividing by the standard deviation
    ! Assumes mean and stddev have already been calculated as part of the data structure
    ! If stddev is 0, then no division occurs (mean is still removed)
    subroutine normalize(var, j)
        implicit none
        class(atm_variable_type), intent(inout) :: var
        integer, intent(in) :: j
        integer :: i, n
        
        n = size(var%mean,1)
        
        do i=1,n
            var%data(:,i,j) = var%data(:,i,j) - var%mean(i,j)
            
            if (var%stddev(i,j) /= 0) then
                var%data(:,i,j) = var%data(:,i,j) / var%stddev(i,j)
            else
                write(*,*) "ERROR Normalizing:", trim(var%name)
                write(*,*) "  For point: ", i, j
            endif
        enddo
        
    end subroutine normalize
    
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
                distances = distances + (input(:,i) - match(i))**2
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
        
        write(*,*) "-------------------"
        write(*,*) "Training ATM"
        call setup_time_indices(training_atm, options)
        write(*,*) "-------------------"
        write(*,*) "Training Obs"
        call setup_time_indices(training_obs, options)
        write(*,*) "-------------------"
        write(*,*) "Predictors"
        call setup_time_indices(predictors, options)

    end subroutine setup_timing

end module downscaling_mod
