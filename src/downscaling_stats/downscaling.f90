module downscaling_mod

    use data_structures
    use string,             only : str
    use regression_mod,     only : compute_regression, compute_logistic_regression
    use analog_mod,         only : find_analogs, compute_analog_mean, compute_analog_error, compute_analog_exceedance
    use model_constants
    use quantile_mapping,   only : develop_qm, apply_qm
    use time_util,          only : setup_time_indices
    implicit none
    
    integer*8, dimension(10) :: master_timers
    real, parameter :: LOG_FILL_VALUE = 1e-30

    
contains
    function downscale(training_atm, training_obs, predictors, options) result(output)
            implicit none
            type(atm),    intent(inout) :: training_atm, predictors
            type(obs),    intent(inout) :: training_obs
            type(config), intent(in)    :: options
            
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
            ! variables to store timing data for profiling code
            integer*8 :: timeone, timetwo, master_timeone, master_timetwo, master_time_post_init
            integer*8, dimension(10) :: timers
            
            integer :: total_number_of_gridcells, current_completed_gridcells
            
            timers = 0
            master_timers = 0
            current_completed_gridcells = 0
            
            call System_Clock(timeone)
            ! should this be done outside of "downscale" ? 
            call setup_timing(training_atm, training_obs, predictors, options)
            ! simple local variables to make code more legible later
            p_start    = predictors%first_time
            p_end      = predictors%last_time
            t_tr_start = training_atm%training_start
            t_tr_stop  = training_atm%training_stop
            o_tr_start = training_obs%training_start
            o_tr_stop  = training_obs%training_stop
            
            ntimes     = size(predictors%variables(1)%data,1)
            ntrain     = size(training_atm%variables(1)%data,1)
            nobs       = size(training_obs%variables(1)%data,1)
            nx         = size(training_obs%variables(1)%data,2)
            ny         = size(training_obs%variables(1)%data,3)
            
            tr_size    = o_tr_stop - o_tr_start + 1
            noutput    = p_end - p_start + 1
            
            n_atm_variables = size(training_atm%variables)
            n_obs_variables = size(training_obs%variables)
            
            ! nx = 128
            ! ny = 128
            
            print*, "=========================================="
            print*, "Running for "
            print*, nx, "   by",ny, " Grid Cells"
            
            total_number_of_gridcells = nx * ny
            ! total_number_of_gridcells = nx * (ny-100)
            
            call System_Clock(timeone)
            call allocate_data(output, n_obs_variables, noutput, nx, ny, n_atm_variables, tr_size, options)
            print*, ""
            
            do v=1,n_obs_variables
                output%variables(v)%name = training_obs%variables(v)%name
                do j=1,ny
                    do i=1,nx
                        call transform_data(options%obs%input_Xforms(v), training_obs%variables(v)%data(:,i,j), 1, nobs)
                    enddo
                enddo
            enddo
            
            call System_Clock(timetwo)
            timers(7) = timetwo - timeone

            call System_Clock(master_timeone)
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!                                       !!
            !!  Begin parallelization                !!
            !!                                       !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$omp parallel default(shared)                  &
            !$omp private(train_data, pred_data, i, j, v, timeone, timetwo)   &
            !$omp firstprivate(nx,ny,n_atm_variables, n_obs_variables, noutput, ntrain, nobs, ntimes) &
            !$omp firstprivate(p_start, p_end, t_tr_start, t_tr_stop, o_tr_start, o_tr_stop, timers)
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Set up the data structures used in downscaling
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            allocate( train_data( ntrain, n_atm_variables+1 ) )
            allocate( pred_data(  ntimes, n_atm_variables+1 ) )
            train_data=0
            pred_data=0
            ! constant coefficient for regressions... might be better to keep this in the point downscaling section? 
            pred_data(:,1) = 1
            train_data(:,1) = 1
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Apply relevant transformations to the input atmospheric data
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            do v=1,n_atm_variables
                !$omp do
                do j=1,size(predictors%variables(v)%data,3)
                    do i=1,size(predictors%variables(v)%data,2)
                        call transform_data(options%prediction%input_Xforms(v), predictors%variables(v)%data(:,i,j), 1, ntimes)
                    enddo
                    ! does not need to be normalized if it is transformed to match training_atm anyway
                    if (options%prediction%transformations(v) /= kQUANTILE_MAPPING) then
                        call normalize(predictors%variables(v), j)
                    endif
                enddo
                !$omp end do

                !$omp do
                do j=1,size(training_atm%variables(v)%data,3)
                    do i=1,size(training_atm%variables(v)%data,2)
                        call transform_data(options%training%input_Xforms(v), training_atm%variables(v)%data(:,i,j), 1, ntrain)
                    enddo
                    call normalize(training_atm%variables(v), j)
                enddo
                !$omp end do
            enddo
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Normalization must be finished before running anything else
            !$omp barrier
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !$omp single
            write(*,*), "Downscaling..."
            call System_Clock(master_time_post_init)
            !$omp end single
            ! parallelization could be over x and y, (do n=1,ny*nx; j=n/nx; i=mod(n,nx)) and use schedule(dynamic)
            !$omp do schedule(static, 1)
            do j=1,ny
                do i=1,nx
                    
                    if (training_obs%mask(i,j)) then
                        
                        ! index 1 is the constant for the regression code... should be moved into regression code...?
                        if (options%debug) then
                            output%variables(1)%predictors(:,i,j,1) = pred_data( p_start   : p_end,    1)
                        endif
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!
                        !!  Set up the data structures to downscale the current point
                        !!
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        do v=1,n_atm_variables
                            call System_Clock(timeone)
                            call read_point( predictors%variables(v)%data,   pred_data(:,v+1),  i,j, predictors%geoLUT,   options%prediction%interpolation_method)
                            call read_point( training_atm%variables(v)%data, train_data(:,v+1), i,j, training_atm%geoLUT, options%training%interpolation_method)
                            call System_Clock(timetwo)
                            timers(4) = timers(4) + (timetwo-timeone)
                            
                            ! perform (e.g.) quantile mapping
                            call System_Clock(timeone)
                            call transform_data(options%prediction%transformations(v),                                         &
                                                pred_data(:,v+1),  predictors%transform_start,   predictors%transform_stop,    &
                                                train_data(:,v+1), training_atm%transform_start, training_atm%transform_stop)
                            
                            if (options%debug) then
                                ! save these data for output debugging / algorithm development while we are at it. 
                                output%variables(1)%predictors(:,i,j,v+1) = pred_data( p_start   : p_end,    v+1)
                                output%variables(1)%training  (:,i,j,v+1) = train_data(t_tr_start: t_tr_stop,v+1)
                            endif
                            call System_Clock(timetwo)
                            timers(6) = timers(6) + (timetwo-timeone)
                        enddo

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!
                        !!  Downscale all variables for this point
                        !!
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        do v=1,n_obs_variables
                            
                            if (options%debug) then
                                !$omp critical (print_lock)
                                write(*,*) ""
                                write(*,*) "-----------------------------"
                                write(*,*) "   Downscaling point: ",i,j
                                write(*,*) "   For variable     : ", trim(training_obs%variables(v)%name)
                                !$omp end critical (print_lock)
                                output%variables(v)%obs(:,i,j) = training_obs%variables(v)%data(o_tr_start:o_tr_stop, i, j)
                            endif
                            
                            ! As tempting as it may be, associate statements are not threadsafe!!!
                            output%variables(v)%data(:,i,j) = downscale_point(                                         &
                                                    pred_data                       (   p_start : p_end,     :),       &
                                                    train_data                      (t_tr_start : t_tr_stop, :),       &
                                                    training_obs%variables(v)%data  (o_tr_start : o_tr_stop, i, j),    &
                                                    output%variables(v)%errors      (           :,           i, j),    &
                                                    output%variables(v)%coefficients(           :,  :,       i, j),    &
                                                    output%variables(v)%logistic    (           :,           i, j),    &
                                                    output%variables(v)%logistic_threshold,                            &
                                                    options, timers )

                            ! if the input data were transformed with e.g. a cube root or log transform, then reverse that transformation for the output
                            call transform_data(options%obs%input_Xforms(v), output%variables(v)%data(:,i,j), 1, noutput, reverse=.True.)


                        enddo
                        
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!
                    !!  Else This is a masked point in the observations
                    !!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    else
                        if (options%debug) then
                            !$omp critical (print_lock)
                            write(*,*) ""
                            write(*,*) "-----------------------------"
                            write(*,*) "   Masked point: ",i,j
                            write(*,*) "-----------------------------"
                            !$omp end critical (print_lock)
                        endif
                        
                        ! store a fill value in the output
                        do v=1,n_obs_variables
                            output%variables(v)%data(:,i,j) = kFILL_VALUE
                        enddo
                        
                    endif
                    
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !!
                    !!  Print a status updates
                    !!
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    !$omp critical (print_lock)
                    current_completed_gridcells = current_completed_gridcells + 1
                    write(*,"(A,f5.1,A$)") char(13), current_completed_gridcells/real(total_number_of_gridcells)*100.0, " %"
                    !$omp end critical (print_lock)
                enddo
            enddo
            !$omp end do
            !$omp critical
            master_timers = timers + master_timers
            !$omp end critical
            !$omp end parallel
            call System_Clock(master_timetwo)
            master_timers(8) = master_timers(8) + (master_time_post_init - master_timeone) * 16
            master_timers(5) = master_timers(5) + (master_timetwo-master_timeone) * 16
            print*, ""
            print*, "---------------------------------------------------"
            print*, "           Time profiling information "
            print*, "---------------------------------------------------"
            print*, "Total Time : ",  nint(100.0),                     "%    "!,  nint(0.01*master_timers(5))
            print*, "---------------------------------------------------"
            print*, "Allocation : ",  nint((100.d0 * master_timers(7))/master_timers(5)),  "%    "!, nint(0.01*master_timers(7))
            print*, "Data Init  : ",  nint((100.d0 * master_timers(8))/master_timers(5)),  "%    "!, nint(0.01*master_timers(9))
            print*, "GeoInterp  : ",  nint((100.d0 * master_timers(4))/master_timers(5)),  "%    "!, nint(0.01*master_timers(4))
            print*, "Transform  : ",  nint((100.d0 * master_timers(6))/master_timers(5)),  "%    "!, nint(0.01*master_timers(6))
            print*, "Analog     : ",  nint((100.d0 * master_timers(1))/master_timers(5)),  "%    "!, nint(0.01*master_timers(1))
            print*, "Regression : ",  nint((100.d0 * master_timers(2))/master_timers(5)),  "%    "!, nint(0.01*master_timers(2))
            print*, "Log.Regres : ",  nint((100.d0 * master_timers(3))/master_timers(5)),  "%    "!, nint(0.01*master_timers(3))
            print*, "Log.Analog : ",  nint((100.d0 * master_timers(9))/master_timers(5)),  "%    "!, nint(0.01*master_timers(9))
            print*, "---------------------------------------------------"
            print*, "Parallelization overhead (?)"
            print*, "Residual   : ",  100-nint((100.0 * (master_timers(1)+master_timers(2)+master_timers(3)+master_timers(4)+master_timers(6)+master_timers(8)+master_timers(9))) / master_timers(5)),"%    "
            print*, "---------------------------------------------------"
            

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
                            t_xf_start, t_xf_stop,      &
                            reverse)
        implicit none
        integer,            intent(in)    :: transform_type
        real, dimension(:), intent(inout) :: pred_data
        integer,            intent(in)    :: p_xf_start, p_xf_stop
        real, dimension(:), intent(in), optional :: train_data
        integer,            intent(in), optional :: t_xf_start, t_xf_stop
        logical,            intent(in), optional :: reverse
        
        
        ! local variables needed by the quantile mapping transform
        type(qm_correction_type) :: qm
        real, dimension(:), allocatable :: temporary
        logical :: reverse_internal
        
        reverse_internal = .False.
        if (present(reverse)) reverse_internal = reverse
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
            
            if (reverse_internal) then
                pred_data = exp(pred_data)
                where(pred_data <= LOG_FILL_VALUE) pred_data = 0
            else
                where(pred_data <= 0) pred_data = LOG_FILL_VALUE
                pred_data = log(pred_data)
            endif
            
        case (kCUBE_ROOT)
            if (reverse_internal) then
                pred_data = pred_data ** 3
            else
                ! this may get optimized as a cube root by most compilers.
                pred_data = pred_data ** (1/3.0)
            endif
            
        case (kFIFTH_ROOT)
            if (reverse_internal) then
                pred_data = pred_data ** 5
            else
                pred_data = pred_data ** (1/5.0)
            endif
        end select
    end subroutine transform_data
    
    function downscale_point(predictor, atm, obs, errors, output_coeff, logistic, logistic_threshold, options, timers) result(output)
        implicit none
        real,    dimension(:,:), intent(inout):: predictor, atm ! (ntimes, nvars)
        real,    dimension(:),   intent(in)   :: obs
        real,    dimension(:),   intent(inout):: errors
        real,    dimension(:,:), intent(inout):: output_coeff
        real,    dimension(:),   intent(inout):: logistic
        real,                    intent(in)   :: logistic_threshold
        type(config),            intent(in)   :: options
        integer(8),dimension(:), intent(inout):: timers
        
        real,    dimension(:),   allocatable :: output
        real(8), dimension(:),   allocatable :: coefficients
        
        integer(8)  :: timeone, timetwo
        integer     :: i, n, nvars
        !integer     :: a, n_analogs, selected_analog, real_analogs
        !real        :: rand, analog_threshold
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!
        !!  Generic Initialization code
        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nvars = size(atm,2)
        n = size(predictor,1)
        
        allocate(coefficients(nvars))
        allocate(output(n))

        ! This just prevents any single points that were WAY out (most likely due to the QM?)
        where(predictor < -10) predictor = -10
        where(predictor >  10) predictor =  10

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!
        !!  Initialization code for pure regression
        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (options%pure_regression) then
            call System_Clock(timeone)
            output(1)  = compute_regression(predictor(1,:), atm, obs, coefficients, errors(1))
            errors(2:) = errors(1)
            if (options%debug) then
                do i=1,nvars
                    output_coeff(i,:) = coefficients(i)
                enddo
            endif
            call System_Clock(timetwo)
            timers(2) = timers(2) + (timetwo-timeone)

            call System_Clock(timeone)
            if (logistic_threshold/=kFILL_VALUE) then
                logistic(1) = compute_logistic_regression(predictor(1,:), atm, obs, coefficients, logistic_threshold)
                
                if (options%debug) then
                    do i = 1,nvars
                        output_coeff(i+nvars,:) = coefficients(i)
                    enddo
                endif
            endif
            call System_Clock(timetwo)
            timers(3) = timers(3) + (timetwo-timeone)

        endif
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!
        !!  Loop through time applying downscaling technique
        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i = 1, n
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Pure Analog downscaling
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (options%pure_analog) then
                call downscale_pure_analog(predictor(i,:), atm, obs, output_coeff(:,i),     &
                                           output(i), errors(i), logistic(i),               &
                                           options, timers)
                
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Analog Regression Downscaling
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (options%analog_regression) then
                call downscale_analog_regression(predictor(i,:), atm, obs, output_coeff(:,i),     &
                                                 output(i), errors(i), logistic(i),               &
                                                 options, timers)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Pure regression downscaling
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (options%pure_regression) then
                ! to test matmul(predictor, output_coeff(:,1)) should provide all time values efficiently?
                call apply_pure_regression(output(i), predictor(i,:), output_coeff(:,1), logistic(i), logistic_threshold, timers)
                
            endif
        end do

    end function downscale_point

    subroutine downscale_analog_regression(x, atm, obs, output_coeff, output, error, logistic, options, timers)
        implicit none
        real,           intent(in),     dimension(:,:)  :: atm
        real,           intent(in),     dimension(:)    :: x, obs
        real,           intent(inout),  dimension(:)    :: output_coeff
        real,           intent(inout)                   :: output, logistic, error
        type(config),   intent(in)                      :: options
        integer*8,      intent(inout),  dimension(:)    :: timers
        
        real,    dimension(:),   allocatable :: obs_analogs
        real,    dimension(:,:), allocatable :: regression_data
        integer, dimension(:),   allocatable :: analogs
        real(8), dimension(:),   allocatable :: coefficients
        integer*8   :: timeone, timetwo
        real        :: analog_threshold, logistic_threshold
        integer     :: real_analogs, selected_analog, nvars, n_analogs
        integer     :: j, a
        
        ! variables for handling packing e.g. precip to only use positive values in computing the amount
        real,    dimension(:,:),allocatable :: threshold_atm 
        real,    dimension(:),  allocatable :: threshold_obs
        logical, dimension(:),  allocatable :: threshold_packing
        integer :: n_packed
        
        n_analogs           = options%n_analogs
        analog_threshold    = options%analog_threshold
        logistic_threshold  = options%logistic_threshold
        nvars               = size(x)
                
        call System_Clock(timeone)
        call find_analogs(analogs, x, atm, n_analogs, analog_threshold)
        real_analogs = size(analogs)

        allocate(coefficients(nvars))
        allocate(regression_data(real_analogs, nvars))
        allocate(obs_analogs(real_analogs))
        if (logistic_threshold/=kFILL_VALUE) then
            allocate(threshold_packing(real_analogs))
        endif
        
        do a=1,real_analogs
            obs_analogs(a) = obs(analogs(a))
            regression_data(a,:) = atm(analogs(a),:)
        enddo
        call System_Clock(timetwo)
        timers(1) = timers(1) + (timetwo-timeone)

        call System_Clock(timeone)
        if (logistic_threshold==kFILL_VALUE) then
            output = compute_regression(x, regression_data, obs_analogs, coefficients, error)
        else
            
            threshold_packing = obs_analogs > logistic_threshold
            n_packed = count(threshold_packing)
            
            if (n_packed > nvars) then
                allocate(threshold_atm(n_packed, nvars))
                allocate(threshold_obs(n_packed))
                do j=1,nvars
                    threshold_atm(:,j) = pack(regression_data(:,j), threshold_packing)
                enddo
                threshold_obs = pack(obs_analogs, threshold_packing)
                
                output = compute_regression(x, threshold_atm, threshold_obs, coefficients, error)
                
            elseif (n_packed > 0) then
                output = sum(pack(obs_analogs, threshold_packing)) / n_packed
            else
                output = compute_regression(x, regression_data, obs_analogs, coefficients, error)
            endif
            
        endif
        if (options%debug) then
            output_coeff(1:nvars) = coefficients
        endif
        call System_Clock(timetwo)
        timers(2) = timers(2) + (timetwo-timeone)

        call System_Clock(timeone)
        if (logistic_threshold/=kFILL_VALUE) then
            if (options%logistic_from_analog_exceedance) then
                
                ! if the user specified using fewer analogs specifically for the logistic component, find them here
                if ((options%n_log_analogs /= real_analogs).and.(options%n_log_analogs > 0)) then
                    
                    ! find the logistic analogs
                    call find_analogs(analogs, x, atm, options%n_log_analogs, -1.0)
                    call System_Clock(timetwo)
                    timers(9) = timers(9) + (timetwo-timeone)
                    
                    ! compute the logistic as the probability of threshold exceedance in the analog population
                    call System_Clock(timeone)
                    logistic = compute_analog_exceedance(obs, analogs(1:options%n_log_analogs), logistic_threshold)
                else
                    ! compute the logistic as the probability of threshold exceedance in the analog population
                    logistic = compute_analog_exceedance(obs, analogs, logistic_threshold)
                endif
                
            else
                logistic = compute_logistic_regression(x, regression_data, obs_analogs, coefficients, logistic_threshold)
                if (options%debug) then
                    do j = 1,nvars
                        output_coeff(j+nvars) = coefficients(j)
                    enddo
                endif
            endif
        endif
        
        call System_Clock(timetwo)
        timers(3) = timers(3) + (timetwo-timeone)
    end subroutine downscale_analog_regression

    subroutine downscale_pure_analog(x, atm, obs, output_coeff, output, error, logistic, options, timers)
        implicit none
        real,           intent(in),     dimension(:,:)  :: atm
        real,           intent(in),     dimension(:)    :: x, obs
        real,           intent(inout),  dimension(:)    :: output_coeff
        real,           intent(inout)                   :: output, logistic, error
        type(config),   intent(in)                      :: options
        integer*8,      intent(inout),  dimension(:)    :: timers
        
        real        :: analog_threshold, logistic_threshold
        integer*8   :: timeone, timetwo
        integer     :: real_analogs, selected_analog, nvars, n_analogs
        integer, dimension(:), allocatable :: analogs
        real        :: rand
        integer     :: j
        
        n_analogs           = options%n_analogs
        analog_threshold    = options%analog_threshold
        logistic_threshold  = options%logistic_threshold
        nvars               = size(x)
        
        if (n_analogs > 0) then
            allocate(analogs(n_analogs))
        endif
        
        call System_Clock(timeone)
        ! find the best n_analog matching analog days in atm to match x
        call find_analogs(analogs, x, atm, n_analogs, analog_threshold)
        real_analogs = size(analogs)
        call System_Clock(timetwo)
        timers(1) = timers(1) + (timetwo-timeone)
        
        call System_Clock(timeone)
        if (options%sample_analog) then
            call random_number(rand)
            selected_analog = floor(rand * real_analogs)+1
            
            output = obs( analogs(selected_analog) )
            if (options%debug) then
                output_coeff(1:nvars) = atm(analogs(selected_analog), :)
            endif
        else
            output = compute_analog_mean(obs, analogs)
            
            if (options%debug) then
                do j=1,nvars
                    output_coeff(j) = compute_analog_mean(atm(:, j), analogs)
                enddo
            endif
        endif
        
        error = compute_analog_error(obs, analogs, output)
        call System_Clock(timetwo)
        timers(2) = timers(2) + (timetwo-timeone)
        
        call System_Clock(timeone)
        if (logistic_threshold/=kFILL_VALUE) then
            logistic = compute_analog_exceedance(obs, analogs, logistic_threshold)
        endif
        call System_Clock(timetwo)
        timers(3) = timers(3) + (timetwo-timeone)
        
    end subroutine downscale_pure_analog

    subroutine apply_pure_regression(output, x, B, logistic, threshold, timers)
        implicit none
        real,       intent(in),     dimension(:)    :: x, B
        real,       intent(inout)                   :: output, logistic
        real,       intent(in)                      :: threshold
        integer*8,  intent(inout),  dimension(:)    :: timers
        
        integer*8   :: timeone, timetwo
        integer     :: nvars
        
        nvars = size(x)
        call System_Clock(timeone)
        output = dot_product(x, B(1:nvars))
        call System_Clock(timetwo)
        timers(2) = timers(2) + (timetwo-timeone)

        call System_Clock(timeone)
        if (threshold/=kFILL_VALUE) then
            logistic = 1.0 / (1.0 + exp(-dot_product(x, B(nvars+1:nvars*2))))
        endif
        call System_Clock(timetwo)
        timers(3) = timers(3) + (timetwo-timeone)

    end subroutine apply_pure_regression
    
    
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
                !$omp critical (print_lock)
                write(*,*) "ERROR Normalizing:", trim(var%name)
                write(*,*) "  For point: ", i, j
                !$omp end critical (print_lock)
            endif
        enddo
        
    end subroutine normalize
    
    
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
    
    subroutine allocate_data(output, n_obs_variables, noutput, nx, ny, n_atm_variables, tr_size, options)
        implicit none
        type(results),  intent(inout)   :: output
        integer,        intent(in)      :: n_obs_variables, noutput, nx, ny, n_atm_variables, tr_size
        type(config),   intent(in)      :: options
        
        integer :: v
        
        write(*,*), "Allocating Memory"
        
        allocate(output%variables(n_obs_variables))
        
        do v = 1,n_obs_variables
            allocate(output%variables(v)%data        (noutput, nx, ny))
            allocate(output%variables(v)%errors      (noutput, nx, ny))
            
            if (options%logistic_threshold/=kFILL_VALUE) then
                allocate(output%variables(v)%logistic    (noutput, nx, ny))
                allocate(output%variables(v)%coefficients((n_atm_variables+1)*2, noutput, nx, ny))
            else
                ! unfortunately logistic has to be allocated even if it is not used because it is indexed 
                ! however, it does not need a complete time sequence because time is never indexed
                allocate(output%variables(v)%logistic    (2, nx, ny))
                allocate(output%variables(v)%coefficients(n_atm_variables+1, noutput, nx, ny))
            endif
            
            output%variables(v)%data            = 0
            output%variables(v)%logistic        = 0
            
            output%variables(v)%logistic_threshold = options%logistic_threshold
            
            if (options%obs%input_Xforms(v) == kLOG_TRANSFORM) then
                if (options%logistic_threshold<=0) then
                    output%variables(v)%logistic_threshold = log(LOG_FILL_VALUE)
                endif
            endif
            output%variables(v)%coefficients    = 0
                        
            if (options%debug) then
                allocate(output%variables(v)%obs         (tr_size, nx, ny))
                allocate(output%variables(v)%training    (tr_size, nx, ny, n_atm_variables+1))
                allocate(output%variables(v)%predictors  (noutput, nx, ny, n_atm_variables+1))
                output%variables(v)%predictors  = 0
                output%variables(v)%training    = 0
                output%variables(v)%obs         = 0
            endif
            
        end do
    end subroutine allocate_data

end module downscaling_mod
