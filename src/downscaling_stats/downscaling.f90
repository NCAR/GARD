module downscaling_mod

    use data_structures
    use string,             only : str
    use regression_mod,     only : compute_regression, compute_logistic_regression
    use analog_mod,         only : find_analogs, compute_analog_mean, compute_analog_error, compute_analog_exceedance
    use model_constants
    use quantile_mapping,   only : develop_qm, apply_qm
    use time_util,          only : setup_time_indices
    use basic_stats_mod,    only : stddev
    implicit none

    integer*8, dimension(10) :: master_timers
    real, parameter :: LOG_FILL_VALUE = 1e-30


contains
    subroutine downscale(training_atm, training_obs, predictors, output, options)
            implicit none
            type(atm),    intent(inout) :: training_atm, predictors
            type(obs),    intent(inout) :: training_obs
            type(results),intent(out)   :: output
            type(config), intent(in)    :: options

            type(qm_correction_type) :: qm
            real, dimension(:,:), allocatable :: train_data, pred_data
            integer :: nx, ny, ntimes, ntrain, nobs, noutput
            integer :: n_obs_variables, n_atm_variables
            integer :: i, j, l, x, y, v
            integer :: Mem_Error
            real :: w
            ! prediction period index variables
            integer :: p_start, p_end
            ! training index variables
            integer :: t_tr_start, t_tr_stop, o_tr_start, o_tr_stop, tr_size
            ! variables to store timing data for profiling code
            integer*8 :: timeone, timetwo, master_timeone, master_timetwo, master_time_post_init
            integer*8, dimension(10) :: timers
            integer*8 :: COUNT_RATE

            integer :: total_number_of_gridcells, current_completed_gridcells

            timers = 0
            master_timers = 0
            current_completed_gridcells = 0

            call System_Clock(timeone, COUNT_RATE)
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
            !$omp parallel default(shared)                                                              &
            !$omp private(train_data, pred_data, i, j, v, timeone, timetwo, Mem_Error)                  &
            !$omp firstprivate(nx,ny,n_atm_variables, n_obs_variables, noutput, ntrain, nobs, ntimes)   &
            !$omp firstprivate(p_start, p_end, t_tr_start, t_tr_stop, o_tr_start, o_tr_stop, timers)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Set up the data structures used in downscaling
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            allocate( train_data( ntrain, n_atm_variables+1 ), STAT = Mem_Error )
            if (Mem_Error /= 0) call memory_error(Mem_Error, "train_data", [ntrain, n_atm_variables+1])

            allocate( pred_data(  ntimes, n_atm_variables+1 ), STAT = Mem_Error )
            if (Mem_Error /= 0) call memory_error(Mem_Error, "pred_data", [ntimes, n_atm_variables+1])

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
                do j=1,size(training_atm%variables(v)%data,3)
                    do i=1,size(training_atm%variables(v)%data,2)
                        call transform_data(options%training%input_Xforms(v), training_atm%variables(v)%data(:,i,j), 1, ntrain)
                    enddo
                    if (options%training%input_Xforms(v) /= kNO_TRANSFORM) then
                        call update_statistics(training_atm%variables(v), j)
                    endif

                    if (options%training%normalization_method == kSELF_NORMALIZE) then
                        call normalize(training_atm%variables(v), j)
                    endif
                enddo
                !$omp end do

                !$omp barrier

                !$omp do
                do j=1,size(predictors%variables(v)%data,3)
                    do i=1,size(predictors%variables(v)%data,2)
                        call transform_data(options%prediction%input_Xforms(v), predictors%variables(v)%data(:,i,j), 1, ntimes)
                    enddo
                    if (options%prediction%input_Xforms(v) /= kNO_TRANSFORM) then
                        call update_statistics(predictors%variables(v), j)
                    endif

                    ! does not need to be normalized if it will be transformed to match training_atm anyway
                    if (abs(options%prediction%transformations(v)) /= kQUANTILE_MAPPING) then
                        if (options%prediction%normalization_method == kSELF_NORMALIZE) then
                            call normalize(predictors%variables(v), j)
                        elseif (options%prediction%normalization_method == kTRAININGDATA) then
                            call normalize(predictors%variables(v), j, other=training_atm%variables(v))
                        endif
                    endif
                enddo
                !$omp end do

                !$omp barrier

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
                                                    options, timers, i,j)

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
            ! should be * omp_num_threads(), but if no OPENMP, it won't work...
            master_timers(8) = master_timers(8) + (master_time_post_init - master_timeone) * 16
            master_timers(5) = master_timers(5) + (master_timetwo-master_timeone) * 16 + master_timers(7)
            print*, ""
            print*, "---------------------------------------------------"
            print*, "           Time profiling information "
            print*, "---------------------------------------------------"
            print*, "Total Time : ",  nint(master_timers(5) / real(COUNT_RATE)),                     " s  (CPU time) "
            print*, "---------------------------------------------------"
            print*, "Allocation : ",  nint((100.d0 * master_timers(7)/16.0)/master_timers(5)),  "%    "
            print*, "Data Init  : ",  nint((100.d0 * master_timers(8))/master_timers(5)),  "%    "
            print*, "GeoInterp  : ",  nint((100.d0 * master_timers(4))/master_timers(5)),  "%    "
            print*, "Transform  : ",  nint((100.d0 * master_timers(6))/master_timers(5)),  "%    "
            print*, "Analog     : ",  nint((100.d0 * master_timers(1))/master_timers(5)),  "%    "
            print*, "Regression : ",  nint((100.d0 * master_timers(2))/master_timers(5)),  "%    "
            print*, "Log.Regres : ",  nint((100.d0 * master_timers(3))/master_timers(5)),  "%    "
            print*, "Log.Analog : ",  nint((100.d0 * master_timers(9))/master_timers(5)),  "%    "
            print*, "---------------------------------------------------"
            print*, "Parallelization overhead (assumes the use of 16 processors)"
            print*, "Residual   : ",  100-nint((100.0 * &
                    (sum(master_timers(1:4))+master_timers(6)+sum(master_timers(8:9)))) &
                    / master_timers(5)),"%    "
            print*, "---------------------------------------------------"


    end subroutine downscale

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
        select case (abs(transform_type))

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

    function downscale_point(predictor, atm, obs, errors, output_coeff, logistic, logistic_threshold, options, timers, xptn, ypnt) result(output)
        implicit none
        real,    dimension(:,:), intent(inout):: predictor, atm ! (ntimes, nvars)
        real,    dimension(:),   intent(in)   :: obs
        real,    dimension(:),   intent(inout):: errors
        real,    dimension(:,:), intent(inout):: output_coeff
        real,    dimension(:),   intent(inout):: logistic
        real,                    intent(in)   :: logistic_threshold
        type(config),            intent(in)   :: options
        integer(8),dimension(:), intent(inout):: timers
        integer, intent(in) :: xptn, ypnt

        real,    dimension(:),   allocatable :: output
        real(8), dimension(:),   allocatable :: coefficients
        real,    dimension(:),   allocatable :: coefficients_r4

        integer(8)  :: timeone, timetwo
        integer     :: i, n, nvars, v
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!
        !!  Generic Initialization code
        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nvars = size(atm,2)
        n = size(predictor,1)

        allocate(coefficients(nvars))
        allocate(coefficients_r4(nvars*2))
        allocate(output(n))


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!
        !!  Just pass through a given predictor variable. 
        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (options%pass_through) then
            output = predictor(:, options%pass_through_var + 1)
            return
        endif
            
        ! This just prevents any single points that were WAY out (most likely due to the QM?)
        ! note that if the data have not been normalized, this is not valid
        ! where(predictor < -10) predictor = -10
        ! where(predictor >  10) predictor =  10

            
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!
        !!  Initialization code for pure regression
        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (options%pure_regression) then
            call System_Clock(timeone)
            output(1)  = compute_regression(predictor(1,:), atm, obs, coefficients, errors(1))
            errors(2:) = errors(1)
            where( coefficients >  1e20 ) coefficients =  1e20
            where( coefficients < -1e20 ) coefficients = -1e20
            coefficients_r4(1:nvars) = coefficients
            if (options%debug) then
                do v=1,nvars
                    output_coeff(v,:) = coefficients(v)
                enddo
            endif
            call System_Clock(timetwo)
            timers(2) = timers(2) + (timetwo-timeone)

            call System_Clock(timeone)
            if (logistic_threshold/=kFILL_VALUE) then
                logistic(1) = compute_logistic_regression(predictor(1,:), atm, obs, coefficients, logistic_threshold)
                where( coefficients >  1e20 ) coefficients =  1e20
                where( coefficients < -1e20 ) coefficients = -1e20
                coefficients_r4(nvars+1:nvars*2) = coefficients
                if (options%debug) then
                    do v = 1,nvars
                        output_coeff(v+nvars,:) = coefficients(v)
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
                call downscale_pure_analog(predictor(i,:), atm, obs, coefficients_r4,     &
                                           output(i), errors(i), logistic(i),               &
                                           options, timers, cur_time=i)

                if (options%debug) then
                    do v=1,nvars
                       output_coeff(v,i) = coefficients_r4(v)
                    enddo
                endif

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Analog Regression Downscaling
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (options%analog_regression) then
                call downscale_analog_regression(predictor(i,:), atm, obs, coefficients_r4,     &
                                                 output(i), errors(i), logistic(i),             &
                                                 options, timers, [xptn,ypnt,i], cur_time=i)

                if (options%debug) then
                    do v=1,nvars
                        output_coeff(v,i) = coefficients_r4(v)
                    enddo
                endif


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!
            !!  Pure regression downscaling
            !!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            elseif (options%pure_regression) then
                ! to test matmul(predictor, output_coeff(:,1)) should provide all time values efficiently?
                call apply_pure_regression(output(i), predictor(i,:), coefficients_r4, logistic(i),  logistic_threshold, timers)

            endif
        end do

    end function downscale_point

    subroutine downscale_analog_regression(x, atm, obs, output_coeff, output, error, logistic, options, timers, cur_point, cur_time)
        implicit none
        real,           intent(in),     dimension(:,:)  :: atm
        real,           intent(in),     dimension(:)    :: x, obs
        real,           intent(inout),  dimension(:)    :: output_coeff
        real,           intent(inout)                   :: output, logistic, error
        type(config),   intent(in)                      :: options
        integer*8,      intent(inout),  dimension(:)    :: timers

        integer,        intent(in),     dimension(3)    :: cur_point
        integer,        intent(in),     optional        :: cur_time
        integer     :: test

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
        real,    dimension(:),  allocatable :: weights, packed_weights
        integer :: n_packed

        n_analogs           = options%n_analogs
        analog_threshold    = options%analog_threshold
        logistic_threshold  = options%logistic_threshold
        nvars               = size(x)

        call System_Clock(timeone)
        if (options%analog_weights) then
            call find_analogs(analogs, x, atm, n_analogs, analog_threshold, weights, skip_analog=cur_time)
        else
            call find_analogs(analogs, x, atm, n_analogs, analog_threshold, skip_analog=cur_time)
        endif

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
            output = compute_regression(x, regression_data, obs_analogs, coefficients, error, weights)
        else

            threshold_packing = obs_analogs > logistic_threshold
            n_packed = count(threshold_packing)

            if (n_packed > (nvars*2)) then
                allocate(threshold_atm(n_packed, nvars))
                allocate(threshold_obs(n_packed))
                do j=1,nvars
                    threshold_atm(:,j) = pack(regression_data(:,j), threshold_packing)
                enddo
                threshold_obs = pack(obs_analogs, threshold_packing)
                if (options%analog_weights) then
                    allocate(packed_weights(n_packed))
                    packed_weights = pack(weights, threshold_packing)
                endif

                output = compute_regression(x, threshold_atm, threshold_obs, coefficients, error, packed_weights)

                if ((output>maxval(threshold_obs)*1.2)                        &
                    .or.(output<(logistic_threshold-2))                       &
                    .or.(abs(coefficients(1))>(maxval(threshold_obs) * 1.5))) then

                    call System_Clock(timetwo)
                    timers(2) = timers(2) + (timetwo-timeone)
                    ! revert to a pure analog approach.  By passing analogs and weights, it will not recompute which analogs to use
                    ! it will just compute the analog mean, pop, and error statistics
                    if (present(cur_time)) then
                        call downscale_pure_analog(x, atm, obs, output_coeff, output, error, logistic, options, timers, analogs, weights, cur_time)
                    else
                        call downscale_pure_analog(x, atm, obs, output_coeff, output, error, logistic, options, timers, analogs, weights)
                    endif
                    coefficients(1:nvars) = output_coeff(1:nvars)
                    coefficients(1) = 1e20
                    call System_Clock(timeone)
                endif
                ! if the regression picked a vaguely valid value <0 then just set it to 0?
                ! But since this is just the mean expected value, it can have a positive error term added to it resulting
                ! in a net positive value, so for now we will accept "negative" precipitation amounts
                ! if (output<0) then
                !     output=0
                ! endif

            elseif (n_packed > 0) then
                output = sum(pack(obs_analogs, threshold_packing)) / n_packed
            else
                ! note, for precip (and a 0 threshold), this could just be setting output, error, logistic, and coefficients to 0
                output = compute_regression(x, regression_data, obs_analogs, coefficients, error, weights)
            endif

        endif

        where( coefficients >  1e20 ) coefficients =  1e20
        where( coefficients < -1e20 ) coefficients = -1e20
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

                    if (options%analog_weights) then
                        logistic = compute_analog_exceedance(obs, analogs, logistic_threshold, weights)
                    else
                        ! compute the logistic as the probability of threshold exceedance in the analog population
                        logistic = compute_analog_exceedance(obs, analogs, logistic_threshold)
                    endif
                endif

            else
                logistic = compute_logistic_regression(x, regression_data, obs_analogs, coefficients, logistic_threshold)
                where( coefficients >  1e20 ) coefficients =  1e20
                where( coefficients < -1e20 ) coefficients = -1e20
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

    subroutine downscale_pure_analog(x, atm, obs, output_coeff, output, error, logistic, options, timers, input_analogs, input_weights, cur_time)
        implicit none
        real,           intent(in),     dimension(:,:)  :: atm
        real,           intent(in),     dimension(:)    :: x, obs
        real,           intent(inout),  dimension(:)    :: output_coeff
        real,           intent(inout)                   :: output, logistic, error
        type(config),   intent(in)                      :: options
        integer*8,      intent(inout),  dimension(:)    :: timers
        integer,        intent(in),     dimension(:),   optional    :: input_analogs
        real,           intent(in),     dimension(:),   optional    :: input_weights
        integer,        intent(in),     optional        :: cur_time

        real        :: analog_threshold, logistic_threshold
        integer*8   :: timeone, timetwo
        integer     :: real_analogs, selected_analog, nvars, n_analogs
        integer, dimension(:), allocatable :: analogs
        real,    dimension(:), allocatable :: weights
        real        :: rand
        integer     :: j

        n_analogs           = options%n_analogs
        analog_threshold    = options%analog_threshold
        logistic_threshold  = options%logistic_threshold
        nvars               = size(x)

        call System_Clock(timeone)
        ! find the best n_analog matching analog days in atm to match x
        if ((.not.present(input_analogs)).or.((options%analog_weights).and.(.not.present(input_weights)))) then
            if (n_analogs > 0) then
                allocate(analogs(n_analogs))
            endif
            if (options%analog_weights) then
                call find_analogs(analogs, x, atm, n_analogs, analog_threshold, weights, skip_analog=cur_time)
            else
                call find_analogs(analogs, x, atm, n_analogs, analog_threshold, skip_analog=cur_time)
            endif
        else
            allocate(analogs(size(input_analogs)))
            analogs = input_analogs
            if (options%analog_weights) then
                allocate(weights(size(input_weights)))
                weights = input_weights
            endif
        endif

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
            if (options%analog_weights) then
                output = compute_analog_mean(obs, analogs, weights)
            else
                output = compute_analog_mean(obs, analogs)
            endif

            if (options%debug) then
                do j=1,nvars
                    output_coeff(j) = compute_analog_mean(atm(:, j), analogs)
                enddo
            endif
        endif

        if (options%analog_weights) then
            error = compute_analog_error(obs, analogs, output, weights)
        else
            error = compute_analog_error(obs, analogs, output)
        endif

        call System_Clock(timetwo)
        timers(2) = timers(2) + (timetwo-timeone)

        call System_Clock(timeone)
        if (logistic_threshold/=kFILL_VALUE) then
            if (options%analog_weights) then
                logistic = compute_analog_exceedance(obs, analogs, logistic_threshold, weights)
            else
                logistic = compute_analog_exceedance(obs, analogs, logistic_threshold)
            endif
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
    subroutine normalize(var, j, other)
        implicit none
        class(atm_variable_type), target, intent(inout) :: var
        integer, intent(in)                     :: j
        class(atm_variable_type), intent(in), target, optional :: other
        integer :: i, n
        class(atm_variable_type), pointer :: norm_data

        n = size(var%mean,1)

        if (present(other)) then
            norm_data => other
        else
            norm_data => var
        endif

        do i=1,n
            var%data(:,i,j) = var%data(:,i,j) - norm_data%mean(i,j)

            if (norm_data%stddev(i,j) /= 0) then
                var%data(:,i,j) = var%data(:,i,j) / norm_data%stddev(i,j)
            else
                if (maxval(abs(var%data(:,i,j))) > 0) then
                    !$omp critical (print_lock)
                    write(*,*) "ERROR Normalizing:", trim(var%name)
                    write(*,*) "  For point: ", i, j
                    write(*,*) "  Should this point have been masked?", var%data(1,i,j)
                    !$omp end critical (print_lock)
                endif
            endif
            ! shift to a 0-based range so that variables such as precip have a testable non-value
            var%min_val(i,j) = minval(var%data(:,i,j))

            var%data(:,i,j) = var%data(:,i,j) - norm_data%min_val(i,j)
            where(abs(var%data(:,i,j)) < 1e-10) var%data(:,i,j)=0
        enddo


    end subroutine normalize

    subroutine update_statistics(input_var, j)
        implicit none
        class(variable_type), intent(inout) :: input_var
        integer, intent(in) :: j

        integer :: ntimes, nx, i

        ntimes = size(input_var%data,1)
        nx     = size(input_var%data,2)

        do i=1,nx
            input_var%mean(i,j) = sum(input_var%data(:,i,j)) / ntimes
            input_var%stddev(i,j) = stddev(input_var%data(:,i,j), input_var%mean(i,j))
            input_var%min_val(i,j) = minval(input_var%data(:,i,j))
        enddo

    end subroutine update_statistics


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

        integer :: v, Mem_Error

        write(*,*), "Allocating Memory"
        write(*,*), "  N output times:", noutput
        write(*,*), "  nx:", nx, "        ny:",ny
        write(*,*), "  Training size:",  tr_size
        write(*,*), "  N atm variables:",n_atm_variables

        allocate(output%variables(n_obs_variables), stat=Mem_Error)
        if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables", [n_obs_variables])

        do v = 1,n_obs_variables
            allocate(output%variables(v)%data        (noutput, nx, ny), stat=Mem_Error)
            if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%data v="//trim(str(v)), [noutput,nx,ny])
            allocate(output%variables(v)%errors      (noutput, nx, ny), stat=Mem_Error)
            if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%errors v="//trim(str(v)), [noutput,nx,ny])

            if (options%logistic_threshold/=kFILL_VALUE) then
                allocate(output%variables(v)%logistic    (noutput, nx, ny), stat=Mem_Error)
                if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%logistic v="//trim(str(v)), [noutput,nx,ny])
            else
                ! unfortunately logistic has to be allocated even if it is not used because it is indexed
                ! however, it does not need a complete time sequence because time is never indexed
                allocate(output%variables(v)%logistic    (2, nx, ny), stat=Mem_Error)
                if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%logistic v="//trim(str(v)), [2,nx,ny])
            endif

            output%variables(v)%data            = 0
            output%variables(v)%logistic        = 0

            output%variables(v)%logistic_threshold = options%logistic_threshold

            if (options%obs%input_Xforms(v) == kLOG_TRANSFORM) then
                if (options%logistic_threshold<=0) then
                    output%variables(v)%logistic_threshold = log(LOG_FILL_VALUE)
                endif
            endif

            if (options%debug) then
                allocate(output%variables(v)%obs         (tr_size, nx, ny), stat=Mem_Error)
                if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%obs v="//trim(str(v)), [tr_size,nx,ny])
                allocate(output%variables(v)%training    (tr_size, nx, ny, n_atm_variables+1), stat=Mem_Error)
                if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%training v="//trim(str(v)), [tr_size,nx,ny, n_atm_variables+1])
                allocate(output%variables(v)%predictors  (noutput, nx, ny, n_atm_variables+1), stat=Mem_Error)
                if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%predictors v="//trim(str(v)), [noutput,nx,ny, n_atm_variables+1])

                if (options%logistic_threshold/=kFILL_VALUE) then
                    allocate(output%variables(v)%coefficients((n_atm_variables+1)*2, noutput, nx, ny), stat=Mem_Error)
                    if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%coefficients v="//trim(str(v)), [(n_atm_variables+1)*2,noutput,nx,ny])
                else
                    allocate(output%variables(v)%coefficients(n_atm_variables+1, noutput, nx, ny), stat=Mem_Error)
                    if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%coefficients v="//trim(str(v)), [n_atm_variables+1,noutput,nx,ny])
                endif

                output%variables(v)%predictors   = 0
                output%variables(v)%training     = 0
                output%variables(v)%obs          = 0
                output%variables(v)%coefficients = 0
            else
                ! unfortunately coefficients has to be allocated even if it is not used because it is indexed
                ! however, it does not need a complete time sequence or variable size because they are never indexed
                allocate(output%variables(v)%coefficients(n_atm_variables+1, 2, nx, ny), stat=Mem_Error)
                if (Mem_Error /= 0) call memory_error(Mem_Error, "out%variables(v)%coefficients v="//trim(str(v)), [2,2,nx,ny])
            endif

        end do
    end subroutine allocate_data

    subroutine memory_error(error, variable_name, dims)
        implicit none
        integer,          intent(in)                :: error
        character(len=*), intent(in)                :: variable_name
        integer,          intent(in), dimension(:)  :: dims

        write(*,*), "Error allocating memory for variable: ", trim(variable_name)
        write(*,*), "  ERROR        = ", error
        write(*,*), "  Dimensions   = ", dims

        stop "MEMORY ALLOCATION ERROR"

    end subroutine memory_error


end module downscaling_mod
