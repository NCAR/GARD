module downscaling_mod

    use data_structures
    use string,             only : str
    use regression_mod,     only : compute_regression, compute_logistic_regression
    use analog_mod,         only : find_analogs, compute_analog_mean, compute_analog_error, compute_analog_exceedance
    use model_constants
    use quantile_mapping,   only : develop_qm, apply_qm
    use time_util,          only : setup_time_indices
    implicit none
    
    double precision, dimension(10) :: master_timers

    
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
            double precision :: timeone, timetwo, master_timeone, master_timetwo
            double precision, dimension(10) :: timers
            
            integer :: total_number_of_gridcells, current_completed_gridcells
            
            timers = 0
            master_timers = 0
            current_completed_gridcells = 0
            
            call CPU_TIME(timeone)
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
            
            nx = 32
            ny = 132
            total_number_of_gridcells = nx * (ny-100)
            
            call allocate_data(output, n_obs_variables, noutput, nx, ny, n_atm_variables, tr_size, options)

            call CPU_TIME(master_timeone)
            write(*,*), "Entering Parallel Region and Normalizing"
            !$omp parallel default(shared)                  &
            !$omp private(train_data, pred_data, i, j, v, timeone, timetwo)   &
            !$omp firstprivate(nx,ny,n_atm_variables, n_obs_variables, noutput, ntrain, nobs, ntimes) &
            !$omp firstprivate(p_start, p_end, t_tr_start, t_tr_stop, o_tr_start, o_tr_stop, timers)
            
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
            !$omp single
            write(*,*), "Downscaling..."
            !$omp end single
            ! parallelization could be over x and y, (do n=1,ny*nx; j=n/nx; i=mod(n,nx)) and use schedule(dynamic)
            !$omp do schedule(static, 1)
            do j=101,ny
                do i=1,nx
                    
                    if (training_obs%mask(i,j)) then
                        
                        ! index 1 is the constant for the regression code... should be moved into regression code...?
                        output%variables(1)%predictors(:,i,j,1) = pred_data( p_start   : p_end,    1)
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
                                !$omp critical (print_lock)
                                write(*,*) ""
                                write(*,*) "-----------------------------"
                                write(*,*) "   Downscaling point: ",i,j
                                write(*,*) "   For variable     : ", trim(training_obs%variables(v)%name)
                                !$omp end critical (print_lock)
                            endif
                            
                            call CPU_TIME(timeone)
                            ! As tempting as it may be, associate statements are not threadsafe!!!
                            output%variables(v)%data(:,i,j) = downscale_point(                                         &
                                                    pred_data                       (   p_start : p_end,     :),       &
                                                    train_data                      (t_tr_start : t_tr_stop, :),       &
                                                    training_obs%variables(v)%data  (o_tr_start : o_tr_stop, i, j),    &
                                                    output%variables(v)%errors      (           :,           i, j),    &
                                                    output%variables(v)%coefficients(           :,  :,       i, j),    &
                                                    output%variables(v)%logistic    (           :,           i, j),    &
                                                    output%variables(v)%logistic_threshold,                            &
                                                    options, timeone, timers )
                                                          
                            output%variables(v)%obs(:,i,j) = training_obs%variables(v)%data(o_tr_start:o_tr_stop, i, j)
                            call CPU_TIME(timetwo)
                            !$omp critical
                            timers(8) = timers(8) + (timetwo-timeone)
                            !$omp end critical

                        enddo
                    else ! this is a masked point in the observations
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
                    !$omp critical (print_lock)
                    current_completed_gridcells = current_completed_gridcells + 1
                    ! write(*,"(A,f5.1,A)", advance="NO") char(13), current_completed_gridcells/real(total_number_of_gridcells)*100.0, " %"
                    write(*,"(A,f5.1,A$)") char(13), current_completed_gridcells/real(total_number_of_gridcells)*100.0, " %"
                    !$omp end critical (print_lock)
                enddo
            enddo
            !$omp end do
            !$omp critical
            master_timers = timers + master_timers
            !$omp end critical
            !$omp end parallel
            call CPU_TIME(master_timetwo)
            master_timers(5) = master_timers(5) + (master_timetwo-master_timeone)*16
            print*, ""
            print*, "---------------------------------------------------"
            print*, "           Time profiling information "
            print*, "---------------------------------------------------"
            print*, "Total Time : ",  nint(100.0),                     "%    ",  nint(master_timers(5))
            print*, "---------------------------------------------------"
            print*, "Total DSpt : ",  nint(100 * master_timers(8)/master_timers(5)),  "%    ", nint(master_timers(8))
            print*, "Internal   : ",  nint(100 * master_timers(4)/master_timers(5)),  "%    ", nint(master_timers(4))
            print*, "Analog     : ",  nint(100 * master_timers(1)/master_timers(5)),  "%    ", nint(master_timers(1))
            print*, "Regression : ",  nint(100 * master_timers(2)/master_timers(5)),  "%    ", nint(master_timers(2))
            print*, "Log.Regres : ",  nint(100 * master_timers(3)/master_timers(5)),  "%    ", nint(master_timers(3))
            print*, "Second 1/2 : ",  nint(100 * master_timers(7)/master_timers(5)),  "%    ", nint(master_timers(7))
            print*, "Sec.1/2 in : ",  nint(100 * master_timers(6)/master_timers(5)),  "%    ", nint(master_timers(6))
            print*, "Log.Analog : ",  nint(100 * master_timers(9)/master_timers(5)),  "%    ", nint(master_timers(9))
            print*, "Setup time : ",  nint(100 * master_timers(10)/master_timers(5)), "%    ", nint(master_timers(10))

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
    
    function downscale_point(predictor, atm, obs, errors, output_coeff, logistic, logistic_threshold, options, start_time, timers) result(output)
        implicit none
        real,    dimension(:,:), intent(in)   :: predictor, atm ! (ntimes, nvars)
        real,    dimension(:),   intent(in)   :: obs
        real,    dimension(:),   intent(inout):: errors
        real,    dimension(:,:), intent(inout):: output_coeff
        real,    dimension(:),   intent(inout):: logistic
        real,                    intent(in)   :: logistic_threshold
        type(config),            intent(in)   :: options
        double precision,        intent(in)   :: start_time
        double precision, dimension(:), intent(inout) :: timers
        
        real,    dimension(:),   allocatable :: output
        real,    dimension(:),   allocatable :: obs_analogs
        real,    dimension(:,:), allocatable :: regression_data
        integer, dimension(:),   allocatable :: analogs
        real(8), dimension(:),   allocatable :: coefficients
        
        double precision :: timeone, timetwo
        double precision :: master_timeone, master_timetwo, half_time, inner_timeone
        integer :: i, j, n, nvars
        integer :: a, n_analogs, selected_analog
        real    :: rand
        integer :: omp_get_thread_num
        
        call CPU_TIME(master_timeone)
        n_analogs = options%n_analogs
        nvars = size(atm,2)
        n = size(predictor,1)
        
        allocate(coefficients(nvars))
        allocate(output(n))
        allocate(analogs(n_analogs))
        
        call CPU_TIME(timeone)
        if (options%analog_regression) then
            allocate(regression_data(n_analogs, nvars))
            allocate(obs_analogs(n_analogs))
        endif
        call CPU_TIME(timetwo)
        timers(1) = timers(1) + (timetwo-timeone)
        timers(10) = timers(10) + (timetwo-start_time)

        if (options%pure_regression) then
            call CPU_TIME(timeone)
            output(1)  = compute_regression(predictor(1,:), atm, obs, coefficients, errors(1))
            errors(2:) = errors(1)
            do i=1,nvars
                output_coeff(i,:) = coefficients(i)
            enddo
            call CPU_TIME(timetwo)
            timers(2) = timers(2) + (timetwo-timeone)

            call CPU_TIME(timeone)
            if (logistic_threshold/=kFILL_VALUE) then
                logistic(1) = compute_logistic_regression(predictor(1,:), atm, obs, coefficients, logistic_threshold)
                
                do i = 1,nvars
                    output_coeff(i+nvars,:) = coefficients(i)
                enddo
            endif
            call CPU_TIME(timetwo)
            timers(3) = timers(3) + (timetwo-timeone)

        endif
        
        call CPU_TIME(half_time)
        do i = 1, n
            if (options%pure_analog) then
                
                analogs = find_analogs(predictor(i,:), atm, n_analogs)
                
                if (options%sample_analog) then
                    call random_number(rand)
                    selected_analog = floor(rand * n_analogs)+1
                    
                    output(i) = obs( analogs(selected_analog) )
                else
                    output(i) = compute_analog_mean(obs, analogs)
                endif
                
                errors(i) = compute_analog_error(obs, analogs, output(i))
                
                output_coeff(1:nvars,i) = atm(analogs(selected_analog), :)

                if (logistic_threshold/=kFILL_VALUE) then
                    logistic(i) = compute_analog_exceedance(obs, analogs, logistic_threshold)
                endif
                
            elseif (options%analog_regression) then
                call CPU_TIME(inner_timeone)
                call CPU_TIME(timeone)
                analogs = find_analogs(predictor(i,:), atm, n_analogs)
                do a=1,n_analogs
                    obs_analogs(a) = obs(analogs(a))
                    regression_data(a,:) = atm(analogs(a),:)
                enddo
                call CPU_TIME(timetwo)
                timers(1) = timers(1) + (timetwo-timeone)

                call CPU_TIME(timeone)
                output(i) = compute_regression(predictor(i,:), regression_data, obs_analogs, coefficients, errors(i))
                output_coeff(1:nvars,i) = coefficients
                call CPU_TIME(timetwo)
                timers(2) = timers(2) + (timetwo-timeone)

                call CPU_TIME(timeone)
                if (logistic_threshold/=kFILL_VALUE) then
                    if (options%logistic_from_analog_exceedance) then
                        
                        ! if the user specified using fewer analogs specifically for the logistic component, find them here
                        if (options%n_log_analogs /= n_analogs) then
                            analogs(1:options%n_log_analogs) = find_analogs(predictor(i,:), atm, options%n_log_analogs)
                            call CPU_TIME(timetwo)
                            timers(9) = timers(9) + (timetwo-timeone)
                            
                            call CPU_TIME(timeone)
                            logistic(i) = compute_analog_exceedance(obs, analogs(1:options%n_log_analogs), logistic_threshold)
                        else
                            logistic(i) = compute_analog_exceedance(obs, analogs, logistic_threshold)
                        endif
                            
                        do j = 1,nvars
                            output_coeff(j+nvars,i) = analogs(j)
                        enddo
                    else
                        logistic(i) = compute_logistic_regression(predictor(i,:), regression_data, obs_analogs, coefficients, logistic_threshold)
                        do j = 1,nvars
                            output_coeff(j+nvars,i) = coefficients(j)
                        enddo
                    endif
                endif
                call CPU_TIME(timetwo)
                timers(3) = timers(3) + (timetwo-timeone)
                
                call CPU_TIME(timetwo)
                timers(6) = timers(6) + (timetwo - inner_timeone)
                if (omp_get_thread_num() == 2) then
                !$omp critical (print_lock)
                    print*, timers(6)-(timers(1)+timers(2)+timers(3)+timers(9)), timers(6), timers(1)+timers(2)+timers(3)+timers(9)
                    print*, timers(1), timers(2), timers(3), timers(9)
                !$omp end critical (print_lock)
                endif
            elseif (options%pure_regression) then
                !  test matmul(predictor, output_coeff(:,1)) should provide all time values efficiently?
                ! output(i) = dot_product(predictor(i,:), output_coeff(:,1))
                output(i) = 0
                do j=1,nvars
                    output(i) = output(i) + predictor(i,j) * output_coeff(j,1)
                enddo

                if (logistic_threshold/=kFILL_VALUE) then
                    logistic(i) = 1.0 / (1.0 + exp(-dot_product(predictor(i,:), output_coeff(nvars+1:nvars*2,1))))
                endif
                
            endif
        end do
        call CPU_TIME(master_timetwo)
        timers(4) = timers(4) + (master_timetwo - master_timeone)
        timers(7) = timers(7) + (master_timetwo - half_time)

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
            
            allocate(output%variables(v)%obs         (tr_size, nx, ny))
            allocate(output%variables(v)%training    (tr_size, nx, ny, n_atm_variables+1))
            allocate(output%variables(v)%predictors  (noutput, nx, ny, n_atm_variables+1))
            
            
            output%variables(v)%logistic_threshold = options%logistic_threshold
            output%variables(v)%data        = 0
            output%variables(v)%predictors  = 0
            output%variables(v)%training    = 0
            output%variables(v)%obs         = 0
        end do

        
        
    end subroutine allocate_data

end module downscaling_mod
