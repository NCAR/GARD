module downscaling_mod

    use data_structures
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
            real, dimension(:), allocatable :: train_data, pred_data
            integer :: nx, ny, ntimes, n_variables, ntrain, nobs
            integer :: i, j, l, x, y
            real :: w
            
            nx = 7
            ny = 1
            ntimes = size(predictors%variables(1)%data,1)
            ntrain = size(training_atm%variables(1)%data,1)
            nobs   = size(training_obs%variables(1)%data,1)
            n_variables = size(training_obs%variables)
            
            print*, "ntrain=",ntrain, "ntimes=",ntimes, "n_variables=",n_variables
            
            allocate( train_data( ntrain ) )
            allocate( pred_data(  ntimes ) )
            allocate(output%variables(n_variables))
            
            call setup_timing(training_atm, training_obs, predictors, options)
            
            do i = 1, n_variables
                allocate(output%variables(i)%data(ntimes, nx, ny))
                
                associate(var        => output%variables(i)%data,       &
                          predict    => predictors%variables(i)%data,   &
                          train      => training_atm%variables(i)%data, &
                          observed   => training_obs%variables(1)%data, &
                          p_geo      => predictors%geoLUT,              &
                          t_geo      => training_atm%geoLUT,            &
                          p_start    => predictors%first_time,          &
                          p_end      => predictors%last_time,           &
                          t_tr_start => training_atm%training_start,    &
                          t_tr_stop  => training_atm%training_stop,     &
                          o_tr_start => training_obs%training_start,    &
                          o_tr_stop  => training_obs%training_stop,     &
                          p_xf_start => predictors%transform_start,     &
                          p_xf_stop  => predictors%transform_stop,      &
                          t_xf_start => training_atm%transform_start,   &
                          t_xf_stop  => training_atm%transform_stop     &
                          )
                
                
                ! get interpolated location (via various possible interpolation schemes...)
                pred_data  = 0
                train_data = 0
                do l = 1,4
                    x = p_geo%x(l,135,100)
                    y = p_geo%y(l,135,100)
                    w = p_geo%w(l,135,100)
                    pred_data = pred_data + predict(:, x, y) * w
                    x = t_geo%x(l,135,100)
                    y = t_geo%y(l,135,100)
                    w = t_geo%w(l,135,100)
                    train_data = train_data + train(:, x, y) * w
                end do

                print*, shape(observed)
                ! transform the relevant data (e.g. quantile map or other)
                call develop_qm(pred_data( p_xf_start:p_xf_stop), &
                                train_data(t_xf_start:t_xf_stop), qm, n_segments = 300)

                call apply_qm(pred_data, var(:,2,1), qm)
                
                var(:,1,1) = pred_data
                var(:ntrain,3,1) = train_data
                var(:,4,1) = pred_data
                var(:nobs,5,1) = observed(:,135,100)
                
                ! compute regressions as necessary
                print*, ""
                print*, "-----------------------------"
                print*, "Downscaling point: ",135, 100
                print*, "Predictor range: ", p_start, p_end
                print*, "Observed range: ", o_tr_start, o_tr_stop, o_tr_stop-o_tr_start
                print*, "Training range: ", t_tr_start, t_tr_stop, t_tr_stop-t_tr_start
                var(p_start:p_end,6,1) = downscale_point( var(p_start:p_end,2,1), train_data(t_tr_start:t_tr_stop), observed(o_tr_start:o_tr_stop, 135,100), options)
                
                call develop_qm(var(p_start:p_end,6,1), &
                                observed(o_tr_start:o_tr_stop, 135,100), qm, n_segments = 300)

                call apply_qm(var(p_start:p_end,6,1), var(p_start:p_end,7,1), qm)

                
                end associate
            end do
            
            
    end function downscale
    
    function downscale_point(predictor, atm, obs, options) result(output)
        implicit none
        real,    dimension(:), intent(in)  :: predictor, atm, obs
        type(config),          intent(in)  :: options
        real,    dimension(:), allocatable :: output
        integer, dimension(:), allocatable :: analogs
        
        integer :: i,n, n_analogs, selected_analog
        real    :: rand
        
        n_analogs = options%n_analogs
        n = size(predictor)
        allocate(output(n))
        allocate(analogs(n_analogs))
        
        do i = 1, n
            analogs = find_analogs(predictor(i), atm, n_analogs)
            call random_number(rand)
            selected_analog = floor(rand * n_analogs)+1
            print*, analogs(selected_analog)
            output(i) = obs( analogs(selected_analog) )
        end do
        
    end function downscale_point
    
    function find_analogs(match, input, n) result(analogs)
        implicit none
        real,    intent(in)     :: match
        real,    intent(in), dimension(:) :: input
        integer, intent(in)     :: n
        integer, dimension(1:n)   :: analogs
        real,    dimension(:), allocatable :: distances
        logical, dimension(:), allocatable :: mask
        integer, dimension(1) :: min_location
        integer :: i, n_inputs
        real    :: large_value
        
        if (match==0) then
            analogs = pick_n_random_zeros(input, n)
        else
            
            n_inputs = size(input)
            allocate(mask(n_inputs))
            allocate(distances(n_inputs))
            
            mask = .True.
            distances = abs(input - match)
            
            large_value = maxval(distances)
            do i = 1, n
                min_location = minloc(distances,mask=mask)
                mask(min_location) = .false.
                
                analogs(i) = min_location(1)
            end do
        endif
        
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
