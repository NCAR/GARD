module downscaling_mod

    use data_structures
    use quantile_mapping,   only : develop_qm, apply_qm
    use sort_mod,           only : sort
    implicit none
    
contains
    function downscale(training_atm, training_obs, predictors, options) result(output)
            implicit none
            type(atm),    intent(in) :: training_atm, predictors
            type(obs),    intent(in) :: training_obs
            type(config), intent(in) :: options
            
            type(results) :: output
            type(qm_correction_type) :: qm
            integer :: nx, ny, ntimes, n_variables
            integer :: i, j
            
            nx = 3
            ny = 1
            ntimes = size(predictors%variables(1)%data,1)
            n_variables = size(training_obs%variables)
            
            allocate(output%variables(n_variables))
            
            do i = 1, n_variables
                associate(var => output%variables(i))
                
                allocate(var%data(ntimes, nx, ny))
                var%data(:,1,1) = predictors%variables(1)%data(:,1,1)
                
                call develop_qm(predictors%variables(1)%data(:,1,1), training_atm%variables(1)%data(:,35,15), qm, n_segments = 3000)

                call apply_qm(predictors%variables(1)%data(:,1,1), var%data(:,2,1), qm)
                var%data(:size(training_atm%variables(1)%data,1),3,1) = training_atm%variables(1)%data(:,35,15)
                
                
                
                end associate
            end do
            
            
    end function downscale

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
