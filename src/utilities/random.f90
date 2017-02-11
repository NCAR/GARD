module random_mod

    use model_constants

    implicit none

    double precision, parameter :: root_pi_over_8 = 0.626657068658 ! sqrt(kPI/8)
    
contains
    
    !>------------------------------------------------
    !! Use the Box-Muller Transform to convert uniform to normal random deviates
    !!
    !! Note random_sample should be an allocated 1D real array
    !! On return, random_sample will be filled with random normal (0,1) data
    !!
    !! Caveat: this transform is unable to generate extremely high values (>6.66)
    !! Having switched to double precision it may do better now.
    !!
    !-------------------------------------------------
    subroutine box_muller_random(random_sample)
        implicit none
        real, intent(inout) :: random_sample(:)
        integer :: n,i

        double precision :: u1, u2, s
        double precision, allocatable :: double_random(:)

        n = size(random_sample)
        allocate(double_random(n))
        call random_number(double_random)

        do i=1,n,2
            u1 = double_random(i)
            if (i<n) then
                u2 = double_random(i+1)
            else
                call random_number(u2)
            endif

            s = sqrt(-2 * log(u1))
            random_sample(i) = s * cos(2 * kPI * u2)
            if (i<n) then
                random_sample(i+1) = s * sin(2 * kPI * u2)
            endif

        enddo

    end subroutine box_muller_random
    
    !>------------------------------------------------
    !! Computes an approximation to the normal CDF based on Aludaat and Alodat (2008)
    !!
    !! Aludaat, K.M. and Alodat, M.T. (2008). A note on approximating the normal distribution function. Applied Mathematical Sciences, Vol 2, no 9, pgs 425-429.
    !!
    !!------------------------------------------------
    function normal_cdf(normal) result(cdf)
        implicit none
        real :: normal
        real :: cdf
        
        cdf = 0.5 * (1+ sqrt((1- exp(0-(root_pi_over_8 * normal**2))) ))
        
        if (normal < 0) then
            cdf = 1 - cdf
        endif
        
    end function normal_cdf

    

    subroutine sample_distribution(output, mean_values, error_term, uniform_random_values, exceedence_probability, threshold)
        implicit none
        real, intent(out) :: output(:)
        real, intent(in)  :: mean_values(:), error_term(:), uniform_random_values(:)
        real, intent(in), optional :: exceedence_probability(:)
        real, intent(in), optional :: threshold
        
        integer :: i,n
        real :: random_normal, poe, threshold_internal
        
        threshold_internal = 0
        if (present(threshold)) threshold_internal = threshold
        
        n = size(mean_values)
        
        do i=1,n
            if (present(exceedence_probability)) then
                if (exceedence_probability(i) > (1-uniform_random_values(i))) then

                    ! random_normal = get_random_from_cdf(uniform_random_values(i) - exceedence_probability(i) / (1-exceedence_probability(i)) )

                endif
            endif
            
        enddo
        output = 1
    end subroutine sample_distribution

end module random_mod
