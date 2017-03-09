module sampling_mod

    use random_mod, only : normal_cdf, box_muller_random, get_normal_from_cdf
    implicit none
contains

    subroutine sample_distribution(output, mean_values, error_term, &
                                   normal_random_values, uniform_random_values, &
                                   exceedence_probability, threshold)
        implicit none
        real, intent(out) :: output(:)
        real, intent(in)  :: mean_values(:), error_term(:)
        real, intent(in), optional :: normal_random_values(:)
        real, intent(in), optional :: exceedence_probability(:)
        real, intent(in), optional :: uniform_random_values(:)
        real, intent(in), optional :: threshold

        integer :: i,n
        real :: random_normal, poe, threshold_internal
        real, allocatable :: uniform(:), normal(:)

        threshold_internal = 0
        if (present(threshold)) threshold_internal = threshold

        n = size(mean_values)
        ! the following section just ensures that optional variables have the necessary internal variables created.
        if (present(exceedence_probability)) then
            allocate(uniform(n))
            if (present(uniform_random_values)) then
                uniform = uniform_random_values
            else
                if (present(normal_random_values)) then
                    ! if we were given normal random values, but not uniform, then create uniform from normal
                    do i=1,n
                        uniform(i) = normal_cdf( normal_random_values(i) )
                    enddo
                else
                    ! While this will generate n uniform random numbers, there will not be any spatio-temporal
                    ! auto-correlation to them.
                    call random_number(uniform)
                endif
            endif
        else
            allocate(normal(n))
            if (present(normal_random_values)) then
                normal = normal_random_values
            else
                if (present(uniform_random_values)) then
                    normal = uniform_random_values
                else
                    call random_number(normal)
                endif
                call box_muller_random(normal)
            endif
        endif

        ! this is the heart of the distribution sampling code
        do i=1,n
            if (present(exceedence_probability)) then
                if (uniform(i) > (1-exceedence_probability(i)) ) then

                    ! if random number is greater than the non-exceedence probability, then we have exceeded it...
                    ! rescale the uniform random number for the region of exceedence and convert to a normal random number
                    random_normal = get_normal_from_cdf((uniform(i) - (1 - exceedence_probability(i))) &
                                    / (exceedence_probability(i)))

                    ! ultimately this is the same equation as below for variables without an exceedence threshold
                    output(i) = mean_values(i) + error_term(i) * random_normal

                    ! also check, if the random term has pushed the value below the threshold anyway, then just set to the threshold
                    if (output(i) < threshold_internal) then
                        output(i) = threshold_internal
                    endif
                else
                    ! if this value did not exceed the threshold, then just set it to the threshold.
                    output(i) = threshold_internal
                endif
            else
                ! standard approach is to take the residuals multiply by a normal random variate and add to the mean prediction
                output(i) = mean_values(i) + error_term(i) * normal(i)
            endif
        enddo

    end subroutine sample_distribution


end module sampling_mod
