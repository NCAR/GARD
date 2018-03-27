module sampling_mod

    implicit none
interface

    module subroutine sample_distribution(output, mean_values, error_term, &
                                   normal_random_values, uniform_random_values, &
                                   exceedence_probability, threshold)
        implicit none
        real, intent(out) :: output(:)
        real, intent(in)  :: mean_values(:), error_term(:)
        real, intent(in), optional :: normal_random_values(:)
        real, intent(in), optional :: exceedence_probability(:)
        real, intent(in), optional :: uniform_random_values(:)
        real, intent(in), optional :: threshold

    end subroutine sample_distribution
end interface

end module sampling_mod
