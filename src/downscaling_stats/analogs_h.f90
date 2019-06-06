module analog_mod

    implicit none

    integer, parameter :: MIN_NUMBER_ANALOGS = 20

interface
    module subroutine find_analogs(analogs, match, input, n, threshold, weights, skip_analog, stochastic_analog_perturbation)
        implicit none
        integer, intent(inout), dimension(:), allocatable  :: analogs
        real,    intent(in),    dimension(:)   :: match
        real,    intent(in),    dimension(:,:) :: input
        integer, intent(in)                    :: n
        real,    intent(in)                    :: threshold
        real,    intent(inout), dimension(:), allocatable, optional :: weights
        integer, intent(in),    optional       :: skip_analog ! specify an element of the input data to skip
        real,    intent(in),    optional       :: stochastic_analog_perturbation

    end subroutine find_analogs

    ! compute the mean of input[analogs]
    module function compute_analog_mean(input, analogs, mask, weights) result(mean)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in), dimension(:), optional :: mask
        real,    intent(in), dimension(:), optional :: weights
        real                              :: mean

    end function compute_analog_mean

    ! compute the root mean square error between input[analogs] and y_hat
    module function compute_analog_error(input, analogs, y_hat, mask, weights) result(error)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in)               :: y_hat
        real,    intent(in), dimension(:), optional :: mask
        real,    intent(in), dimension(:), optional :: weights
        real                              :: error

    end function compute_analog_error


    module function compute_analog_exceedance(input, analogs, weights) result(probability)
        implicit none
        real,    intent(in), dimension(:) :: input ! this is a mask 1 = exceeded, 0 = not exceeded
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in), dimension(:), optional :: weights
        real                              :: probability

    end function compute_analog_exceedance

end interface
end module analog_mod
