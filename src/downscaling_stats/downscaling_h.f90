module downscaling_mod

    use data_structures

    implicit none

interface
    module subroutine downscale(training_atm, training_obs, predictors, output, options)
            implicit none
            type(atm),    intent(inout) :: training_atm, predictors
            type(obs),    intent(inout) :: training_obs
            type(results),intent(out)   :: output
            type(config), intent(inout) :: options

    end subroutine downscale

end interface
end module downscaling_mod
