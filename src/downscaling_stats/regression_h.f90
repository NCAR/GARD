!> -------------------------------------------
!! Supply regression utilities to GARD
!!
!! Provides least squares regression and logistic regression
!! includes output of both the coefficients and the residual error, optionally computed with supplied weights
!!
!! --------------------------------------------
module  regression_mod
interface
    module function compute_regression(x, training_x, training_y, coefficients, y_test, error, weights, used_vars) result(y)
        implicit none
        real,    intent(in),    dimension(:)   :: x
        real,    intent(in),    dimension(:,:) :: training_x
        real,    intent(in),    dimension(:)   :: training_y
        real(8), intent(inout), dimension(:)   :: coefficients
        real,    intent(inout),                             optional :: error
        real,    intent(in),    dimension(:),               optional :: y_test
        real,    intent(inout), dimension(:),  allocatable, optional :: weights
        integer, intent(inout), dimension(:),               optional :: used_vars
        real :: y
    end function compute_regression

    module function compute_logistic_regression(x, training_x, training_y, coefficients, threshold) result(output)
        implicit none
        real,    intent(in) :: x(:)
        real,    intent(in) :: training_x(:,:)
        real,    intent(in) :: training_y(:)
        real(8), intent(out):: coefficients(:)
        real,    intent(in) :: threshold
        real :: output
    end function compute_logistic_regression

end interface

end module regression_mod
