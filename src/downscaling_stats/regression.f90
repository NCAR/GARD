module  regression_mod
contains
    
    function compute_regression(x, training_x, training_y, coefficients, error, weights) result(y)
        implicit none
        real,    intent(in),    dimension(:)   :: x
        real,    intent(in),    dimension(:,:) :: training_x
        real,    intent(in),    dimension(:)   :: training_y
        real(8), intent(inout), dimension(:)   :: coefficients
        real,    intent(inout),                             optional :: error
        real,    intent(inout), dimension(:),  allocatable, optional :: weights
        
        real,    allocatable,   dimension(:,:) :: training_x_lp
        real,    allocatable,   dimension(:)   :: training_y_lp
        
        real :: y
        
        integer :: nvars, ntimes, i, j, info, useable_vars, thisvar, nonzero_count
        integer, dimension(:), allocatable :: varlist
        real,    dimension(:), allocatable :: useful_x
        
        if (maxval(training_y) == minval(training_y)) then
            y = training_y(1)
            coefficients = 0
            if (present(error)) error=0
            return
        endif
        
        nvars  = size(x)
        ntimes = size(training_x,1)
        
        allocate(varlist(nvars))
        varlist = -1
        useable_vars = 0
        
        ! find all non-zero vars
        do i=1,nvars
            nonzero_count = 0
            do j=1,ntimes
                if (abs(training_x(j,i))>0) then
                    nonzero_count = nonzero_count + 1
                endif
            enddo
            if ( nonzero_count > nvars ) then
                useable_vars = useable_vars + 1
                varlist(i) = i
            end if
        enddo

        ! if there are no useable training variables other than the constant, just return (this should never happen?)
        if (useable_vars <= 1) then
            y = 1e20 ! error value so it will recompute from analogs if possible
            coefficients = 0
            if (present(error)) error=0
            return
        endif
        
        allocate(useful_x(useable_vars))
        allocate(training_x_lp(ntimes, useable_vars))
        
        thisvar=1
        do i=1,nvars
            if (varlist(i) /= -1) then
                training_x_lp(:,thisvar) = training_x(:,i)
                useful_x(thisvar) = x(i)
                thisvar = thisvar + 1
            endif
        enddo
        
        allocate(training_y_lp(ntimes))
        training_y_lp = training_y
        
        
        call lapack_least_squares(training_x_lp, training_y_lp, coefficients, info=info)
        if (info/=0) then
            y = 1e20 ! error value so it will recompute from analogs
            ! !$omp critical (print_lock)
            ! print*, "Regression error, using value:",y
            ! !$omp end critical (print_lock)
            coefficients = 0
            if (present(error)) error=0
            return
        endif
        
        if (maxval(abs(coefficients)) == 0) then
            y = 1e20
            if (present(error)) error=0
            return
        endif
        
        y = 0
        do i=1,useable_vars
            y = y + coefficients(i) * useful_x(i)
        end do
        
        if (present(error)) then
            training_y_lp = 0
            thisvar       = 1
            do i=1, nvars
                if (varlist(i) /= -1) then
                    training_y_lp = training_y_lp + coefficients(thisvar) * training_x(:,i)
                    thisvar       = thisvar + 1
                endif
            enddo
            ! root mean square error
            ! this is the regression error
            
            ! optionally use the weights computed in analog_weights to weight the error calculation
            ! weights should sum to 1, but just in case, we divide by sum(weights)
            if (present(weights)) then
                if (allocated(weights)) then
                    if (size(weights)/=size(training_y_lp)) then
                        write(*,*) "ERROR size of weights /= data"
                        write(*,*), shape(weights), shape(training_y_lp)
                    endif
                    error = sqrt( sum(((training_y_lp - training_y)**2)*weights) / sum(weights) )
                else
                    error = sqrt( sum((training_y_lp - training_y)**2) / ntimes )
                endif
            else
                error = sqrt( sum((training_y_lp - training_y)**2) / ntimes )
            endif
        endif
        
    end function compute_regression

    !
    ! Solve linear equation for x (Ax = b => x = bA^-1) using the lapack SGELS library call
    ! Input:
    !   X  = An m by n array.
    !   Y  = An m-element vector containing the right-hand side of the linear system Ax = b. 
    ! Output: 
    !   B  = An n-element vector.
    !
    subroutine lapack_least_squares(X, Y, B, working_space, info)
        implicit none
        real,    intent(inout), dimension(:,:) :: X
        real,    intent(inout), dimension(:)   :: Y
        real(8), intent(inout), dimension(:)   :: B
        real,    intent(inout), dimension(:), optional :: working_space ! temporary working space
        integer, intent(inout),               optional :: info
        
        ! it could be more efficient to keep one space allocated and pass it in every time (working_space), 
        ! in practice, having a modest "WORK" allocated on the stack seems to be at least as fast. 
        real, dimension(1000) :: WORK
        integer :: n, m, nrhs, LDX, LDY, LWORK, innerINFO
        
        m = size(X,1)
        LDX = m
        
        n = size(X,2)
        LDY = m
        
        nrhs = 1

        if (present(working_space)) then
            ! First find the optimal work size (LWORK)
            LWORK = -1
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, working_space, LWORK, innerINFO )
            LWORK = MIN( size(working_space), INT( working_space( 1 ) ) )
            
            ! Now solve the equations X*B = Y
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, working_space, LWORK, innerINFO )
        else
            ! First find the optimal work size (LWORK)
            LWORK = -1
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, WORK, LWORK, innerINFO )
            LWORK = MIN( 1000, INT( WORK( 1 ) ) )
            
            ! Now solve the equations X*B = Y
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, WORK, LWORK, innerINFO )
        endif
        if (innerINFO/=0) then
            ! !$omp critical (print_lock)
            ! if (innerINFO<0) then
            !     write(*,*) "ERROR in SGELS with argument:", 0-innerINFO
            ! else
            !     write(*,*) "ERROR in SGELS A does not have full rank, position:",innerINFO
            ! endif
            ! !$omp end critical (print_lock)
            Y(1) = sum(Y)/size(Y)
            Y(2:)= 0
        endif
        
        if (present(info)) info = innerINFO
        B(1:n) = Y(1:n)

    end subroutine lapack_least_squares

    subroutine weighted_least_squares(X, Y, B, W)
        implicit none
        real, intent(inout), dimension(:,:) :: X
        real, intent(inout), dimension(:)   :: Y, W
        real(8), intent(inout), dimension(:)   :: B
        real, dimension(:,:), allocatable :: X_weighted, W_full
        real, dimension(:),   allocatable :: Y_weighted, Bs, output_y
        
        real, dimension(10000) :: WORK
        integer :: i, n, m, nrhs, LDX, LDY, LDW, LWORK, INFO, p
        
        ! ntimes
        m = size(X,1)
        LDX = m
        LDY = m
        LDW = m
        p   = m
        
        ! nvars
        n = size(X,2)
        
        nrhs = 1

        allocate(output_y(p))
        allocate(Bs(n))
        allocate(W_full(m,m))
        W_full = 0
        do i=1,m
            W_full(i,i) = 1/W(i)
        enddo
        
        ! print*, "entering sggglm"
        LWORK = 10000
        call sggglm(m, n, p, X, LDX, W_full, LDW, Y, Bs, output_y, WORK, LWORK, INFO)
        
        B(1:n) = Bs(1:n)

    end subroutine weighted_least_squares


    
    function compute_logistic_regression(x, training_x, training_y, coefficients, threshold) result(output)
        implicit none
        real,    intent(in) :: x(:)
        real,    intent(in) :: training_x(:,:)
        real,    intent(in) :: training_y(:)
        real(8), intent(out):: coefficients(:)
        real,    intent(in) :: threshold
        real :: output
        
        real, dimension(:), allocatable :: binary_values
        integer :: n
        
        n = size(training_y)
        allocate(binary_values(n))
        binary_values = 0
        where(training_y > threshold) binary_values = 1
        
        call logistic_regression(training_x, binary_values, coefficients)
        
        output = dot_product(x, coefficients)
        ! at this point exp(-output) is foolishly large, much higher will cause an overflow
        output = max(-80.0, output)
        
        output = 1.0 / (1.0 + exp(-output))
        
    end function compute_logistic_regression

    subroutine logistic_regression(X, Y, B)
        implicit none
        real,    intent(in) :: X(:,:)
        real,    intent(in) :: Y(:)
        real(8), intent(out) :: B(:)

        real,    allocatable :: Yp(:), P(:), V(:,:)
        real(8), allocatable :: BN(:)
        real,    allocatable :: YN(:), XV(:,:)
        
        integer :: nvars, ntimes, i, t, f, it
        real :: d
        integer :: info
        
        if (all(Y > 0)) then
            B(1) = 88
            B(2:)= 0
            return
        elseif (all(Y==0)) then
            B(1) = -88
            B(2:)= 0
            return
        endif

        nvars = size(X,2) -1
        ntimes = size(Y)

        allocate(BN(nvars+1))
        allocate(Yp(ntimes))
        allocate(YN(ntimes))
        allocate(P(ntimes))
        allocate(V(ntimes,ntimes))
        allocate(XV(ntimes, nvars+1))

        do t = 1, ntimes, 1
            if(Y(t) .ne. 0.0) then
                Yp(t) = 1.0d0
            else
                Yp(t) = 0.0d0
            end if
        end do

        B = 0.0d0
        f = 0
        i = 0
        it = 0
        do while (f /= 1)
            
            P = max(-80.0, matmul(X, B))
            P = 1.0 / (1.0 + exp(-P))
            if (ANY(P > 0.9999999)) then
                ! !$omp critical (print_lock)
                !print*, "WARNING: logistic regression diverging"
                ! !$omp end critical (print_lock)
                f = 1
            else

                YN = Yp - P
                V = 0.0
                do t = 1, ntimes, 1
                    V(t,t) = P(t)*(1.0-P(t))
                enddo
                XV = matmul(V,X)
                call lapack_least_squares(XV, YN, BN, info=info)
                if (info /= 0) then
                    B(1) = -1 * log( (1.0 / (sum(Y)/size(Y))) - 1.0)
                    B(2:)= 0
                    ! !$omp critical (print_lock)
                    ! print*, B(1)
                    ! !$omp end critical (print_lock)
                    return
                endif
                    

                f = 1
                do i = 1, nvars+1, 1
                    if(BN(i) .GT. 1.0E-04 .OR. BN(i) .LT. -1.0E-04) then
                        f = 0
                    end if
                end do
                if(it > 8) then
                    ! !$omp critical (print_lock)
                    !print*, "WARNING: logistic regression failed to converge"
                    ! !$omp end critical (print_lock)
                    f = 1
                endif

                B = B + BN

            endif
            it = it+1
        end do

    end subroutine logistic_regression


end module regression_mod
