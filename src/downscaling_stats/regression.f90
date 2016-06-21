module  regression_mod
contains
    
    function compute_regression(x, training_x, training_y, coefficients, error) result(y)
        implicit none
        real,    intent(in),    dimension(:)   :: x
        real,    intent(in),    dimension(:,:) :: training_x
        real,    intent(in),    dimension(:)   :: training_y
        real(8), intent(inout), dimension(:)   :: coefficients
        real,    intent(inout), optional       :: error
        real,    allocatable,   dimension(:,:) :: training_x_lp
        real,    allocatable,   dimension(:)   :: training_y_lp
        
        real :: y
        
        integer :: nvars, ntimes, i
        
        nvars  = size(x)
        ntimes = size(training_x,1)
        
        allocate(training_x_lp(ntimes, nvars))
        training_x_lp = training_x
        
        allocate(training_y_lp(ntimes))
        training_y_lp = training_y
        
        if (maxval(training_y) == 0) then
            y=0
            return
        endif
        
        if (any(maxval(abs(training_x),dim=1) == 0.0)) then
            y=0
            return
        endif
        
        call lapack_least_squares(training_x_lp, training_y_lp, coefficients)
        
        y = 0
        if (maxval(abs(coefficients)) == 0) then
            y = sum(training_y) / size(training_y)
        endif
        do i=1,nvars
            y = y + coefficients(i) * x(i)
        end do
        
        if (present(error)) then
            training_y_lp = 0
            do i=1,nvars
                training_y_lp = training_y_lp + coefficients(i) * training_x(:,i)
            enddo
            ! root mean square error
            error = sqrt( sum((training_y_lp - training_y)**2) / ntimes )
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
    subroutine lapack_least_squares(X, Y, B, working_space)
        implicit none
        real, intent(inout), dimension(:,:) :: X
        real, intent(inout), dimension(:)   :: Y
        real(8), intent(inout), dimension(:)   :: B
        real, intent(inout), dimension(:), optional :: working_space ! temporary working space
        ! it could be more efficient to keep one space allocated and pass it in every time, 
        ! in practice, having a modest "WORK" allocated on the stack seems to be at least as fast. 
        real, dimension(1000) :: WORK
        integer :: n, m, nrhs, LDX, LDY, LWORK, INFO
        
        m = size(X,1)
        LDX = m
        
        n = size(X,2)
        LDY = m
        
        nrhs = 1

        if (present(working_space)) then
            ! First find the optimal work size (LWORK)
            LWORK = -1
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, working_space, LWORK, INFO )
            LWORK = MIN( size(working_space), INT( working_space( 1 ) ) )
            
            ! Now solve the equations X*B = Y
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, working_space, LWORK, INFO )
        else
            ! First find the optimal work size (LWORK)
            LWORK = -1
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, WORK, LWORK, INFO )
            LWORK = MIN( 1000, INT( WORK( 1 ) ) )
            
            ! Now solve the equations X*B = Y
            CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, WORK, LWORK, INFO )
        endif
        if (INFO/=0) then
            !$omp critical (print_lock)
            write(*,*) "ERROR in SGELS with argument:", 0-INFO
            !$omp end critical (print_lock)
            Y = 0
        endif

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
                !print *, "WARNING: logistic regression diverging"
                f = 1
            else

                YN = Yp - P
                V = 0.0
                do t = 1, ntimes, 1
                    V(t,t) = P(t)*(1.0-P(t))
                enddo
                XV = matmul(V,X)
                call lapack_least_squares(XV, YN, BN)

                f = 1
                do i = 1, nvars+1, 1
                    if(BN(i) .GT. 1.0E-04 .OR. BN(i) .LT. -1.0E-04) then
                        f = 0
                    end if
                end do
                if(it > 8) then
                    !print *, "WARNING: logistic regression failed to converge"
                    f = 1
                endif

                B = B + BN

            endif
            it = it+1
        end do

    end subroutine logistic_regression


end module regression_mod
