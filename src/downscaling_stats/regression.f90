module  regression_mod
    use nr
contains
    
    function nr_compute_regression(x, training_x, training_y, coefficients) result(y)
        implicit none
        real,    intent(in),    dimension(:)   :: x
        real,    intent(in),    dimension(:,:) :: training_x
        real,    intent(in),    dimension(:)   :: training_y
        real(8), intent(inout), dimension(:)   :: coefficients
        
        real :: y
        
        integer :: nvars, ntimes, i
        real, dimension(:,:), allocatable :: transpose_x
        
        nvars  = size(x)
        ntimes = size(training_x,1)
        
        allocate( transpose_x(nvars,ntimes) )
        
        if (maxval(training_y) == 0) then
            y=0
            return
        endif
        
        if (any(maxval(abs(training_x),dim=1) == 0.0)) then
            y=0
            return
        endif

        do i=1,nvars
            transpose_x(i,:) = training_x(:,i)
        enddo
        
        call least_squares(training_x, training_y, transpose_x, coefficients)
        
        y = 0
        if (maxval(abs(coefficients)) == 0) then
            y = sum(training_y)/size(training_y)
        endif
        do i=1,nvars
            y = y + coefficients(i) * x(i)
        end do
        
    end function nr_compute_regression
    
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
            write(*,*) "ERROR in SGELS with argument:", 0-INFO
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


    !
    ! Solve linear equation for x (Ax = b => x = bA^-1) using LU decomposition and back substitution.
    ! Input:
    !   X  = An m by n array.
    !   TX = Precalculated transpose array of X, size n by m
    !   Y  = An m-element vector containing the right-hand side of the linear system Ax = b. 
    ! Output: 
    !   B  = An n-element vector.
    !
    subroutine least_squares(X, Y, TX, B)
        implicit none

        real,    intent(in)     :: X(:,:)
        real,    intent(in)     :: Y(:)
        real,    intent(in)     :: TX(:,:)
        real(8), intent(inout)  :: B(:)

        real(8), allocatable :: A(:,:)
        integer, allocatable :: indx(:)
        integer :: nvars, ntimes
        real(8) :: d

        nvars = size(X,2) -1
        ntimes = size(Y)

        ! if (.not.allocated(B)) then
        !     allocate(B(nvars+1))
        ! endif
        allocate(A(nvars+1,nvars+1))
        allocate(indx(nvars+1))

        
        B = matmul(TX, Y)
        A = matmul(TX, X)

        if (any(maxval(abs(A),dim=2) == 0.0)) then
            !$omp critical
            print*, "ERROR"
            !$omp end critical
            B = 0
            return
        endif
        call ludcmp(A, indx, d )

        if (ANY(ABS(A) < 9.99999968E-15)) then
            B(:) = 0.0d0
            !$omp critical
            print*, "Warning, LUdcmp produced a zero."
            !$omp end critical
            return
        endif

        call lubksb(A, indx, B)
        
        deallocate(A)
        deallocate(indx)

    end subroutine least_squares
    

    
    subroutine nr_logistic_regression(X, Y, TX, B)
        implicit none
        real, intent(in) :: X(:,:)
        real, intent(in) :: Y(:)
        real, intent(in) :: TX(:,:)
        real(8), allocatable, intent(out) :: B(:)

        real(8), allocatable :: Yp(:), P(:), BN(:), V(:,:)
        real, allocatable :: YN(:), XV(:,:)
        integer :: nvars, ntimes, i, t, f, it
        real :: d

        nvars = size(X,2) -1
        ntimes = size(Y)

        allocate(B(nvars+1))
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
            P = 1.0d0 / (1.0d0 + exp(-matmul(X, B))) 
            if (ANY(P > 0.97)) then
                !print *, "WARNING: logistic regression diverging"
                f = 1
            else

                YN = Yp - P
                V = 0.0d0
                do t = 1, ntimes, 1
                    V(t,t) = P(t)*(1.0d0-P(t))
                enddo
                XV = matmul(V,X)
                call least_squares(XV, YN, TX, BN)

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
                deallocate(BN)

            endif
            it = it+1
        end do

    end subroutine nr_logistic_regression

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
            P = 1.0 / (1.0 + exp(-matmul(X, B))) 
            if (ANY(P > 0.97)) then
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
