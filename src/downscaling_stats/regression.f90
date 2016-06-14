module  regression_mod
    use nr
contains
    
    function compute_regression(x, training_x, training_y, coefficients) result(y)
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
        
    end function compute_regression
    
    subroutine lapack_least_squares(X, Y, B)
        implicit none
        real, intent(inout), dimension(:,:) :: X
        real, intent(inout), dimension(:)   :: Y
        real(8), intent(inout), dimension(:)   :: B
        
        real, dimension(1000) :: WORK
        integer :: n, m, nrhs, LDX, LDY, LWORK, INFO
        
        m = size(X,1)
        LDX = m
        
        n = size(X,2)
        LDY = m
        
        nrhs = 1
        
        ! First find the optimal work size (LWORK)
        LWORK = -1
        CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, WORK, LWORK, INFO )
        LWORK = MIN( 1000, INT( WORK( 1 ) ) )
        
        ! Now solve the equations X*B = Y
        CALL SGELS( 'N', M, N, NRHS, X, LDX, Y, LDY, WORK, LWORK, INFO )
        
        ! print*, "=============================="
        ! print*, X
        
        print*, "=============================="
        print*, Y(1:n)
        B(1:n) = Y(1:n)

    end subroutine lapack_least_squares


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
    

    
    subroutine logistic_regression(X, Y, TX, B)
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

    end subroutine logistic_regression


end module regression_mod
