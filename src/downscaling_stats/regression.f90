module  regression_mod
    use nr
contains
    
    function compute_regression(x, training_x, training_y) result(y)
        implicit none
        real, intent(in), dimension(:)   :: x
        real, intent(in), dimension(:,:) :: training_x
        real, intent(in), dimension(:)   :: training_y
        
        real :: y
        
        integer :: nvars, ntimes, i
        real(8), dimension(:),allocatable :: coefficients
        real, dimension(:,:), allocatable :: transpose_x
        
        nvars  = size(x)
        ntimes = size(training_x,2)
        
        allocate( transpose_x(ntimes,nvars) )
        
        do i=1,nvars
            transpose_x(:,i) = training_x(i,:)
        enddo
        
        call least_squares(training_x, training_y, transpose_x, coefficients)
        
        y = 0
        do i=1,nvars
            y = y + coefficients(i) * x(i)
        end do
        
    end function compute_regression

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

        real, intent(in) :: X(:,:)
        real, intent(in) :: Y(:)
        real, intent(in) :: TX(:,:)
        real(8), allocatable, intent(out) :: B(:)

        real(8), allocatable :: A(:,:)
        integer, allocatable :: indx(:)
        integer :: nvars, ntimes
        real(8) :: d

        nvars = size(X,2) -1
        ntimes = size(Y)

        if (.not.allocated(B)) then
            allocate(B(nvars+1))
        endif
        allocate(A(nvars+1,nvars+1))
        allocate(indx(nvars+1))

        
        B = matmul(TX, Y)
        A = matmul(TX, X)

        call ludcmp(A, indx, d )

        if (ANY(ABS(A) < 9.99999968E-15)) then
            B(:) = 0.0d0
            print *, "Warning, LUdcmp produced a zero."
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
