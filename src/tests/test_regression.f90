program test_regression
    
    use regression_mod
    
    real,    dimension(:,:), allocatable :: X, TX
    real,    dimension(:),   allocatable :: Y
    real(8), dimension(:),   allocatable :: B, B_lapack
    
    integer :: random_size
    real, parameter :: MAX_ERROR = 1e-4
    integer, parameter :: n=2, m=10
    real :: random
    integer :: i
    logical :: pass_fail

    allocate(X(m, n))
    allocate(TX(n,m))
    allocate(Y(m))
    allocate(B(n))
    allocate(B_lapack(n))
    
    do i=1,m
        X(i,1) = 1
        X(i,2) = i
        call random_number(random)
        Y(i) = i + (random-0.5)/2
    enddo
    
    TX = transpose(X)
    
    call least_squares(X, Y, TX, B)
    
    print*, "--------------------------"
    print*, "      X             Y      "
    print*, "--------------------------"
    do i=1,m
        print*, X(i,2), Y(i)
    enddo
    print*, "--------------------------"
    print*, " NR based Coefficients"
    print*, B

    Y = Y*2 + 5
    call least_squares(X, Y, TX, B)
    print*, "--------------------------"
    print*, "--------------------------"
    print*, " Y = Y * 2 + 5"
    print*, "--------------------------"
    print*, " NR based Coefficients new Y values"
    print*, "--------------------------"
    print*, B
    
    call lapack_least_squares(X, Y, B_lapack)
    print*, "--------------------------"
    print*, " LAPACK based regression "
    print*, "--------------------------"
    print*, B_lapack
    
    print*, "--------------------------"
    pass_fail = .True.
    do i=1,n
        if (abs(B_lapack(i)-B(i)) > MAX_ERROR) then
            pass_fail = .False.
        endif
    enddo
    print*, "--------------------------"
    if (pass_fail) then
        print*, " PASSED"
    else
        print*, " FAILED!!!!!!!!!!!"
    endif
    print*, "--------------------------"
    
end program test_regression
