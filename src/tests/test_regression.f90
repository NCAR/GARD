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
    logical :: pass_fail, verbose

    allocate(X(m, n))
    allocate(TX(n,m))
    allocate(Y(m))
    allocate(B(n))
    allocate(B_lapack(n))
    
    verbose = read_verbosity()
    
    do i=1,m
        X(i,1) = 1
        X(i,2) = i
        call random_number(random)
        Y(i) = i + (random-0.5)/2
    enddo
    
    TX = transpose(X)
    
    call least_squares(X, Y, TX, B)
    
    print*, "--------------------------"
    print*, "  Testing Regression "
    if (verbose) then
        print*, "--------------------------"
        print*, "      X             Y      "
        print*, "--------------------------"
        do i=1,m
            print*, X(i,2), Y(i)
        enddo
        print*, "--------------------------"
        print*, " NR based Coefficients"
        print*, B
    endif

    Y = Y*2 + 5
    call least_squares(X, Y, TX, B)
    if (verbose) then
        print*, "--------------------------"
        print*, "--------------------------"
        print*, " Y = Y * 2 + 5"
        print*, "--------------------------"
        print*, " NR based Coefficients new Y values"
        print*, "--------------------------"
        print*, B
    endif
    
    call lapack_least_squares(X, Y, B_lapack)
    if (verbose) then
        print*, "--------------------------"
        print*, " LAPACK based regression "
        print*, "--------------------------"
        print*, B_lapack
        print*, "--------------------------"
    endif
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
    
    do i=1,m
        X(i,1) = 1
        X(i,2) = i
        call random_number(random)
        Y(i) = nint(random/2 + i/real(m*1.5) )
    enddo
    ! The random numbers above are a better way to create a test for automated large sets, 
    ! but it is difficult to get it to look like a clean test for small numbers, so Y is created manually here
    if (m==10) Y = [0,0,1,0,1,1,0,1,1,1]
    print*, "--------------------------"
    print*, "Testing logistic regression"
    print*, "--------------------------"
    if (verbose) then
        print*, "      X             Y      "
        print*, "--------------------------"
        do i=1,m
            print*, X(i,2), Y(i)
        enddo
        print*, "--------------------------"
        call logistic_regression(X, Y, B)
        print*, " Logistic regression coefficients"
        print*, "--------------------------"
        print*, B
        print*, "--------------------------"
        print*, "logistic regression results"
        print*, "--------------------------"
        print*, "      X             Y          Probability"
        TX(1,:) = 1.0 / (1.0 + exp(-matmul(X, B)))
        do i=1,m
            print*, X(i,2), Y(i), TX(1,i)
        enddo
    endif
    
    if (m==10) Y = [0,0,0,0,0,1,1,1,1,1]
    call logistic_regression(X, Y, B)
    TX(1,:) = 1.0 / (1.0 + exp(-matmul(X, B)))
    if (verbose) then
        print*, "--------------------------"
        print*, " Another set of coefficients"
        print*, "--------------------------"
        print*, B
        print*, "--------------------------"
        print*, "      X             Y          Probability"
        do i=1,m
            print*, X(i,2), Y(i), TX(1,i)
        enddo
    endif
    if ((TX(1,5)<0.5).and.(TX(1,6)>0.5)) then
        print*, " PASSED"
    else
        print*, " FAILED!!!!!!!!"
    endif

contains
    function read_verbosity() result(verbosity)
        implicit none
        
        integer :: i, count, j
        character(len=1024) :: arg
        logical :: verbosity
        
        verbosity = .False.
        count = command_argument_count()
        if (count>0) then
            do j=1,count
                call get_command_argument(j, arg)
                
                if (len_trim(arg)>0) then
                    do i=1,len_trim(arg)
                        if (i<=len("verbose")) then
                            if (arg(i:i)=="verbose"(i:i)) then
                                verbosity = .True.
                            else
                                verbosity = .False.
                            endif
                        else
                            verbosity = .False.
                        endif
                    end do
                endif
                if (verbosity) exit
            end do
        endif
    end function read_verbosity
end program test_regression
