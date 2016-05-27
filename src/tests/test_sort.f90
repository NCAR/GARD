program test_sort
    use sort_mod
    
    implicit none

    real, parameter :: valid_error = 1e-5
    integer, parameter :: n = 10
    
    real, allocatable, dimension(:) :: to_sort, sorted, correct
    integer :: i
    real :: maxerr
    
    allocate(to_sort(n))
    allocate(sorted(n))
    allocate(correct(n))
    
    to_sort = [3.14159, 4., 2., 2., 40., 9., 1., 1.1, 9., -1.0]
    correct = [-1.0, 1.0, 1.1, 2.0, 2.0, 3.14159, 4.0, 9.0, 9.0, 40.]
    
    print*, "---------------------------"
    do i=1,n
        print*, to_sort(i)
    end do
    print*, "---------------------------"
    
    to_sort = [3.14159, 4., 2., 2., 40., 9., 1., 1.1, 9., -1.0]
    call sort(to_sort, sorted)
    
    print*, "---------------------------"
    maxerr = 0
    do i=1,n
        maxerr = max(maxerr, abs(correct(i) - sorted(i)))
        ! print*, correct(i), sorted(i)
    end do
    print*, "HEAPSORT: "
    if ( maxerr < valid_error ) then
        print*, "  PASSED"
    else
        print*, "  FAILED"
    endif
    print*, "---------------------------"

    to_sort = [3.14159, 4., 2., 2., 40., 9., 1., 1.1, 9., -1.0]
    sorted = 0
    call quick_sort(to_sort, sorted)
    
    print*, "---------------------------"
    maxerr = 0
    do i=1,n
        maxerr = max(maxerr, abs(correct(i) - sorted(i)))
        print*, correct(i), sorted(i)
    end do
    print*, "Quicksort: "
    if ( maxerr < valid_error ) then
        print*, "  PASSED"
    else
        print*, "  FAILED"
    endif
    print*, "---------------------------"


end program test_sort
