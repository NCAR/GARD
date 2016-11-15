!>----------------------------------------
!! Simple program to test and time sorting routines
!!
!! Generates n random numbers, uses heap_sort to sort them,
!!  Then tests that the sort succeeded and prints timing info
!!
!! Next this is repeated for quick_sort.
!!
!! Initial data show that this implementation of quick_sort is ~2x faster than
!!  the (somewhat naive) implementation of heap_sort
!!
!!----------------------------------------
program test_sort
    use sort_mod

    implicit none

    real, parameter :: valid_error = 1e-5
    integer*8, parameter :: n = 1000000

    real, allocatable, dimension(:) :: to_sort, sorted
    integer :: i, start_time, end_time, COUNT_RATE, COUNT_MAX

    allocate(to_sort(n))
    allocate(sorted(n))

    !-------------------------------
    ! SETUP sorting
    !-------------------------------
    call random_number(to_sort)
    sorted=0
    ! test heap_sort (system_clock is used to time how long it takes)

    call system_clock(start_time)
    !-------------------------------
    ! HEAP SORT
    !-------------------------------
    call heap_sort(to_sort, sorted)

    call system_clock(end_time, COUNT_RATE, COUNT_MAX)
    if (start_time>end_time) end_time=end_time+COUNT_MAX

    call show_results(sorted, "HEAP_SORT", (end_time-start_time) / real(COUNT_RATE))


    !-------------------------------
    ! SETUP sorting
    !-------------------------------
    ! now test quick_sort
    call random_number(to_sort)
    sorted = 0

    call system_clock(start_time)
    !-------------------------------
    ! QUICK SORT
    !-------------------------------
    call quick_sort(to_sort, sorted)

    call system_clock(end_time, COUNT_RATE, COUNT_MAX)
    if (start_time>end_time) end_time=end_time+COUNT_MAX

    call show_results(sorted, "QUICK_SORT", (end_time-start_time) / real(COUNT_RATE))

    !-------------------------------
    ! SETUP sorting
    !-------------------------------
    ! now test quick_sort
    call random_number(to_sort)

    call system_clock(start_time)
    !-------------------------------
    ! QUICK SORT in place
    !-------------------------------
    call sort(to_sort)

    call system_clock(end_time, COUNT_RATE, COUNT_MAX)
    if (start_time>end_time) end_time=end_time+COUNT_MAX

    call show_results(to_sort, "QUICK_SORT (in place)", (end_time-start_time) / real(COUNT_RATE))

contains
    !>--------------------------
    !! Print out results from a given sort test
    !!
    !!--------------------------
    subroutine show_results(data, name, time)
        implicit none
        real, dimension(:), intent(in) :: data
        character(len=*),   intent(in) :: name
        real,               intent(in) :: time

        integer :: i, n, err

        n = size(data)

        print*, "---------------------------"
        err = 0
        do i=1,n-1
            if (sorted(i+1)<sorted(i)) then
                err = 1
            endif
        end do
        print*, trim(name)
        if ( err < valid_error ) then
            print*, "  PASSED"
        else
            print*, "  FAILED:", err
        endif
        print*, "---------------------------"
        print*, "   Elapsed time = ", time

    end subroutine show_results

end program test_sort
