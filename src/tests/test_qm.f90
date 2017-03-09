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
program test_qm
    use quantile_mapping
    use data_structures
    use basic_stats_mod

    implicit none
    logical :: verbose
    integer*8, parameter :: n = 10000
    real, parameter :: acceptable_error = 1 ! 1e-3 <- use this with more QM segments and train on the complete dataset

    real, allocatable, dimension(:) :: input_data, matching_data

    allocate(input_data(n))
    allocate(matching_data(n))

    verbose = read_verbosity()

    !-------------------------------
    ! SETUP input data
    !-------------------------------
    call random_number(input_data)

    ! first test a simple gain and offset
    call random_number(matching_data)
    matching_data = matching_data * 10 + 3
    call test_mapping(input_data, matching_data, "Gain + Offset")

    ! test a negative gain and offset
    call random_number(matching_data)
    matching_data = matching_data * -3 - 3
    call test_mapping(input_data, matching_data, "Negative Values")

    ! Then test a fourth root transformation
    call random_number(matching_data)
    matching_data = sqrt(sqrt(matching_data)) * 2 + 3.14159
    call test_mapping(input_data, matching_data, "Double sqrt")

    ! Finally test a twentieth power transformation
    call random_number(matching_data)
    matching_data = (matching_data+1) ** 20
    call test_mapping(input_data, matching_data, "To The Twentieth")

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

    subroutine test_mapping(input_data, matching_data, name)
        implicit none
        real, intent(in), dimension(:) :: input_data, matching_data
        character(len=*), intent(in) :: name

        real, allocatable, dimension(:) :: output_data
        integer :: i, start_time, end_time, COUNT_RATE, COUNT_MAX
        integer*8 :: n

        type(qm_correction_type) :: qm

        n = size(matching_data)
        allocate(output_data(n))

        print*, ""
        print*, ""
        print*, "---------------------------"
        print*, trim(name)
        if (verbose) then
            print*, "Input:",sum(input_data)/size(input_data)
            print*, "   min/max",minval(input_data),maxval(input_data)
            print*, "To match:",sum(matching_data)/n
            print*, "   min/max",minval(matching_data),maxval(matching_data)
            print*, "---------------------------"
            if (n < 20) then
                print*, input_data
                print*, "---------------------------"
                print*, matching_data
                print*, "---------------------------"
            endif
        endif

        call system_clock(start_time)
        !-------------------------------
        ! Develop the Quantile mapping
        !-------------------------------
        call develop_qm(input_data(1:size(input_data)/2), matching_data, qm, min(300,min(size(input_data)/2, size(matching_data)/2)))

        call system_clock(end_time, COUNT_RATE, COUNT_MAX)
        if (start_time>end_time) end_time=end_time+COUNT_MAX
        if (verbose) then
            print*, "Develop_qm timing", (end_time-start_time) / real(COUNT_RATE)
        endif

        !-------------------------------
        ! Now test the quantile mapping
        !-------------------------------
        call system_clock(start_time)
        call apply_qm(input_data, output_data, qm)

        call system_clock(end_time, COUNT_RATE, COUNT_MAX)
        if (start_time>end_time) end_time=end_time+COUNT_MAX

        call show_results(output_data, matching_data, name, (end_time-start_time) / real(COUNT_RATE))
    end subroutine test_mapping

    subroutine pass_fail(error, name)
        implicit none
        real, intent(in) :: error
        character(len=*), intent(in) :: name

        if (abs(error) < acceptable_error) then
            print*, trim(name), " PASSED"
        else
            print*, trim(name), " FAILED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        endif

    end subroutine pass_fail

    !>--------------------------
    !! Print out results from a given sort test
    !!
    !!--------------------------
    subroutine show_results(qmdata, original, name, time)
        implicit none
        real, dimension(:), intent(in) :: qmdata, original
        character(len=*),   intent(in) :: name
        real,               intent(in) :: time
        real :: qm_data_stat, original_stat, error
        integer :: i, n, err

        n = size(qmdata)

        qm_data_stat = sum(qmdata)/n
        original_stat = sum(original)/size(original)
        error = 100.0 * (original_stat - qm_data_stat) / original_stat

        if (verbose) then
            print*, "---------------------------"
            print*, trim(name)
            print*, "Q-Mapped Mean = ",qm_data_stat, "Original Mean = ",original_stat
            print*, "   Error % = ", error
        endif
        call pass_fail(error, "  Mean:")

        qm_data_stat = stddev(qmdata, mean_in=qm_data_stat)
        original_stat = stddev(original, mean_in=original_stat)
        error = 100.0 * (original_stat - qm_data_stat) / original_stat

        if (verbose) then
            print*, "Q-Mapped Stddev = ",qm_data_stat, "Original stddev = ",original_stat
            print*, "   Error % = ", error
            print*, "---------------------------"
            print*, "Q-Mapped Min/Max = ",minval(qmdata),   maxval(qmdata)
            print*, "Original Min/Max = ",minval(original), maxval(original)
            print*, "---------------------------"
            if (n < 20) then
                print*, qmdata
                print*, "---------------------------"
            endif
            print*, "   Elapsed time = ", time
        endif
        call pass_fail(error, "  Stddev:")

    end subroutine show_results

end program test_qm
