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

    integer*8, parameter :: n = 100000
    
    real, allocatable, dimension(:) :: input_data, matching_data
    
    allocate(input_data(n))
    allocate(matching_data(n))
    
    !-------------------------------
    ! SETUP input data
    !-------------------------------
    call random_number(input_data)

    ! first test a simple gain and offset
    call random_number(matching_data)
    matching_data = matching_data * 10 + 3
    call test_mapping(input_data, matching_data, "Gain + Offset")

    ! Then test a fourth root transformation
    call random_number(matching_data)
    matching_data = sqrt(sqrt(matching_data)) * 2 + 3.14159
    call test_mapping(input_data, matching_data, "Double sqrt")

    ! Finally test a twentieth power transformation
    call random_number(matching_data)
    matching_data = (matching_data+1) ** 20
    call test_mapping(input_data, matching_data, "To The Twentieth")
    
contains
    
    subroutine test_mapping(input_data, matching_data, name)
        implicit none
        real, allocatable, intent(in), dimension(:) :: input_data, matching_data
        character(len=*), intent(in) :: name
        
        real, allocatable, dimension(:) :: output_data
        integer :: i, start_time, end_time, COUNT_RATE, COUNT_MAX
        integer*8 :: n
        
        type(qm_correction_type) :: qm

        n = size(matching_data)
        allocate(output_data(n))
        
        print*, "---------------------------"
        print*, trim(name)
        print*, "Input:",sum(input_data)/size(input_data)
        print*, "To match:",sum(matching_data)/n
        print*, "---------------------------"
        
        call system_clock(start_time)
        !-------------------------------
        ! Develop the Quantile mapping
        !-------------------------------
        call develop_qm(input_data, matching_data, qm, 10000)
        
        call system_clock(end_time, COUNT_RATE, COUNT_MAX)
        if (start_time>end_time) end_time=end_time+COUNT_MAX
        print*, "Develop_qm timing", (end_time-start_time) / real(COUNT_RATE)

        !-------------------------------
        ! Now test the quantile mapping
        !-------------------------------
        call system_clock(start_time)
        !-------------------------------
        ! QUICK SORT
        !-------------------------------
        call apply_qm(input_data, output_data, qm)
        ! print*, output_data
        
        call system_clock(end_time, COUNT_RATE, COUNT_MAX)
        if (start_time>end_time) end_time=end_time+COUNT_MAX

        call show_results(output_data, matching_data, name, (end_time-start_time) / real(COUNT_RATE))
    end subroutine test_mapping
    
    !>--------------------------
    !! Print out results from a given sort test
    !!
    !!--------------------------
    subroutine show_results(qmdata, original, name, time)
        implicit none
        real, dimension(:), intent(in) :: qmdata, original
        character(len=*),   intent(in) :: name
        real,               intent(in) :: time
        
        integer :: i, n, err
        
        n = size(qmdata)
        
        print*, "---------------------------"
        print*, trim(name)
        print*, "Q-Mapped Mean = ",sum(qmdata)/n, "Original Mean = ",sum(original)/size(original)
        print*, "   Error % = ", 100.0 * (sum(original)/size(original) - sum(qmdata)/n) / (sum(original)/size(original))
        print*, "Q-Mapped Stddev = ",stddev(qmdata), "Original stddev = ",stddev(original)
        print*, "   Error % = ", 100.0 * (stddev(original) - stddev(qmdata)) / stddev(original)
        print*, "---------------------------"
        print*, "   Elapsed time = ", time
        
    end subroutine show_results

end program test_qm
