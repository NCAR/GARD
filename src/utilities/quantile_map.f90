!>------------------------------------------------
!! Perform a basic quantile mapping from one variable to another
!!
!! Develop or apply quantile mapping.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module quantile_mapping
    use sort_mod, only : sort
    use data_structures

    implicit none
    private
    public :: develop_qm
    public :: apply_qm
    public :: reverse_qm

    integer, parameter :: DEFAULT_QM_SEGMENTS = 100
    integer, parameter :: EXTREME_WINDOW      = 5
    real,    parameter :: SMALL_VALUE         = 1e-10

contains

    subroutine develop_qm(input_data, data_to_match, qm, n_segments)
        implicit none
        real, intent(in), dimension(:) :: input_data, data_to_match
        integer, intent(in), optional :: n_segments
        type(qm_correction_type), intent(out) :: qm

        integer :: ni, nm, nseg
        integer :: i, top, bottom
        real    :: input_i, input_step, match_i, match_step ! these are real to handle more bins than elements
        real    :: denominator
        real, dimension(:), allocatable :: input_data_sorted, data_to_match_sorted

        if ( present(n_segments) ) then
            nseg = n_segments
        else
            nseg = DEFAULT_QM_SEGMENTS
        end if

        ! setup the sorted input data to process
        ni = size(input_data)
        allocate(input_data_sorted(ni))
        call sort(input_data, input_data_sorted)

        ! setup the sorted matching data to process
        nm = size(data_to_match)
        allocate(data_to_match_sorted(nm))
        call sort(data_to_match, data_to_match_sorted)

        if (nseg<10) then
            write(*,*), "Quantile mapping with <10 segments not recommended, using 10 segments. "
            nseg=10
        endif
        input_step = ni/real(nseg-1)
        match_step = nm/real(nseg-1)

        input_i = 1
        match_i = 1

        if (allocated(qm%start_idx)) then
            deallocate(qm%start_idx,qm%end_idx,qm%slope,qm%offset)
        endif
        allocate(qm%start_idx(nseg))
        allocate(qm%end_idx(nseg))
        allocate(qm%slope(nseg))
        allocate(qm%offset(nseg))

        ! first handle the bottom of the QM segment
        bottom = 1
        top = max(1, min(EXTREME_WINDOW, int(input_step/4)))
        qm%start_idx(1) = sum(input_data_sorted( bottom:top )) / (top - bottom + 1)

        bottom = top
        top = 1 + input_step
        qm%end_idx(1)   = sum(input_data_sorted( bottom:top )) / (top - bottom + 1)

        bottom = 1
        top = max(1, min(EXTREME_WINDOW, int(match_step/4)) )
        qm%offset(1)    = sum(data_to_match_sorted( bottom:top )) / (top - bottom + 1)
        denominator     = qm%end_idx(1) - qm%start_idx(1)
        if (denominator < SMALL_VALUE) then
            qm%slope(1) = 0
        else
            bottom = top
            top = 1 + match_step
            qm%offset(2) = sum(data_to_match_sorted( bottom:top ))/ (top - bottom + 1)
            qm%slope(1)  = ( qm%offset(2) - qm%offset(1)) &
                           / denominator
        endif

        input_i = input_step
        match_i = match_step
        do i=2,nseg-1
            qm%start_idx(i) = qm%end_idx(i-1)
            bottom          = int(input_i)
            top             = int(input_i + input_step)
            qm%end_idx(i)   = sum(input_data_sorted( bottom:top )) / (top - bottom + 1)
            bottom          = int(match_i)
            top             = int(match_i + match_step)
            qm%offset(i+1)  = sum(data_to_match_sorted( bottom:top ))/ (top - bottom + 1)

            denominator     = qm%end_idx(i) - qm%start_idx(i)
            if (denominator < SMALL_VALUE) then
                qm%slope(i) = 0
            else
                qm%slope(i) = (qm%offset(i+1) - qm%offset(i)) &
                             / denominator
            endif

            input_i = input_i + input_step
            match_i = match_i + match_step

        end do

        ! For the last segment, just fit to an average over the top 2-5 values so we can get a more extreme fit
        ! that isn't dependant just on the absolute highest value.
        i = nseg
        qm%start_idx(i) = qm%end_idx(i-1)
        bottom          = min(max(int(input_i - input_step/4), ni - EXTREME_WINDOW), ni-1)
        top             = ni
        qm%end_idx(i)   = sum(input_data_sorted( bottom:top )) / (top - bottom + 1)

        bottom          = min(max(int(match_i - match_step/4), nm - EXTREME_WINDOW), nm-1)
        top             = nm
        input_i         = sum(data_to_match_sorted( bottom:top ))/ (top - bottom + 1)
        denominator     = max(SMALL_VALUE, qm%end_idx(i) - qm%start_idx(i))
        qm%slope(i)     = (input_i - qm%offset(i)) &
                         / denominator
        qm%end_idx(i)   = input_data_sorted( ni )

    end subroutine develop_qm

    ! reverse the quantile map for a given input value
    function qm_value_reverse(input, qm) result(output)
        implicit none
        real, intent(in) :: input
        type(qm_correction_type), intent(in) :: qm

        real :: output
        integer :: i, n, step
        logical :: found

        found=.False.
        n = size(qm%start_idx)
        i = n/2
        step = n/4+1

        ! search for the first point in offset with start > input
        ! Uses a log(n) binary search because we "know" the qm data are sorted.
        do while (.not.found)
            i = max(1,min(n,i))
            if (qm%offset(i) > input) then
                if (qm%offset(max(1,i-1)) <= input) then
                    found=.True.
                else if (i==1) then
                    found=.True.
                else
                    i = i-step
                    step = ceiling(step/2.0)
                endif
            else if (i==n) then
                found=.True.
            else
                i = i+step
                step = ceiling(step/2.0)
            endif
        end do
        if (i>1) then
            if (qm%offset(i) > input) then
                i=i-1
            endif

        endif

        if (i==n) then
            output = qm%end_idx(i)
        elseif (qm%slope(i) /= 0) then
            output = (input - qm%offset(i)) / qm%slope(i) + qm%start_idx(i)
        else
            output = qm%start_idx(i)
        endif

    end function qm_value_reverse


    ! apply a quantile mapping scheme to a given input value
    function qm_value(input, qm) result(output)
        implicit none
        real, intent(in) :: input
        type(qm_correction_type), intent(in) :: qm

        real :: output
        integer :: i, n, step
        logical :: found

        found=.False.
        n = size(qm%start_idx)
        i = n/2
        step = n/4+1

        ! search for the first point in start_idx with start > input
        ! Uses a log(n) binary search because we "know" the qm data are sorted.
        do while (.not.found)
            i = max(1,min(n,i))
            if (qm%start_idx(i) > input) then
                if (qm%start_idx(max(1,i-1)) <= input) then
                    found=.True.
                else if (i==1) then
                    found=.True.
                else
                    i = i-step
                    step = ceiling(step/2.0)
                endif
            else if (i==n) then
                found=.True.
            else
                i = i+step
                step = ceiling(step/2.0)
            endif
        end do
        if (i>1) then
            if (qm%start_idx(i) > input) then
                i=i-1
            endif
        endif

        output = qm%offset(i) + qm%slope(i) * (input - qm%start_idx(i))

    end function qm_value

    !>----------------------------------------------
    !! Apply a previously developed quantile mapping to input and store the results in output
    !!
    !!----------------------------------------------
    subroutine apply_qm(input, output, qm)
        implicit none
        real, intent(in),  dimension(:) :: input
        real, intent(out), dimension(:) :: output
        type(qm_correction_type), intent(in) :: qm

        integer :: i, n

        n = size(input)
        !$omp parallel default(shared) &
        !$omp firstprivate(i,n)
        !$omp do
        do i = 1, n
            output(i) = qm_value(input(i), qm)
        end do
        !$omp end do
        !$omp end parallel

    end subroutine apply_qm

    !>----------------------------------------------
    !! Reverse a previously applied quantile mapping
    !!
    !!----------------------------------------------
    subroutine reverse_qm(input, qm)
        implicit none
        real, intent(inout),  dimension(:) :: input
        type(qm_correction_type), intent(in) :: qm

        integer :: i, n

        n = size(input)
        !$omp parallel default(shared) &
        !$omp firstprivate(i,n)
        !$omp do
        do i = 1, n
            input(i) = qm_value_reverse(input(i), qm)
        end do
        !$omp end do
        !$omp end parallel

    end subroutine reverse_qm

end module quantile_mapping
