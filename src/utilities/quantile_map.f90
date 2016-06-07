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
    
    integer, parameter :: DEFAULT_QM_SEGMENTS = 100
    real, parameter :: SMALL_VALUE = 1e-20

contains

    subroutine develop_qm(input_data, data_to_match, qm, n_segments)
        implicit none
        real, intent(in), dimension(:) :: input_data, data_to_match
        integer, intent(in), optional :: n_segments
        type(qm_correction_type), intent(out) :: qm
        
        integer :: ni, nm, nseg
        integer :: i
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
            print*, "Quantile mapping with <10 segments not recommended, using 10 segments. "
            nseg=10
        endif
        input_step = ni/real(nseg)
        match_step = nm/real(nseg)
        print*, input_step, match_step

        input_i = 1
        match_i = 1

        if (allocated(qm%start_idx)) then 
            deallocate(qm%start_idx,qm%end_idx,qm%slope,qm%offset)
        endif
        allocate(qm%start_idx(nseg))
        allocate(qm%end_idx(nseg))
        allocate(qm%slope(nseg))
        allocate(qm%offset(nseg))
        
        do i=1,nseg-1
            qm%start_idx(i) = input_data_sorted( int(input_i) )
            qm%end_idx(i)   = input_data_sorted( int(input_i + input_step) )
            qm%offset(i)    = data_to_match_sorted( int(match_i) )
            denominator     = max(SMALL_VALUE, qm%end_idx(i) - qm%start_idx(i))
            qm%slope(i)     = (data_to_match_sorted( int(match_i + match_step) ) - qm%offset(i)) &
                            / denominator
            
            input_i = input_i + input_step
            match_i = match_i + match_step
        end do
        
        ! For the last segment, make sure it finishes at the final value
        ! this will not happen otherwise because nseg invariably will not divide nm and ni evenly
        i = nseg
        qm%start_idx(i) = input_data_sorted( input_i )
        qm%end_idx(i)   = input_data_sorted( ni )
        qm%offset(i)    = data_to_match_sorted( match_i )
        denominator = max(SMALL_VALUE, qm%end_idx(i) - qm%start_idx(i))
        qm%slope(i)     = (data_to_match_sorted( nm ) - qm%offset(i)) &
                        / denominator
        
        
    end subroutine develop_qm
    
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

end module quantile_mapping
