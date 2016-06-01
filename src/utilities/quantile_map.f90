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

contains

    subroutine develop_qm(input_data, data_to_match, qm, n_segments)
        implicit none
        real, intent(in), dimension(:) :: input_data, data_to_match
        integer, intent(in), optional :: n_segments
        type(qm_correction_type), intent(out) :: qm
        
        integer :: ni, nm, nseg
        integer :: i, input_i, input_step, match_i, match_step
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
        
        input_step = ni/nseg
        match_step = nm/nseg

        input_i = 1
        match_i = 1

        if (allocated(qm%start_idx)) then 
            deallocate(qm%start_idx,qm%end_idx,qm%slope,qm%offset)
        endif
        allocate(qm%start_idx(nseg))
        allocate(qm%end_idx(nseg))
        allocate(qm%slope(nseg))
        allocate(qm%offset(nseg))
        
        do i=1,nseg
            qm%start_idx(i) = input_data_sorted( input_i )
            qm%end_idx(i)   = input_data_sorted( min(ni, input_i + input_step))
            qm%offset(i)    = data_to_match_sorted( match_i )
            qm%slope(i)     = (data_to_match_sorted( min(nm, match_i + match_step)) - qm%offset(i)) &
                            / (qm%end_idx(i) - qm%start_idx(i))
            
            input_i = input_i + input_step
            match_i = match_i + match_step
        end do
        
    end subroutine develop_qm
    
    ! apply a quantile mapping scheme to a given input value
    function qm_value(input, qm) result(output)
        implicit none
        real, intent(in) :: input
        type(qm_correction_type), intent(in) :: qm
        
        real :: output
        integer :: i
        logical :: found
        
        found=.False.
        i=1
        
        ! search for the first point in start_idx with start > input
        ! if none then apply the last 
        do while (.not.found)
            if (qm%start_idx(i) > input) then
                found=.True.
            else
                i = i+1
            endif
            if (i>size(qm%start_idx,1)) then
                found=.True.
            endif
        end do
        if (i>1) then
            i=i-1
        endif
        
        output = qm%offset(i) + qm%slope(i) * (input - qm%start_idx(i))
        
    end function qm_value
    
    subroutine apply_qm(input, output, qm)
        implicit none
        real, intent(in),  dimension(:) :: input
        real, intent(out), dimension(:) :: output
        type(qm_correction_type), intent(in) :: qm
        
        integer :: i, n
        
        n = size(input)
        do i = 1, n
            output(i) = qm_value(input(i), qm)
        end do
        
    end subroutine apply_qm

end module quantile_mapping
