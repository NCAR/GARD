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

    implicit none
    private
    public :: develop_qm
    public :: apply_qm
    
    integer, parameter :: DEFAULT_QM_SEGMENTS = 100

contains
    ! type qm_correction_type
    !     integer, allocatable, dimension(:) :: start_idx, end_idx
    !     real, allocatable, dimension(:) :: slope, offset
    ! end type qm_correction_type

    function develop_qm(input_data, data_to_match, n_segments) result(qm)
        implicit none
        real, intent(in), dimension(:) :: input_data, data_to_match
        integer, intent(in), optional :: n_segments
        type(qm_correction_type) :: qm
        
        integer :: ni, nm, nseg
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
        match_i = 0

        allocate(qm%start_idx(nseg))
        allocate(qm%end_idx(nseg))
        allocate(qm%slope(nseg))
        allocate(qm%offset(nseg))
        
        do i=1,nseg
            qm%start_idx(i) = input_data_sorted( input_i )
            qm%end_idx(i)   = input_data_sorted( input_i + input_step)
            qm%offset(i)    = data_to_match_sorted( match_i )
            qm%slope(i)     = (data_to_match_sorted( match_i + match_step) - qm%offset(i)) &
                            / (qm%end_idx(i) - qm%end_idx(i))
        end do
        
    end function develop_qm
    
    ! apply a quantile mapping scheme to a given input value
    function apply_qm(input, qm) result(output)
        implicit none
        real, intent(in) :: input
        type(qm_correction_type), intent(in) :: qm
        
        real :: output
        integer :: i
        logical :: found
        
        found=.False.
        i=0
        
        ! search for the first point in start_idx with start > input
        ! if none then apply the last 
        do while (.not.found)
            if (qm%start_idx(i) <= input) then
                found=.True.
            else
                i = i+1
            endif
            if (i>size(qm%start_idx,1)) then
                found=.True.
            endif
        end do
        i=i-1
        
        output = qm%offset(i) + qm%slope(i) * (input - qm%start_idx(i))
        
    end function apply_qm

end module quantile_mapping
