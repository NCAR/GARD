module analog_mod

    implicit none
    
    
contains

    function find_analogs(match, input, n) result(analogs)
        implicit none
        real,    intent(in), dimension(:)   :: match
        real,    intent(in), dimension(:,:) :: input
        integer, intent(in)     :: n
        integer, dimension(1:n)   :: analogs
        
        real,    dimension(:), allocatable :: distances
        logical, dimension(:), allocatable :: mask
        integer, dimension(1) :: min_location
        integer :: i, n_inputs, nvars
        
        ! if (match==0) then
        !     analogs = pick_n_random_zeros(input, n)
        ! else
            
            n_inputs = size(input,1)
            nvars    = size(input,2)

            allocate(mask(n_inputs))
            allocate(distances(n_inputs))
            distances = 0
            mask = .True.
            
            do i=1,nvars
                distances = distances + (input(:,i) - match(i))**2
            enddo
            
            ! fast way (O(nlogn) worst case, typically closer to O(n))
            call top_n_analogs(distances, n, analogs)
            
            ! foolish way (O(n^2))
            ! do i = 1, n
            !     min_location = minloc(distances,mask=mask)
            !     mask(min_location) = .false.
            !     
            !     analogs(i) = min_location(1)
            ! end do
        ! endif
        
    end function find_analogs
    
    ! compute the root mean square error between input[analogs] and y_hat
    function compute_analog_error(input, analogs, y_hat) result(error)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in)               :: y_hat
        real                              :: error
        
        double precision :: mean
        integer :: i, n
        
        n = size(analogs)
        mean = 0
        
        do i=1,n
            mean = mean + (input( analogs(i) ) - y_hat)**2
        enddo
        
        error = sqrt(mean / n)
        
    end function compute_analog_error
    
    
    function compute_analog_exceedance(input, analogs, threshold) result(probability)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in)               :: threshold
        real                              :: probability
        integer :: i, n
        
        n = size(analogs)
        probability = 0
        do i = 1, n
            if (input(analogs(i)) > threshold) then
                probability = probability + 1
            endif
        end do
        
        probability = probability / n
        
    end function compute_analog_exceedance

    ! private?
    subroutine top_n_analogs(distances, n, analogs)
        implicit none
        real,    intent(in),    dimension(:):: distances
        integer, intent(in)                 :: n
        integer, intent(inout), dimension(n):: analogs
        
        integer :: i, n_dists
        real, dimension(n) :: heap
        real :: inv_dist
        
        n_dists = size(distances)
        heap    = -9999
        analogs = 0
        
        do i = 1, n_dists
           inv_dist = 0-distances(i)
           ! top_dists is a minheap, if inv_dist < top_dists(1) then we need to replace it with inv_dist
           ! analogs is an index to the positions of each element of distances, and is
           ! maintained along with top_dsts by heapify
           if (inv_dist > heap(1)) then
               ! adding inv_dist(i) to the heap
               heap(1) = inv_dist
               analogs(1)   = i
               ! heapify re-heapifies the heap and the analogs arrays
               call heapify(heap, analogs, 0, n)
           endif
        enddo

    end subroutine top_n_analogs

    
    ! private? 
    function pick_n_random_zeros(input, n) result(analogs)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in)     :: n
        integer, dimension(1:n)   :: analogs
        integer, dimension(:), allocatable :: zero_locations
        
        integer :: i, n_inputs, n_zeros, selected_analog
        real :: rand

        n_inputs = size(input)
        allocate(zero_locations(n_inputs))
        n_zeros = 0
        do i=1,n_inputs
            if (input(i)==0) then
                n_zeros = n_zeros + 1
                zero_locations(n_zeros) = i
            end if
        end do
        
        do i=1,n
            call random_number(rand)
            selected_analog = floor(rand * n_zeros)+1
            analogs(i) = zero_locations(selected_analog)
        end do

    end function pick_n_random_zeros

    
    ! private? 
    subroutine heapify(a, ind, start, bottom)
        implicit none

        real,    intent(inout)  :: a(0:)
        integer, intent(inout)  :: ind(0:)
        integer, intent(in)     :: start, bottom
        
        integer :: child, root
        real    :: temp
        integer :: tempind

        root = start
        do while(root*2 + 1 < bottom)
            child = root * 2 + 1

            if (child + 1 < bottom) then
                if (a(child) > a(child+1)) child = child + 1
            end if

            if (a(root) > a(child)) then
                tempind     = ind(child)
                ind(child)  = ind(root)
                ind(root)   = tempind

                temp        = a(child)
                a(child)    = a (root)
                a(root)     = temp
                root        = child
            else
                return
            end if  
        end do      
    end subroutine heapify


    

end module analog_mod
