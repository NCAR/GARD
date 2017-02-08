module analog_mod

    use string, only: str
    implicit none

    integer, parameter :: MIN_NUMBER_ANALOGS = 20


contains

    subroutine find_analogs(analogs, match, input, n, threshold, weights, skip_analog)
        implicit none
        integer, intent(inout), dimension(:), allocatable  :: analogs
        real,    intent(in),    dimension(:)   :: match
        real,    intent(in),    dimension(:,:) :: input
        integer, intent(in)                    :: n
        real,    intent(in)                    :: threshold
        real,    intent(inout), dimension(:), allocatable, optional :: weights
        integer, intent(in),    optional       :: skip_analog ! specify an element of the input data to skip

        real,    dimension(:), allocatable :: distances
        ! logical, dimension(:), allocatable :: mask
        integer, dimension(1) :: min_location
        integer :: i, n_inputs, nvars

        n_inputs = size(input,1)
        nvars    = size(input,2)

        ! allocate(mask(n_inputs))
        allocate(distances(n_inputs))
        distances = 0
        ! mask = .True.

        do i=1,nvars
            distances = distances + (input(:,i) - match(i))**2
        enddo
        if (present(skip_analog)) then
            if ((skip_analog > 0).and.(skip_analog < size(distances))) then
                distances(skip_analog) = distances(skip_analog) + 1000 ! 1 would be a large value
            endif
        endif

        ! fast selection (O(n_inputs log(nanalogs)) worst case, typically closer to O(n_inputs))
        if (n>0) then
            if (.not.allocated(analogs)) allocate(analogs(n))
            call top_n_analogs(distances, n, analogs)
        elseif (threshold>0) then
            distances = distances / nvars ! normalize by the number of variables so the supplied threshold doesn't have to changes
            call analogs_from_threshold(distances, threshold, analogs)
        else
            !$omp critical (print_lock)
            write(*,*) "ERROR: find_analogs, nanalogs=",trim(str(n))," analog_threshold=",trim(str(threshold))
            stop "No number of or threshold for analog selection supplied"
            !$omp end critical (print_lock)
        endif

        if (minval(analogs)<1) then
            !$omp critical (print_lock)
            write(*,*) ""
            write(*,'(A,I4)'         ) "   n analogs:",n
            write(*,'(A,I7,A,I7)'    ) "   Minimum analog index:",   minval(analogs),   "      Maximum analog index:",maxval(analogs)
            write(*,'(A,4I7)'        ) "   First 4 analogs:", analogs(1:4)
            write(*,'(A,F9.1,A,F9.1)') "   Maximum analog distance:",maxval(distances), "      Minimum analog distance:",minval(distances)
            write(*,*) "   Predictor values to match:", match
            write(*,*) ""
            write(*,*) " Are training and predictor transforms set appropriately?"
            write(*,*) ""
            stop "ERROR selecting analogs"
            !$omp end critical (print_lock)
        endif

        if (present(weights)) then
            if (allocated(weights)) then
                if (size(weights)/=size(analogs)) then
                    deallocate(weights)
                    allocate( weights(size(analogs)) )
                endif
            else
                allocate( weights(size(analogs)) )
            endif
            do i=1,size(analogs)
                if (distances(analogs(i)) > (1e-2)) then
                    weights(i) = 1 / distances(analogs(i))
                else
                    weights(i) = 1 / (1e-2)
                endif
            enddo
            weights = weights / sum(weights)
        endif
    end subroutine find_analogs

    ! compute the mean of input[analogs]
    function compute_analog_mean(input, analogs, weights) result(mean)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in), dimension(:), optional :: weights
        real                              :: mean

        double precision :: internal_mean ! use double precision internally to avoid errors when summing many numbers
        integer :: i, n

        n = size(analogs)
        internal_mean = 0

        if (present(weights)) then
            do i=1, n
                internal_mean = internal_mean + input( analogs(i) ) * weights(i)
            enddo
            mean = internal_mean
        else
            do i=1, n
                internal_mean = internal_mean + input( analogs(i) )
            enddo
            mean = internal_mean / n
        endif


    end function compute_analog_mean

    ! compute the root mean square error between input[analogs] and y_hat
    function compute_analog_error(input, analogs, y_hat, weights) result(error)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in)               :: y_hat
        real,    intent(in), dimension(:), optional :: weights
        real                              :: error

        double precision :: mean
        integer :: i, n

        n = size(analogs)
        mean = 0

        if (present(weights)) then
            do i=1, n
                mean = mean + (input( analogs(i) ) - y_hat)**2 * weights(i)
            enddo
            error = sqrt(mean)
        else
            do i=1,n
                mean = mean + (input( analogs(i) ) - y_hat)**2
            enddo

            error = sqrt(mean / n)
        endif

    end function compute_analog_error


    function compute_analog_exceedance(input, analogs, threshold, weights) result(probability)
        implicit none
        real,    intent(in), dimension(:) :: input
        integer, intent(in), dimension(:) :: analogs
        real,    intent(in)               :: threshold
        real,    intent(in), dimension(:), optional :: weights
        real                              :: probability
        integer :: i, n

        n = size(analogs)
        probability = 0
        if (present(weights)) then
            do i = 1, n
                if (input(analogs(i)) > threshold) then
                    probability = probability + weights(i)
                endif
            end do
        else
            do i = 1, n
                if (input(analogs(i)) > threshold) then
                    probability = probability + 1
                endif
            end do

            probability = probability / n
        endif

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
           ! top_dists is a minheap, if inv_dist > top_dists(1) then we need to replace it with inv_dist
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

    ! private? search for all analog days that match a given threshold
    subroutine analogs_from_threshold(distances, initial_threshold, analogs)
        implicit none
        real,    intent(in),    dimension(:)               :: distances
        real,    intent(in)                                :: initial_threshold
        integer, intent(inout), dimension(:), allocatable  :: analogs

        real :: threshold
        integer :: i, n_distances, n_analogs, current_analog
        n_distances = size(distances)

        threshold = initial_threshold
        n_analogs = 0
        do while (n_analogs < MIN_NUMBER_ANALOGS)
            n_analogs = 0
            do i=1,n_distances
                if (distances(i) < threshold) then
                    n_analogs = n_analogs + 1
                endif
            enddo
            ! in case we did not find enough analogs, increase the search threshold before trying again
            threshold = threshold * 2
        enddo

        if (allocated(analogs)) then
            deallocate(analogs)
        endif
        allocate(analogs(n_analogs))
        threshold = threshold / 2 ! bump the threshold back to where it was when we found these analogs
        current_analog = 1
        do i=1,n_distances
            if (distances(i) < threshold) then
                analogs(current_analog) = i
                current_analog = current_analog + 1
            endif
        enddo

    end subroutine analogs_from_threshold

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
