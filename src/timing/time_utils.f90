module time_util
    use data_structures

    implicit none

contains
    subroutine setup_time_indices(input, options)
        implicit none
        class(base_data_type), intent(inout) :: input
        class(config), intent(in) :: options

        input%first_time      = find_nearest_time(input%times, options%first_time)
        input%last_time       = find_nearest_time(input%times, options%last_time)

        input%training_start  = find_nearest_time(input%times, options%training_start)
        input%training_stop   = find_nearest_time(input%times, options%training_stop)

        input%transform_start = find_nearest_time(input%times, options%transform_start)
        input%transform_stop  = find_nearest_time(input%times, options%transform_stop)

        input%post_start      = find_nearest_time(input%times, options%post_start)
        input%post_end        = find_nearest_time(input%times, options%post_end)

    end subroutine setup_time_indices

    function find_nearest_time(time_list, match_time) result(nearest)
        implicit none
        type(Time_type), dimension(:), intent(in) :: time_list
        type(Time_type),               intent(in) :: match_time

        integer :: nearest, lower, upper, search, n
        logical :: found

        n = size(time_list)
        ! lower and upper define the bounds of the search region
        lower = 1
        upper = n
        found = .False.
        nearest = -1

        do while (.not.found)
            ! the current search is the middle of the search region
            search  = ((lower+upper)/2)

            ! perform a basic binary search
            if (time_list(search) < match_time) then
                lower = search+1
            else if (time_list(search) > match_time) then
                upper = search-1
            else if (time_list(search) == match_time) then
                found = .True.
                nearest = search
            else
                write(*,*) "-----------------------------------------------------"
                write(*,*) "broken..."
                write(*,*) "Searching for ",trim(match_time%as_string())
                write(*,*) "On element: ",trim(time_list(search)%as_string())
                write(*,*) "Seems to be neither < or > or == ..."
                write(*,*) "-----------------------------------------------------"
            endif
            if (lower>=upper) found = .True.
        end do

        if (nearest == -1) then
            ! we didn't actually find it, we exited because our search region decreased to 0
            ! this is likely to happen if the times don't line up exactly.
            if (time_list(search) > match_time) then
                search = max(1,search - 1)
                do while ((time_list(search) > match_time).and.(search /= 1))
                    search = max(1,search - 1)
                end do
            endif
            ! keep this separate from the above if (instead of as an else if) to ensure that the found time is always
            ! after the match time (or equal to it)
            if (time_list(search) < match_time) then
                search = min(n,search + 1)
                do while ((time_list(search) < match_time).and.(search /= n))
                    search = min(n,search + 1)
                end do
            endif
        endif
        nearest = search

    end function find_nearest_time


end module time_util
