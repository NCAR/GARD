module time_io

    use data_structures
    use string
    use time
    use io_routines, only : io_read, io_read_attribute

    implicit none

contains

    function time_gain_from_units(units) result(gain)
        implicit none
        character(len=MAXSTRINGLENGTH), intent(in) :: units
        double precision :: gain

        if ((units(1:4)=="days").or.(units(1:4)=="Days")) then
            gain = 1.0
        else if ((units(1:4)=="hour").or.(units(1:4)=="Hour")) then
            gain = 1/24.0D0
        else if ((units(1:3)=="sec").or.(units(1:3)=="Sec")) then
            gain = 1/86400.0D0
        else if ((units(1:3)=="min").or.(units(1:3)=="Min")) then
            gain = 1/1440.0D0
        else
            write(*,*) trim(units)
            stop "Error: unknown units"
        endif

    end function time_gain_from_units

    function year_from_units(units) result(year)
        implicit none
        character(len=*), intent(in) :: units
        integer :: year

        integer :: since_loc, year_loc

        since_loc = index(units,"since")

        year_loc = index(units(since_loc:)," ")
        year_loc = year_loc+since_loc

        year = get_integer(units(year_loc:year_loc+3))

    end function year_from_units


    function month_from_units(units) result(month)
        implicit none
        character(len=*), intent(in) :: units
        integer :: month

        integer :: since_loc, month_loc

        since_loc = index(units,"since")

        month_loc = index(units(since_loc:)," ") + 5
        month_loc = month_loc + since_loc

        month = get_integer(units(month_loc:month_loc+1))

    end function month_from_units

    function day_from_units(units) result(day)
        implicit none
        character(len=*), intent(in) :: units
        integer :: day

        integer :: since_loc, day_loc

        since_loc = index(units,"since")

        day_loc = index(units(since_loc:)," ") + 8
        day_loc = day_loc + since_loc

        day = get_integer(units(day_loc:day_loc+1))

    end function day_from_units

    function hour_from_units(units, error) result(hour)
        implicit none
        character(len=*), intent(in) :: units
        integer, intent(out), optional :: error
        integer :: hour

        integer :: since_loc, hour_loc
        if (present(error)) error = 0

        since_loc = index(units,"since")

        hour_loc = index(units(since_loc:)," ") + 11
        hour_loc = hour_loc + since_loc

        ! default return value if hours can't be read from the units attribute (e.g. they aren't present)
        hour = 0

        if( hour_loc+1 <= len(units) ) then
            if (trim(units(hour_loc:hour_loc+1)) /= "") then
                hour = get_integer(units(hour_loc:hour_loc+1))
            endif
        else
            if (present(error)) error = 1
        endif

    end function hour_from_units

    function get_selected_time(options) result(selected_time)
        implicit none
        class(input_config), intent(in) :: options
        integer :: selected_time

        select type(options)
        class is (atm_config)
            if (options%selected_time/=-1) then
                selected_time = options%selected_time
            elseif (options%time_indices(1)/=-1) then
                selected_time = options%time_indices( ceiling(size(options%time_indices) / 2.0) )
            else
                selected_time = -1
            endif

        class default
            if (options%selected_time/=-1) then
                selected_time = options%selected_time
            else
                selected_time = -1
            endif
        end select

    end function get_selected_time

    subroutine read_times(options, times, timezone_offset)
        implicit none
        class(input_config), intent(in) :: options
        type(Time_type), intent(inout), dimension(:) :: times
        double precision, optional :: timezone_offset

        double precision, allocatable, dimension(:) :: temp_times
        integer :: ntimes, file_idx, cur_time, time_idx, error
        integer :: start_year, start_month, start_day, start_hour
        character(len=MAXSTRINGLENGTH) :: calendar, units
        double precision :: calendar_gain
        integer :: selected_time

        ntimes = size(times,1)
        cur_time = 1

        selected_time = get_selected_time(options)

        do file_idx = 1, size(options%file_names,1)

            ! first read the time variable (presumebly a 1D double precision array)
            call io_read(options%file_names(file_idx, options%time_file), options%time_name, temp_times)
            ! attempt to read the calendar attribute from the time variable
            call io_read_attribute(options%file_names(file_idx, options%time_file),"calendar", &
                                   calendar, var_name=options%time_name, error=error)
            ! if time attribute it not present, set calendar to one specified in the config file
            if (error/=0) then
                calendar = options%calendar
            endif

            ! attempt to read the units for this time variable
            call io_read_attribute(options%file_names(file_idx, options%time_file), "units", &
                                   units, var_name=options%time_name, error=error)

            if (error==0) then
                start_year = year_from_units(units)
                start_month   = month_from_units(units)
                start_day     = day_from_units(units)
                start_hour    = hour_from_units(units)
                calendar_gain = time_gain_from_units(units)
            else
                start_year = options%calendar_start_year
                calendar_gain = options%time_gain
            endif

            ! puts the units to days since ...
            ! in case it is in units of e.g. "hours since" or "seconds since"
            temp_times = temp_times * calendar_gain
            if (present(timezone_offset)) then
                temp_times = temp_times + timezone_offset / 24.0
            endif

            if (selected_time == -1) then
                do time_idx = 1, size(temp_times,1)

                    call times(cur_time)%init(calendar, start_year, start_month, start_day, start_hour)
                    call times(cur_time)%set(temp_times(time_idx))

                    cur_time = cur_time + 1
                end do
            else
                call times(cur_time)%init(calendar, start_year, start_month, start_day, start_hour)
                call times(cur_time)%set(temp_times(selected_time))

                cur_time = cur_time + 1
            endif

            deallocate(temp_times)
        end do

    end subroutine read_times

end module time_io
