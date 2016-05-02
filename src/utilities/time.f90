!>------------------------------------------------------------
!!  date / time module
!!
!!  Defines a Time_type object
!!  
!!  Contains various utilities for working with dates and times. 
!!  Utilities are capable of handling multiple model calendars. 
!!  
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module time
    implicit none
    private

    integer, parameter, public :: GREGORIAN=0, NOLEAP=1, THREESIXTY=2
    
    integer, parameter :: MAXSTRINGLENGTH = 1024
    
    type, public :: Time_type
        integer :: year_zero = 1800  ! starting year
        integer :: calendar
        integer, dimension(13) :: month_start
        integer :: year, month, day, hour, minute, second
        
        double precision :: current_date_time = 0

      contains
        procedure, public  :: init        => time_init
        procedure, public  :: date        => calendar_date
        procedure, public  :: mjd         => get_mjd
        procedure, public  :: day_of_year => calc_day_of_year
        procedure, public  :: date_to_mjd => date_to_mjd
        procedure, public  :: as_string   => as_string
        
        generic,   public  :: set         => set_from_string, set_from_date, set_from_mjd
        procedure, private :: set_from_string, set_from_date, set_from_mjd
    end type Time_type
    
contains

    subroutine time_init(this, calendar_name, year_zero)
        implicit none
        class(Time_type) :: this
        character(len=*), intent(in) :: calendar_name
        integer, intent(in), optional :: year_zero
        
        integer :: i
        
        ! zero based month_starts (will have 1 added below)
        this%month_start = [0,31,59,90,120,151,181,212,243,273,304,334,365]
        
        select case (calendar_name)
            case("gregorian")
                this%calendar = GREGORIAN
            case("standard")
                this%calendar = GREGORIAN
            case("365-day")
                this%calendar = NOLEAP
            case("noleap")
                this%calendar = NOLEAP
            case("360-day")
                this%calendar = THREESIXTY
            case default
                write(*,*) "Unknown Calendar: ", trim(calendar_name)
                write(*,*) "Acceptable names = "
                write(*,*) "  gregorian, standard, 365-day, noleap, 360-day"
                write(*,*) " "
                stop
        end select
        
        if (this%calendar == THREESIXTY) then
            do i=0,12
                this%month_start(i+1) = i*30
            end do
        endif
        do i=0,12
            this%month_start(i+1) = this%month_start(i+1) + 1
        end do
        
        if ( present(year_zero) ) then
            this%year_zero = year_zero
        endif
        
    end subroutine time_init
    
    function get_mjd(this) result(mjd)
        implicit none
        class(Time_type) :: this
        double precision :: mjd
        
        mjd = this%current_date_time
    end function get_mjd
    
    !   algorithms from Wikipedia: http://en.wikipedia.org/wiki/Julian_day
    !   originally from Richards, E. G. (2013). Calendars. In S. E. Urban & P. K. Seidelmann, eds. 
    !                   Explanatory Supplement to the Astronomical Almanac, 3rd ed. (pp. 585–624). 
    !                   Mill Valley, Calif.: University Science Books. ISBN 978-1-89138-985-6
    !                   p617-9
    function date_to_mjd(this, year, month, day, hour, minute, second)
        implicit none
        class(Time_type), intent(in) :: this
        integer, intent(in) :: year, month, day, hour, minute, second
        double precision :: date_to_mjd
        double precision :: d,m,y
        integer :: a,b

        if (this%calendar==GREGORIAN) then
            a = (14-month)/12
            y = year+4800-a
            m = month+12*a-3
            ! Gregorian calendar
            b = day + floor((153*m+2)/5) + 365*y + floor(y/4) - floor(y/100) + floor(y/400) - 32045
            ! Julian calendar
            ! b = day + floor(153*m+2/5) + 365*y + floor(y/4) - 32083
            date_to_mjd = b + (((second/60d+0)+minute)/60d+0 + hour-12)/24.0 - 2400000.5
        else if (this%calendar==NOLEAP) then
            date_to_mjd = (year-this%year_zero)*365 + this%month_start(month)-1 + day-1 + (hour + (minute+second/60d+0)/60d+0)/24d+0
        else if (this%calendar==THREESIXTY) then
            date_to_mjd = (year-this%year_zero)*360 + this%month_start(month)-1 + day-1 + (hour + (minute+second/60d+0)/60d+0)/24d+0
        end if

    end function date_to_mjd

    ! compute the year, month, day, hour, minute, second corresponding
    ! to the input modified julian day (mjd)
    ! note mjd for NOLEAP and 360day calendars is not a true MJD
    ! arguably, seconds should be a real number, not an integer...
    subroutine calendar_date(this, year, month, day, hour, minute, second)
        implicit none
        class(Time_type), intent(in) :: this
        integer, intent(out) :: year, month, day, hour, minute, second
        
        integer :: y=4716,j=1401,m=2,n=12,r=4,p=1461
        integer :: v=3,u=5,s=153,w=2,B=274277,C=-38
        integer ::f,e,g,h, jday
        double precision :: day_fraction, mjd
        
        mjd = this%current_date_time+1d-5 ! add less than one second
        
        if (this%calendar==GREGORIAN) then
            jday=nint(mjd+2400000.5)
            f=jday+j+(((4*jday+B)/146097)*3)/4+C
            e=r*f+v
            g=mod(e,p)/r
            h=u*g+w
            day=mod(h,s)/u+1
            month=mod(h/s+m,n)+1
            year=e/p-y+(n+m-month)/n
            
        else if (this%calendar==NOLEAP) then
            year=floor(mjd/365)
            day_fraction=mjd - year*365+1
            do f=1,12
                if (day_fraction>this%month_start(f)) then
                    month=f
                endif
            end do
            day = floor(day_fraction - this%month_start(month))+1
            year=year+this%year_zero
            
        else if (this%calendar==THREESIXTY) then
            year=floor(mjd/360)
            day_fraction=mjd - year*360+1
            do f=1,12
                if (day_fraction>this%month_start(f)) then
                    month=f
                endif
            end do
            day = floor(day_fraction - this%month_start(month))+1
            year=year+this%year_zero
        end if
        
        day_fraction=mod(mjd,1.0)
        hour=floor(day_fraction*24)
        
        day_fraction=day_fraction*24-hour
        minute=floor(day_fraction*60)
        
        day_fraction=day_fraction*60-minute
        second = nint((day_fraction-(24d0*60*1d-5))*60)
        
    end subroutine calendar_date

    ! calculate the day of the year from a "modified julian day"
    ! note mjd for NOLEAP and 360day calendars is not a true MJD
    function calc_day_of_year(this)
        implicit none
        class(Time_type) :: this
        real :: calc_day_of_year
        
        integer :: year, month, day, hour, minute, second
        
        select case (this%calendar)
            case (GREGORIAN)
                call this%date(year, month, day, hour, minute, second)
                calc_day_of_year = this%current_date_time - this%date_to_mjd(year, 1,1,0,0,0)
            case (NOLEAP)
                calc_day_of_year = mod(this%current_date_time,365.0)
            case (THREESIXTY)
                calc_day_of_year = mod(this%current_date_time,360.0)
        end select
        
    end function calc_day_of_year

    ! convert an input date string in the form YYYY/MM/DD or YYYY/MM/DD hh:mm:ss
    ! into integer variables
    subroutine set_from_string(this, date)
        implicit none
        class(Time_type), intent(inout) :: this
        character (len=*), intent(in) :: date

        read(date(1:4), *) this%year
        read(date(6:7), *) this%month
        read(date(9:10),*) this%day
        
        if(len_trim(date) <= 18) then
            this%second  = 0
            this%minute  = 0
            this%hour    = 0
        else
            read(date(12:13), *) this%hour
            read(date(15:16), *) this%minute
            read(date(18:19), *) this%second
        endif
        
        this%current_date_time = this%date_to_mjd(this%year, this%month, this%day, this%hour, this%minute, this%second)
    end subroutine set_from_string
    
    subroutine set_from_date(this, year, month, day, hour, minute, second)
        implicit none
        class(Time_type), intent(inout) :: this
        integer, intent(in) :: year, month, day, hour, minute, second
                
        this%current_date_time = this%date_to_mjd(year, month, day, hour, minute, second)
        this%year   = year
        this%month  = month
        this%day    = day
        this%hour   = hour
        this%minute = minute
        this%second = second
        
    end subroutine set_from_date

    subroutine set_from_mjd(this, mjd)
        implicit none
        class(Time_type), intent(inout) :: this
        double precision, intent(in) :: mjd
        integer :: year, month, day, hour, minute, second
        
        this%current_date_time = mjd
        call this%date(year, month, day, hour, minute, second)
        this%year   = year
        this%month  = month
        this%day    = day
        this%hour   = hour
        this%minute = minute
        this%second = second
        
    end subroutine set_from_mjd
    
    ! Convert the date object into a string in the 0-filled format : "YYYY/MM/DD hh:mm:ss"
    function as_string(this, input_format) result(pretty_string)
        implicit none
        class(Time_type), intent(in) :: this
        character(len=*), intent(in), optional :: input_format
        character(len=MAXSTRINGLENGTH) :: pretty_string
        character(len=MAXSTRINGLENGTH) :: format
        integer :: i
        
        associate(year  => this%year,   &
                  month => this%month,  &
                  day   => this%day,    &
                  hour  => this%hour,   &
                  minute=> this%minute, &
                  second=> this%second  &
                  )
        
        if (present(input_format)) then
            format = input_format
        else
            ! this is the default format string to generate "YYYY/MM/DD hh:mm:ss"
            format = '(I4,A1,I2,A1,I2,A1,I2,A1,I2,A1,I2)'
        endif
        
        ! this and the format statement above are the important bits
        write(pretty_string, format) year,'/',month,'/',day,'=',hour,':',minute,":",second
        
        ! fill missing digits with '0'
        do i=1,len_trim(pretty_string)
            select case (pretty_string(i:i))
                case (' ')
                    pretty_string(i:i) = '0'
                case ('=')
                    pretty_string(i:i) = ' '
            end select
        end do
        
        end associate

    end function as_string
    
end module time
