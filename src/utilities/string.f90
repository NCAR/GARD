!>------------------------------------------------------------
!!  Various functions to convert a number to a string and a string
!!  to a number. 
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------------------
module string

    implicit none
    interface str
        module procedure str_d
        module procedure str_r
        module procedure str_i
    end interface
    
    integer,parameter::MAXSTRINGLENGTH=100
contains
    
    ! could this be done with index(name(start:),token)?
    ! function find_token(name, start, token) result(i)
    !     implicit none
    !     character(len=*), intent(in) :: name
    !     integer,          intent(in), optional :: start
    !     character(len=1), intent(in), optional :: token
    !     integer :: i
    !     logical :: done
    !     character(len=1) :: search
    !     
    !     if ( present(start) ) then
    !         i = start
    !     else
    !         i = 1            
    !     end if
    !     
    !     if ( present(token) ) then
    !         search = token
    !     else
    !         search = "{"
    !     end if
    !     
    !     done=.False.
    !     if (name(i:i)/="{") then
    !         
    !         do while ( .not. done )
    !             i=i+1
    !             if (name(i:i)=="{") then
    !                 done=.True.
    !             endif
    !         end do
    !         
    !     endif
    ! 
    ! end subroutine find_token
    ! 
    function get_double(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        double precision :: get_double
        read(str_in,*) get_double
    end function get_double
    
    function get_real(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        real :: get_real
        read(str_in,*) get_real
    end function get_real

    function get_integer(str_in)
        implicit none
        character(len=*), intent(in) :: str_in
        integer :: get_integer
        read(str_in,*) get_integer
    end function get_integer
    
    
    
    function str_d(value,fmt) result(output_string)
        implicit none
        double precision :: value
        character(len=*), optional :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_d

    function str_r(value,fmt) result(output_string)
        implicit none
        real :: value
        character(len=*), optional :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
    end function str_r

    function str_i(value,fmt, length, pad) result(output_string)
        implicit none
        integer :: value
        character(len=*), optional :: fmt
        character(len=MAXSTRINGLENGTH) :: output_string
        integer, intent(in), optional :: length
        character(len=1), intent(in), optional :: pad
        integer :: extra
        character(len=MAXSTRINGLENGTH) :: temporary
        
        if (present(fmt)) then
            write(output_string,fmt) value
        else
            write(output_string,*) value
        endif
        output_string=adjustl(output_string)
        
        if (present(length)) then
            if (len_trim(output_string)/=length) then
                extra = length - len_trim(output_string)
                temporary(1:extra) = pad
                temporary(extra:length) = trim(output_string)
                output_string = temporary
            endif
        endif
                
    end function str_i


end module string
