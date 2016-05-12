program test_master_file
    use master_file
    use time
    use string
    use data_structures
    
    implicit none
    
    type(Time_type) :: date
    character(len=MAXVARLENGTH) :: name
    integer :: j,n
    type(master_file_type) :: master
    
    name = " This is a test {Y} {m}"
    call master%init(name)
    
    call date%init("gregorian")
    call date%set("1999-01-20 10:01:03")
    
    print*, trim(date%as_string())
    print*, trim(str(date%mjd()))

    print*, trim(name)
    print*, "Accessing object..."
    print*, trim(master%as_string())
    print*, trim(master%get_file(date))
    
    ! call date%set("2020-01-20 10:00:03")

end program test_master_file
