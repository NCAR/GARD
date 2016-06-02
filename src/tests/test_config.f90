program test_config
    use config_mod,     only : read_files_list, read_data_type, get_options_file
    use data_structures
    use model_constants
    
    implicit none
    
    character(len=MAXFILELENGTH) :: opt_file
    character(len=MAXFILELENGTH), dimension(MAX_NUMBER_FILES,3) :: file_list
    
    character(len=MAXFILELENGTH) :: filename
    integer, parameter :: correct_nfiles = 23
    character(len=MAXFILELENGTH), dimension(correct_nfiles) :: correct_file_list
    
    integer :: i, nfiles
    logical :: passing
    
    correct_file_list = [character(len=MAX_NUMBER_FILES) :: &
                            "config/configuration.f90",     &
                            "io/atmosphere_io.f90",         &
                            "io/gcm_io.f90",                &
                            "io/gefs_io.f90",               &
                            "io/io_routines.f90",           &
                            "io/obs_io.f90",                &
                            "io/time_io.f90",               &
                            "main/init.f90",                &
                            "main/main.f90",                &
                            "objects/data_structures.f90",  &
                            "objects/model_constants.f90",  &
                            "objects/reanalysis.f90",       &
                            "tests/test_calendar.f90",      &
                            "tests/test_config.f90",        &
                            "tests/test_gcm_io.f90",        &
                            "tests/test_qm.f90",            &
                            "tests/test_sort.f90",          &
                            "utilities/basic_stats.f90",    &
                            "utilities/geo_reader.f90",     &
                            "utilities/quantile_map.f90",   &
                            "utilities/sorting.f90",        &
                            "utilities/string.f90",         &
                            "utilities/time.f90"            &
    ]
    

    passing = .True.
    write(*,*) "--------------------------"
    write(*,*) "Testing read_file_list:"
    filename = "../test_data/file_list.txt"
    inquire(file=trim(filename), exist=passing)
    if (passing) then
        nfiles = read_files_list(filename, file_list(:,2))
        if (nfiles /= correct_nfiles) passing = .False.
        do i = 1, nfiles
            if (trim(file_list(i,2)) /= trim(correct_file_list(i))) passing = .False.
        end do
    else
        write(*,*) "Unable to test because the test file does not exist: '", trim(filename),"'"
    endif
    call check_passing(passing)
    write(*,*) "--------------------------"
    

    passing = .True.
    write(*,*) "--------------------------"
    write(*,*) "Testing read_data_type:"
    if (read_data_type("GEFS") /= kGEFS_TYPE) passing=.False. 
    if (read_data_type("GCM") /= kGCM_TYPE)   passing=.False. 
    if (read_data_type("obs") /= kobs_TYPE)   passing=.False. 
    call check_passing(passing)
    write(*,*) "--------------------------"

    write(*,*) "--------------------------"
    write(*,*) " Please specify a valid filename (for a file that exists) e.g. 'makefile' on the commandline"
    opt_file = get_options_file()
    write(*,*) " Read commandline options file to be : '", trim(opt_file), "'"
    write(*,*) " If this is what you specified on the commandline, then the test passed. "
    write(*,*) "--------------------------"
    
contains
    subroutine check_passing(passing)
        implicit none
        logical :: passing
        
        if (passing) then
            write(*,*) "  PASSED"
        else
            write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            write(*,*) "                              FAILED"
            write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        endif
    end subroutine check_passing
    
end program test_config
