module output_mod
    use io_routines, only : io_write
    use data_structures
    
    implicit none
contains
    subroutine write_output(output, options)
        implicit none
        type(config),  intent(in)   :: options
        type(results), intent(in)   :: output
        character(len=MAXFILELENGTH):: filename
        
        filename = options%output_file
        
        call io_write(filename, "data", output%variables(1)%data)
        
    end subroutine write_output

end module output_mod
