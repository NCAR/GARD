module output_mod

    use data_structures

    implicit none

interface
    module subroutine write_output(output, options)
        implicit none
        type(config),  intent(in)   :: options
        type(results), intent(in)   :: output

    end subroutine write_output

end interface
end module output_mod
