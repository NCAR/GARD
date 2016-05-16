module init_mod
    use data_structures
    use model_constants
    use gcm_mod, only : init_gcm_io
    
    implicit none
    private
    public :: model_init
contains
    
    subroutine print_model_init(options)
        implicit none
        type(config), intent(in) :: options
        
        write(*,*) ""
        write(*,*) " ------------------------------------------------------ "
        write(*,*) ""
        write(*,*) "          Downscaling Model           "
        write(*,*) "          Version : "//trim(options%version)
        write(*,*) ""
        write(*,*) "  Developed at NCAR: "
        write(*,*) "   The National Center for Atmospheric Research "
        write(*,*) "   NCAR is supported by the National Science Foundation "
        write(*,*) " ------------------------------------------------------ "
        write(*,*) ""

        if (options%debug) then
            write(*,*) "Downscaling for the period : ", trim(options%first_time%as_string())
            write(*,*) "                        to : ", trim(options%last_time%as_string())
            write(*,*) "   Training for the period : ", trim(options%training_start%as_string())
            write(*,*) "                        to : ", trim(options%training_stop%as_string())
            write(*,*) ""
        endif
        
    end subroutine print_model_init
    
    
    subroutine model_init(options)
        implicit none
        type(config), intent(in) :: options
        
        call print_model_init(options)
        
        call init_gcm_io(options)
        
    end subroutine model_init

end module init_mod
    
    
