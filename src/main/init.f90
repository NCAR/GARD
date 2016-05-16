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
        
    end subroutine print_model_init
    
    
    subroutine model_init(options)
        implicit none
        type(config), intent(in) :: options
        
        call print_model_init(options)
        
        call init_gcm_io(options)
        
    end subroutine model_init

end module init_mod
    
    
