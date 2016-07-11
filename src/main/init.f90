!>------------------------------------------------
!! Model Initialization
!! 
!! model_init is called after the configuration file has been read
!! This provides a place to print basic model information before running
!! This is also used to call initializers for any other modules that need it. 
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module init_mod
    use data_structures
    use model_constants
    use gefs_mod, only : init_gefs_io
    use gcm_mod, only : init_gcm_io
    use obs_mod, only : init_obs_io
    
    implicit none
    private
    public :: model_init
contains
    
    !>------------------------------------------------
    !! Prints model configuration info before running
    !!
    !! Prints a welcome and version string as well. 
    !! Only prints configuration info if debug == true
    !!
    !!------------------------------------------------
    subroutine print_model_init(options)
        implicit none
        type(config), intent(in) :: options
        
        write(*,*) ""
        write(*,*) " ====================================================== "
        write(*,*) ""
        write(*,*) "          Downscaling Model           "
        write(*,*) "          Version : "//trim(options%version)
        write(*,*) ""
        write(*,*) "  Developed at NCAR: "
        write(*,*) "   The National Center for Atmospheric Research "
        write(*,*) "   NCAR is supported by the National Science Foundation "
        write(*,*) " ------------------------------------------------------ "
        write(*,*) " ------------------------------------------------------ "
        write(*,*) ""
        write(*,*) "  WARNING WARNING WARNING WARNING WARNING WARNING "
        write(*,*) ""
        write(*,*) "    This is pre-release not-even-beta code. "
        write(*,*) ""
        write(*,*) "  WARNING WARNING WARNING WARNING WARNING WARNING "
        write(*,*) ""
        write(*,*) " ------------------------------------------------------ "
        write(*,*) " ------------------------------------------------------ "
        write(*,*) " ====================================================== "
        write(*,*) ""

        if (options%debug) then
            write(*,*) "Downscaling for the period : ", trim(options%first_time%as_string())
            write(*,*) "                        to : ", trim(options%last_time%as_string())
            write(*,*) "   Training for the period : ", trim(options%training_start%as_string())
            write(*,*) "                        to : ", trim(options%training_stop%as_string())
            write(*,*) ""
        endif
        
    end subroutine print_model_init
    
    
    !>------------------------------------------------
    !! Initialize the model
    !!
    !! Calls routines to print startup message and call module initializers
    !!
    !!------------------------------------------------
    subroutine model_init(options)
        implicit none
        type(config), intent(in) :: options
        
        call print_model_init(options)
        
        call init_gefs_io(options)
        call init_gcm_io(options)
        call init_obs_io(options)
        
    end subroutine model_init

end module init_mod
    
    
