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
    public :: model_init, print_model_init
contains

    !>------------------------------------------------
    !! Prints model configuration info before running
    !!
    !! Prints a welcome and version string as well.
    !! Only prints configuration info if debug == true
    !!
    !!------------------------------------------------
    subroutine print_model_init()
        implicit none

        write(*,*) ""
        write(*,*) " ====================================================== "
        write(*,*) ""
        write(*,*) "          Downscaling Model           "
        write(*,*) "          Version : "//trim(kVERSION_STRING)
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

        call init_gefs_io(options)
        call init_gcm_io(options)
        call init_obs_io(options)

    end subroutine model_init

end module init_mod
