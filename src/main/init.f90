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

        write(*,*) "Generalized Analog Regression Downscaling (GARD)"
        write(*,*) "-----------------------------------------------------------"
        write(*,*) "GARD Version : "//trim(kVERSION_STRING)
        ! TODO: Add compile time options
        write(*,*) ""
        write(*,*) "  The Generalized Analog Regression Downscaling (GARD)"
        write(*,*) "  downscaling tool, version "//trim(kVERSION_STRING)//", Copyright (C) 2017 The"
        write(*,*) "  National Center for Atmospheric Research. GARD comes with"
        write(*,*) "  ABSOLUTELY NO WARRANTY. This is free software, you may "
        write(*,*) "  redistribute it under certain conditions; see LICENSE.txt"
        write(*,*) "  for details."
        write(*,*) ""
        write(*,*) "  Online Documentation      : http://gard.readthedocs.io"
        write(*,*) "  Report Bugs and Issues to : https://github.com/NCAR/GARD/issues"
        write(*,*) ""
        write(*,*) "-----------------------------------------------------------"
        write(*,*) "Usage: gard [-h] [--version] options_file"
        write(*,*) ""
        write(*,*) "-h              Help information for GARD"
        write(*,*) "--version       Print the version number"
        write(*,*) "options_file    Input options file name"
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
