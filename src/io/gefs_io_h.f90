!>------------------------------------------------
!! Handle all IO for GEFS data
!!
!! Loops through GEFS variables reading data and computing basic statistics
!! Reads time and lat/lon variables from the GEFS files as well.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module gefs_mod

    use data_structures

    implicit none

interface
    !>------------------------------------------------
    !! Initialize the GEFS module
    !!
    !!------------------------------------------------
    module subroutine init_GEFS_io(options)
        implicit none
        type(config), intent(in) :: options

    end subroutine init_GEFS_io

    !>------------------------------------------------
    !! Read in the GEFS data
    !!
    !! Uses variable names and filenames defined in the options data structure
    !!
    !! Loops through all input variables, then reads lat, lon, and time variables
    !!
    !!------------------------------------------------
    module function read_GEFS(options) result(GEFS_data)
        implicit none
        class(atm_config), intent(in) :: options
        type(atm) :: GEFS_data

    end function read_GEFS

end interface
end module GEFS_mod
