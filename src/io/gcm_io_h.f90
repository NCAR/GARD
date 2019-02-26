!>------------------------------------------------
!! Handle all IO for GCM data
!!
!! Loops through GCM variables reading data and computing basic statistics
!! Reads time and lat/lon variables from the GCM files as well.
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!!------------------------------------------------
module gcm_mod

    use data_structures

    implicit none

interface

    !>------------------------------------------------
    !! Initialize the GCM module
    !!
    !!------------------------------------------------
    module subroutine init_gcm_io(options)
        implicit none
        type(config), intent(in) :: options

    end subroutine init_gcm_io

    !>------------------------------------------------
    !! Read in the GCM data
    !!
    !! Uses variable names and filenames defined in the options data structure
    !!
    !! Loops through all input variables, then reads lat, lon, and time variables
    !!
    !!------------------------------------------------
    module function read_gcm(options) result(gcm_data)
        implicit none
        class(atm_config), intent(in) :: options
        type(atm) :: gcm_data

    end function read_gcm

end interface
end module gcm_mod
