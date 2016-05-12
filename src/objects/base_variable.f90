module base_variable_mod

    use data_structures
    use io_routines, only   : io_read
    use geo, only           : geo_LUT, geointerp
    use master_file, only   : master_file_type
    
    implicit none
    
    type, public :: variable
        private
        character(len=MAXFILELENGTH) :: name ! variable name to read from netCDF files
        
        logical :: have_files
        
        integer :: nfiles
        character(len=MAXFILELENGTH), dimension(:), allocatable :: file_list
        
        type(master_file_type) :: base_file
        
      contains
        procedure, public  :: init        => init
        procedure, public  :: as_string   => as_string
    end type variable
    interface variable
        module procedure new_var
    end interface variable

contains
    function new_var(options)
        implicit none
        class(var_config_type), intent(in)  :: options
        type(variable)                      :: new_var ! function result
        
        call new_var%init( options )
        
    end function new_var
    
    subroutine init(this, options)
        implicit none
        class(variable),     intent(inout)  :: this
        class(var_config_type), intent(in)  :: options
        
        call this%name = options%name
        
        ! setup either the base file naming scheme or the list of files
        if (options%use_filelist) then
            allocate( this%file_list(options%nfiles) )
            this%file_list = options%file_list
            this%nfiles    = options%nfiles
            this%have_files = .True.
        else
            this%base_file = master_file_type( options%base_file )
            this%have_files = .False.
        endif
        
        this%lat_var = options%lat_var
        this%lon_var = options%lon_var
        
        
        
    end subroutine init

    function as_string(this) result(output)
        implicit none
        class(variable), intent(in) :: this
        character(len=MAXVARLENGTH) :: output
        
        output = this%name
        
    end function as_string

    
end module base_variable_mod
