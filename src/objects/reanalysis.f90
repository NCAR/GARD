module reanalysis_mod

    use data_structures
    use geo
    
    implicit none
    
    interface reanalysis
        module procedure :: new_reanalysis
    end interface reanalysis
    
    type, public :: reanalysis
        private
        character(len=MAXFILELENGTH) :: name
      contains
        procedure, public  :: init        => init
        procedure, public  :: as_string   => as_string
        procedure, public  :: set_name    => set_name
    end type reanalysis

contains
    function new_reanalysis()
        implicit none
        type(reanalysis) :: new_reanalysis
        
        call new_reanalysis%set_name( "default-name" )
        
    end function new_reanalysis

    subroutine init(this, options)
        implicit none
        class(reanalysis), intent(inout) :: this
        class(options_type), intent(in)  :: options
        
        call this%set_name( options%training%name )
        
    end subroutine init
    
    function as_string(this) result(output)
        implicit none
        class(reanalysis), intent(in) :: this
        character(len=MAXVARLENGTH) :: output
        
        output = this%name
        
    end function as_string

    subroutine set_name(this, name)
        implicit none
        class(reanalysis), intent(inout) :: this
        character(len=*) :: name
        
        this%name = name
        
    end subroutine set_name
    
end module reanalysis_mod
