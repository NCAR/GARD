module gcm_mod

    use data_structures
    use geo
    
    implicit none
    
    interface gcm
        module procedure :: new_gcm
    end interface gcm
    private new_gcm
    
    type, public :: gcm
        private
        character(len=MAXFILELENGTH) :: name
      contains
        procedure, public  :: init        => init
        procedure, public  :: as_string   => as_string
        procedure, public  :: set_name    => set_name
    end type gcm

contains
    function new_gcm()
        implicit none
        type(gcm) :: new_gcm
        
        call new_gcm%set_name( "Test name" )
        
    end function new_gcm  
    
    subroutine init(this, options)
        implicit none
        class(gcm) :: this
        class(options_type), intent(in) :: options
        
        call this%set_name( options%gcm%name )
        
    end subroutine init

    function as_string(this) result(output)
        implicit none
        class(gcm), intent(in) :: this
        character(len=MAXVARLENGTH) :: output
        
        output = this%name
        
    end function as_string

    subroutine set_name(this, name)
        implicit none
        class(gcm), intent(inout) :: this
        character(len=*) :: name
        
        this%name = name
        
    end subroutine set_name
    
end module gcm_mod
