module basic_stats_mod

    implicit none
    
contains
    
    !>-------------------------------------
    !! Compute the mean over the first dimension
    !!
    !! precondition: nt>0
    !!
    !!-------------------------------------
    subroutine time_mean(input, mean)
        implicit none
        real, intent(in),   dimension(:,:,:) :: input
        real, intent(inout),dimension(:,:)   :: mean
        
        integer :: nx, ny, nt
        integer :: x, y
        
        nt = size(input,1)
        nx = size(input,2)
        ny = size(input,3)
        
        do y=1,ny
            do x=1,nx
                mean(x,y) = sum(input(:,x,y)) / nt
            end do
        end do
        
    end subroutine time_mean
    
    !>-------------------------------------
    !! Compute the standard deviation over the first dimension
    !!
    !! precondition: nt>1
    !!
    !!-------------------------------------
    subroutine time_stddev(input, stddev, mean_in)
        implicit none
        real, intent(in),    dimension(:,:,:) :: input
        real, intent(inout), dimension(:,:)   :: stddev
        real, intent(in),    dimension(:,:),  optional :: mean_in
        real, dimension(:,:), allocatable :: mean
        
        integer :: nx, ny, nt
        integer :: x, y
        
        nt = size(input,1)
        nx = size(input,2)
        ny = size(input,3)
        
        ! There are one pass stddev algorithms, this is clear and provides a nice option to use an existing mean
        allocate(mean(nx,ny))
        if (.not.present(mean_in)) then
            call time_mean(input,mean)
        else
            mean = mean_in
        endif
        
        ! note, for nt==1, stddev is undefined.  Here return the 10x the mean just to show it is BIG
        if (nt<=1) then
            stddev = abs(mean) * 10
            write(*,*) "WARNING: attempting to compute a standard deviation of one value!", shape(input)
            return
        endif
        
        do y=1,ny
            do x=1,nx
                stddev(x,y) = sqrt( sum( (input(:,x,y) - mean(x,y))**2 ) / (nt-1) )
            end do
        end do
        
    end subroutine time_stddev
    
    function stddev(input, mean_in)
        implicit none
        real, dimension(:) :: input
        real, intent(in), optional :: mean_in
        real :: stddev
        
        real :: mean
        integer :: x, n
        
        n = size(input)
        
        if (present(mean_in)) then
            mean = mean_in
        else
            mean = sum(input) / n
        endif
        
        stddev = sqrt( sum( (input - mean)**2 ) / (n-1) )
        
    end function stddev


end module basic_stats_mod
