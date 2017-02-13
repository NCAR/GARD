module random_mod

    use model_constants

    implicit none

    private
    double precision, parameter :: root_pi_over_8 = 0.626657068658 ! sqrt(kPI/8)

    integer, parameter :: size_of_normal_cdf_lookup_table = 1000000
    double precision, allocatable :: fast_cdf_lookup_table(:), cdf_lookup_table(:,:)
    
    interface normal_cdf
        module procedure real_normal_cdf, double_normal_cdf
    end interface
    
    public :: normal_cdf, box_muller_random, get_normal_from_cdf
    
contains
    
    !>------------------------------------------------
    !! Use the Box-Muller Transform to convert uniform to normal random deviates
    !!
    !! Note random_sample should be an allocated 1D real array
    !! On return, random_sample will be filled with random normal (0,1) data
    !!
    !! Caveat: this transform is unable to generate extremely high values (>6.66)
    !! Having switched to double precision it may do better now.
    !!
    !-------------------------------------------------
    subroutine box_muller_random(random_sample)
        implicit none
        real, intent(inout) :: random_sample(:)
        integer :: n,i

        double precision :: u1, u2, s
        double precision, allocatable :: double_random(:)

        n = size(random_sample)
        allocate(double_random(n))
        call random_number(double_random)

        do i=1,n,2
            u1 = double_random(i)
            if (i<n) then
                u2 = double_random(i+1)
            else
                call random_number(u2)
            endif

            s = sqrt(-2 * log(u1))
            random_sample(i) = s * cos(2 * kPI * u2)
            if (i<n) then
                random_sample(i+1) = s * sin(2 * kPI * u2)
            endif

        enddo

    end subroutine box_muller_random
    
    !>------------------------------------------------
    !! Computes an approximation to the normal CDF based on Aludaat and Alodat (2008)
    !!
    !! Aludaat, K.M. and Alodat, M.T. (2008). A note on approximating the normal distribution function. Applied Mathematical Sciences, Vol 2, no 9, pgs 425-429.
    !!
    !!------------------------------------------------
    function real_normal_cdf(normal) result(cdf)
        implicit none
        real :: normal
        real :: cdf
        
        cdf = 0.5 * (1 + sqrt((1 - exp(0 - (root_pi_over_8 * normal**2))) ))
        
        if (normal < 0) then
            cdf = 1 - cdf
        endif
        
    end function real_normal_cdf

    !>------------------------------------------------
    !! Computes an approximation to the normal CDF based on Aludaat and Alodat (2008)
    !!
    !! Aludaat, K.M. and Alodat, M.T. (2008). A note on approximating the normal distribution function. Applied Mathematical Sciences, Vol 2, no 9, pgs 425-429.
    !!
    !!------------------------------------------------
    function double_normal_cdf(normal) result(cdf)
        implicit none
        double precision :: normal
        double precision :: cdf
        
        cdf = 0.5 * (1 + sqrt((1.d0 - exp(0.d0 - (root_pi_over_8 * normal**2))) ))
        
        if (normal < 0) then
            cdf = 1 - cdf
        endif
        
    end function double_normal_cdf

    !>------------------------------------------------
    !! Create a look up table to convert from Uniform random number to Normal random
    !!
    !!------------------------------------------------
    subroutine create_cdf_lookup_table()
        implicit none
        
        integer :: i
        real :: normal_random_value, normal_step
        
        ! in case another thread beat this one to the lock after the if
        if (.not.allocated(fast_cdf_lookup_table)) then
            
            allocate(cdf_lookup_table(size_of_normal_cdf_lookup_table, 2))
            
            normal_random_value = (-5.5)
            normal_step = 11.0 / (size_of_normal_cdf_lookup_table - 1)
            
            do i = 1, size_of_normal_cdf_lookup_table
                cdf_lookup_table(i,1) = normal_random_value
                cdf_lookup_table(i,2) = normal_cdf(normal_random_value)

                normal_random_value = normal_random_value + normal_step
            enddo
            
            ! now invert this slow lookup table to create a fast lookup table
            allocate(fast_cdf_lookup_table(size_of_normal_cdf_lookup_table))
            
            do i = 1,size_of_normal_cdf_lookup_table
                fast_cdf_lookup_table(i) = look_up(dble(i) / size_of_normal_cdf_lookup_table, cdf_lookup_table)
            enddo
            
            deallocate(cdf_lookup_table)
        endif
        
    end subroutine create_cdf_lookup_table

    !>------------------------------------------------
    !! Perform a binary search through LUT to find value, then interpolate between the two closest entries. 
    !!
    !!------------------------------------------------
    function look_up(value, LUT) result(output)
        implicit none
        double precision, intent(in)  :: value
        double precision, intent(in)  :: LUT(:,:)
        real :: output
        
        real :: fraction, denominator
        integer :: i,n
        integer :: step
        
        n = size(LUT, 1)
        i = n/2
        step = n/2
        
        ! perform a binary search
        do while(step/2.0 >= 1)
            step = (step+1) / 2
            if (LUT(i,2) > value) then
                i = max(1, i - step)
            else
                i = min(n, i + step)
            endif
            
        enddo
        
        ! figure out which side of this index we fall on and interpolate
        if (value > LUT(i,2)) then
            if (i==n) then
                output = LUT(i,1)
            else
                denominator = max(1e-20, (LUT(i+1,2) - LUT(i,2)))
                fraction = (value - LUT(i,2)) / denominator
                output = fraction * LUT(i+1,1) + (1-fraction) * LUT(i,1)
            endif
        else
            if (i==1) then
                output = LUT(i,1)
            else
                denominator = max(1e-20, (LUT(i,2) - LUT(i-1,2)))
                fraction = (LUT(i,2) - value) / denominator
                
                output = fraction * LUT(i-1,1) + (1-fraction) * LUT(i,1)
            endif
        endif
    end function look_up
        

    !>------------------------------------------------
    !! Calculate index into precomputed LUT to find cdf_value, then interpolate between the two closest entries. 
    !!
    !!------------------------------------------------
    function get_normal_from_cdf(cdf_value) result(random_normal)
        implicit none
        real, intent(in)  :: cdf_value
        real :: random_normal
        
        integer :: index
        double precision :: real_index
        double precision :: fraction
        
        if (.not.allocated(fast_cdf_lookup_table)) then
            !$omp critical(cdf_lookup)
            call create_cdf_lookup_table()
            !$omp end critical(cdf_lookup)
        endif

        ! in the fast look up table we can just compute the optimal index
        real_index = cdf_value * 1.0d0 * size_of_normal_cdf_lookup_table
        index = min(size_of_normal_cdf_lookup_table, max(1, nint( real_index)) )
        
        ! then interpolate between this index and the next closest
        if (index > real_index) then
            if (index == 1) then
                random_normal = fast_cdf_lookup_table(index)
            else
                fraction = (index - real_index)
                random_normal = (1-fraction) * fast_cdf_lookup_table(index)         &
                    + fraction * fast_cdf_lookup_table(max(1,index-1))
            endif
        else
            if (index == size_of_normal_cdf_lookup_table) then
                random_normal = fast_cdf_lookup_table(index)
            else
                fraction = (real_index - index)
                random_normal = (1-fraction) * fast_cdf_lookup_table(index)         &
                    + fraction * fast_cdf_lookup_table(min(size_of_normal_cdf_lookup_table,index+1))
            endif
        endif
        
    end function get_normal_from_cdf


end module random_mod
