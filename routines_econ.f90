!=======================================================================
! CONTENT: Routines for computational economics 
! AUTHOR: Manuel V. Montesinos
! DATE: December 2nd, 2019
! REFERENCES: Fehr, H., and Kindermann, F. (2018): "Introduction to 
! Computational Economics using Fortran", Oxford University Press
!-----------------------------------------------------------------------
! List of subroutines included:
! "error": Throws error message and stops program.
! "init_random_seed": To ensure that each random draw is a new sequence
! "simulate_normal_1": Simulates one draw from a normal distribution 
!                      using Box-Muller transformation.
! "simulate_normal_n": Simulates draws from a normal distribution using 
!                      Box-Muller transformation.
! "simulate_uniform_1": Simulates one draw from a uniform distribution.
! "simulate_uniform_n": Simulates draws from a uniform distribution.
!=======================================================================

module routines_econ
    
    implicit none
    
    !===================================================================
    ! Variable declaration
    
    ! declare everything as public by default (originally, it was set to
    ! private)
    public
    
    ! should the random tbox_seed be set
    logical, private :: tbox_seed = .true.
    
    !-------------------------------------------------------------------
    ! Define public access points
    ! simulation
    public :: simulate_uniform
    public :: simulate_normal
    
    !===================================================================
    ! Interface declarations
    !
    ! INTERFACE simulate_uniform
    !
    ! Simulates uniformly distributed random variables.
    
    interface simulate_uniform

        module procedure simulate_uniform_1, simulate_uniform_n

    end interface
    !--------------------------------------------------------------------
    ! INTERFACE simulate_normal
    !
    ! Simulates normallly distributed random variables.
    
    interface simulate_normal

        module procedure simulate_normal_1, simulate_normal_n

    end interface
    
    !-------------------------------------------------------------------
    ! Separate variable declarations from subroutine and functions
    contains
        
        !===============================================================
        ! SUBROUTINE error
        !
        ! Throws error message and stops program.

        subroutine error(routine, message)
    
            !-----------------------------------------------------------
            ! Input/Output variables
        
            ! routine in which error occured
            character(len=*), intent(in) :: routine
        
            ! error message
            character(len=*), intent(in) :: message
        
            !-----------------------------------------------------------
            ! Routine code
        
            ! write error message
            write(*,'(/a,a,a,a/)')'ERROR ',routine,': ',message
        
            ! stop program
            stop
    
        end subroutine error
        
        !===============================================================
        ! SUBROUTINE init_random_seed
        !
        ! To ensure that each random draw is a new sequence

        subroutine init_random_seed(fixed)
     
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, dt(8)
            integer, parameter :: int64 = selected_int_kind(16)
            integer(kind=int64) :: t
            logical, optional :: fixed
     
            call random_seed(size = n)
            allocate(seed(n))
     
            call system_clock(t)
            if (t == 0) then
                call date_and_time(values=dt)
                t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                    + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                    + dt(3) * 24_int64 * 60 * 60 * 1000 &
                    + dt(5) * 60 * 60 * 1000 &
                    + dt(6) * 60 * 1000 + dt(7) * 1000 &
                    + dt(8)
            endif
        
            if(present(fixed))then
                if(fixed)t = 0
            endif
            do i = 1, n
                seed(i) = lcg(t)
            enddo
        call random_seed(put=seed)
     
        contains
        
            !-----------------------------------------------------------
            ! This simple PRNG might not be good enough for real work, 
            ! but is sufficient for seeding a better PRNG.
            function lcg(s)
     
                implicit none
                integer :: lcg
                integer(int64) :: s
     
                if (s == 0) then
                    s = 104729
                else
                    s = mod(s, 4294967296_int64)
                endif
                s = mod(s * 279470273_int64, 4294967291_int64)
                lcg = int(mod(s, int(huge(0), int64)), kind(0))
     
            end function lcg
     
        end subroutine init_random_seed
        
        !===============================================================
        ! SUBROUTINE simulate_uniform_1
        !
        ! Simulates one draw from a uniform distribution.
        
        subroutine simulate_uniform_1(x, a, b, fixed)
     
        implicit none
     
            !-----------------------------------------------------------     
            ! Input/Output variables 
     
            ! point into which the draw should be saved
            real*8, intent(out) :: x
     
            ! left end of the distribution
            real*8, optional :: a
     
            ! right end of the distribution
            real*8, optional :: b
        
            ! should the random seed be initialized at a fixed values
            logical, optional :: fixed
     
            !---------------------------------------------------------------
            ! Other variables 
     
            real*8 :: a_c, b_c
     
            !---------------------------------------------------------------
            ! Routine code
     
            ! initialize parameters
            a_c = 0d0
            if(present(a))a_c = a
            b_c = 1d0
            if(present(b))b_c = b
     
            ! initialize the random seed
            if(tbox_seed)then
                if(present(fixed))then
                    call init_random_seed(fixed)
                else
                    call init_random_seed()
                endif
                tbox_seed = .false.
            endif
     
            ! draw the random number
            call random_number(x)
     
            x = a_c + (b_c-a_c)*x
     
        end subroutine simulate_uniform_1
        
        !===============================================================
        ! SUBROUTINE simulate_uniform_n
        !
        ! Simulates a series draw from a uniform distribution.

        subroutine simulate_uniform_n(x, a, b, fixed)
     
            implicit none
        
            !-----------------------------------------------------------
            ! Input/Output variables
     
            ! point into which the draw should be saved
            real*8, intent(out) :: x(:)
     
            ! left end of the distribution
            real*8, optional :: a
     
            ! right end of the distribution
            real*8, optional :: b
        
            ! should the random seed be initialized at a fixed values
            logical, optional :: fixed
     
            ! Other variables
     
            real*8 :: a_c, b_c
            
            !-----------------------------------------------------------
            ! Routine Code
     
            ! initialize parameters
            a_c = 0d0
            if(present(a))a_c = a
            b_c = 1d0
            if(present(b))b_c = b
     
            ! initialize the random seed
            if(tbox_seed)then
                if(present(fixed))then
                    call init_random_seed(fixed)
                else
                    call init_random_seed()
                endif
                tbox_seed = .false.
            endif
     
            call random_number(x)
     
            x = a_c + (b_c-a_c)*x
     
        end subroutine simulate_uniform_n
        
        !===============================================================
        ! SUBROUTINE simulate_normal_1
        !
        ! Simulates one draw from a normal distribution using
        ! Box-Muller tranformation.
        
        subroutine simulate_normal_1(x, mu, sigma, fixed)
     
        implicit none
            
            !-----------------------------------------------------------
            ! Input/Output variables
     
            ! point into which the draw should be saved
            real*8, intent(out) :: x
     
            ! expectation of distribution
            real*8, optional :: mu
     
            ! variance of distribution
            real*8, optional :: sigma
        
            ! should the random seed be initialized at a fixed values
            logical, optional :: fixed
            
            !-----------------------------------------------------------
            ! Other variables
            real*8 :: uni1, uni2, mu_c, sigma_c
            real*8 :: pi = 3.141592653589793d0
            
            !-----------------------------------------------------------
            ! Routine code
     
            ! initialize parameters
            mu_c = 0d0
            if(present(mu))mu_c = mu
            sigma_c = 1d0
            if(present(sigma))sigma_c = sigma
        
            if(sigma_c <= 0d0)then
                call error('simulate_normal','sigma has zero or negative value')
            endif
     
            ! simulate a uniform variable draw
            if(present(fixed))then
                call simulate_uniform_1(uni1, fixed=fixed)
                call simulate_uniform_1(uni2, fixed=fixed)
            else
                call simulate_uniform_1(uni1)
                call simulate_uniform_1(uni2)
            endif
     
            ! transform by Box-Muller transformation
            x = sqrt(-2d0*log(uni1))*cos(2d0*pi*uni2)
     
            ! transform to mean and variance
            x = mu_c + sqrt(sigma_c)*x
            
        end subroutine simulate_normal_1
        
        !===============================================================
        ! SUBROUTINE simulate_normal_n
        !
        ! Simulates draws from a normal distribution using
        ! Box-Muller tranformation.
        
        subroutine simulate_normal_n(x, mu, sigma, fixed)
     
            implicit none     
        
            !-----------------------------------------------------------
            ! Input/Output variables
     
            ! point into which the draw should be saved
            real*8, intent(out) :: x(:)
     
            ! expectation of distribution
            real*8, optional :: mu
     
            ! variance of distribution
            real*8, optional :: sigma
        
            ! should the random seed be initialized at a fixed values
            logical, optional :: fixed
            
            !-----------------------------------------------------------
            ! Other variables
     
            real*8 :: uni1(size(x, 1)/2), uni2(size(x, 1)/2), mu_c, sigma_c
            integer :: n, in, n2
            real*8 :: pi = 3.141592653589793d0
            
            !-----------------------------------------------------------
            ! Routine code
     
            ! initialize parameters
            mu_c = 0d0
            if(present(mu))mu_c = mu
            sigma_c = 1d0
            if(present(sigma))sigma_c = sigma
        
            if(sigma_c <= 0d0)then
                call error('simulate_normal','sigma has zero or negative value')
            endif
     
            ! get size of x
            n = size(x, 1)
            n2 = n/2
     
            ! simulate a uniform variable draw
            if(present(fixed))then
                call simulate_uniform_n(uni1, fixed=fixed)
                call simulate_uniform_n(uni2, fixed=fixed)
            else
                call simulate_uniform_n(uni1)
                call simulate_uniform_n(uni2)
            endif
     
            ! transform by Box-Muller transformation
            do in = 1, n2
                x(2*(in-1)+1) = sqrt(-2d0*log(uni1(in)))*cos(2d0*pi*uni2(in))
                x(2*(in-1)+2) = sqrt(-2d0*log(uni1(in)))*sin(2d0*pi*uni2(in))
            enddo
     
            if(2*n2 /= n)then
                call simulate_normal_1(x(n), 0d0, 1d0)
            endif
     
            ! transform to mean and variance
            x = mu_c + sqrt(sigma_c)*x
     
        end subroutine simulate_normal_n
    !===================================================================
    
end module
