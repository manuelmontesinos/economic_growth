!=======================================================================
! CONTENT: Solow-Swan Model (stochastic version)
! AUTHOR: Manuel V. Montesinos
! DATE: October 2nd, 2019
! REFERENCE: Fernandez, Novales and Ruiz (2009), "Economic Growth:
! Theory and Numerical Solutions", Springer 
!=======================================================================

include "routines_econ.f90"

program solow_stochastic
    
    use routines_econ
        
    implicit none

    ! STRUCTURAL PARAMETERS:
    ! alpha -> output elasticity of capital
    ! n -> population growth
    ! s -> savings rate
    ! B -> aggregate productivity parameter
    ! delta -> depreciation rate
    ! sigma -> inverse of the elasticity of substitution 
    ! beta -> time discount factor
    ! gammaA -> technology growth
    ! rho -> autoregressive parameter in productivity shock
    ! sigmae -> std dev of the innovation for technology shock
    
    real*8, parameter :: alpha = 0.36, s = 0.36, delta = 0.075, &
    beta = 0.9, rho = 0.9, sigmae = 0.01
    integer, parameter :: n = 0, B = 5, sigma = 3, T = 260
    real*8, parameter :: gammaA = 0
    
    ! DECLARE VARIABLES:
    real*8 :: kss, yss, css, iss, pmgss
    real*8 :: k0, theta0, y0, c0, i0, pmg0
    real*8 :: mu, D1, D2
    real*8 :: thetat, kt, yt, ct, it, pmgt, gyt, gylt
    real*8, dimension(T+1) :: theta, k, y, c, i, pmg, kl, yl, cl, il, &
     pmgl, lnYL, lnYLl, errork, errorc, errory
    real*8, dimension(T) :: gy, gyl
    real*8, dimension(T+100) :: innov0
    real*8, dimension(T) :: innov
    integer, dimension(T+1) :: trend
    integer, parameter :: thetass = 1 
    integer :: tt
    
    !-------------------------------------------------------------------
    ! STEADY STATE (variables per unit of effective labor):
    
    ! Stock of physical capital 
    kss = (s*B*thetass/((1+n)*(1+gammaA)-(1-delta)))**(1/(1-alpha))
    ! Output
    yss = B*thetass*(kss**alpha)
    ! Consumption
    css = (1-s)*yss
    ! Investment
    iss = s*yss
    ! Marginal productivity
    pmgss = alpha*B*thetass*kss**(alpha-1)
    
    !-------------------------------------------------------------------
    ! TRANSITION PATHS. EXACT SOLUTION
    
    ! Initial conditions    
    k0 = kss
    theta0 = thetass
    y0 = B*theta0*(k0**alpha)
    c0 = (1-s)*y0
    i0 = s*y0
    pmg0 = alpha*B*theta0*k0**(alpha-1)
        
    k(1) = k0
    c(1) = c0
    y(1) = y0
    i(1) = i0
    pmg(1) = pmg0
    theta(1) = theta0
    
    mu = 0
    call simulate_normal(innov0, mu, sigmae**2, .true.)
    innov = innov0(101:T+100)

    do tt = 1,T
        thetat = exp(rho*log(theta0)+innov(tt))
        theta0 = thetat
        theta(tt+1) = thetat
        kt = (s/((1+n)*(1+gammaA)))*B*thetat*k0**alpha + &
         ((1-delta)/((1+n)*(1+gammaA)))*k0
        k0 = kt
        k(tt+1) = kt
        yt = B*thetat*kt**alpha
        y(tt+1) = yt
        ct = (1-s)*yt
        c(tt+1) = ct
        it = s*yt
        i(tt+1) = it
        pmgt = alpha*B*thetat*kt**(alpha-1)
        pmg(tt+1) = pmgt
    enddo 
    
    do tt = 1,T+1
        trend(tt) = tt
    enddo
    
    ! Log output per capita
    do tt = 1,T+1
        lnYL(tt) = log(y(tt))+log(1+gammaA)*trend(tt)
    enddo
    
    !===================================================================
    ! LINEAR APPROXIMATION
    ! Initial conditions
    k0 = kss
    theta0 = thetass
    
    y0 = B*theta0*k0**alpha
    c0 = (1-s)*y0
    i0 = s*y0
    pmg0 = alpha*B*theta0*k0**(alpha-1)
    
    kl(1) = k0
    cl(1) = c0
    yl(1) = y0
    il(1) = i0
    pmgl(1) = pmg0
    
    D1 = alpha+(1-alpha)*(1-delta)/((1+n)*(1+gammaA))
    D2 = (1-(1-delta)/((1+n)*(1+gammaA)))*(kss/thetass)
    
    do tt = 1,T
        kt = kss+D1*(k0-kss)+D2*(theta(tt)-thetass)
        k0 = kt
        kl(tt+1) = kt
        yt = B*theta(tt+1)*kt**alpha
        yl(tt+1) = yt
        ct = (1-s)*yt
        cl(tt+1) = ct
        it = s*yt
        il(tt+1) = it
        pmgt = alpha*B*theta(tt+1)*kt**(alpha-1)
        pmgl(tt+1) = pmgt
    enddo
    
    ! Log output per capita
    do tt = 1,T+1
        lnYLl(tt) = log(yl(tt))+log(1+gammaA)*trend(tt)
    enddo
    
    do tt = 1,T+1
        gyt = ((y(tt+1)-y(tt))/y(tt))*100
        gy(tt) = gyt
        gylt = ((yl(tt+1)-yl(tt))/yl(tt))*100
        gyl(tt) = gylt
    enddo
    
    errork(:) = k(:)-kl(:)
    errorc(:) = c(:)-cl(:)
    errory(:) = y(:)-yl(:)
        
    !===================================================================
    ! SAVE OUTPUT (exact solution vs linear approximation)
    
    ! Save the output in a txt file
    10 format(i3.1, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, & 
    f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, & 
    f10.4, 5x, f10.4, 5x)
    
    open(unit=1,file="solow_stochastic_output.txt")
    do tt=1,T+1
        write(unit=1, fmt=10) tt, k(tt), kl(tt), c(tt), cl(tt), y(tt), &
        yl(tt), gy(tt), gyl(tt), lnYL(tt), lnYLl(tt), pmg(tt), pmgl(tt)
    end do
    close(unit=1)
    
    !===================================================================
    
end program solow_stochastic
