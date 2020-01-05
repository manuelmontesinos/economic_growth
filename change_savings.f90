!=======================================================================
! CONTENT: Solow-Swan Model (stochastic version). Short- and long-run
!          effects of a permanent increase in the value of the savings
!          rate        
! AUTHOR: Manuel V. Montesinos
! DATE: January 5th, 2020
! REFERENCE: Fernandez, Novales and Ruiz (2009), "Economic Growth:
!            Theory and Numerical Solutions", Springer 
!=======================================================================

program change_savings
    
    !-------------------------------------------------------------------
    ! Set the structural parameters
    !
    ! alpha -> output elasticity of capital
    ! n -> population growth
    ! s0 -> initial savings rate
    ! s1 -> new savings rate
    ! B -> aggregate productivity parameter
    ! delta -> depreciation rate
    ! sigma -> inverse of the elasticity of substitution
    ! beta -> discount factor
    ! gammap -> growth rate of technology
    ! T -> transition paths
    
    real*8, parameter :: alpha=0.36, n=0.01, delta=0.075, beta=0.95, &
        gammap=0, s0=0.2, s1=0.3 
        
    integer, parameter :: B=3, sigma=3, T=223
    
    ! Declare the variables used below
    real*8 :: kss0, yss0, css0, iss0, uss0, kss1, yss1, css1, iss1, &
        uss1, w, ka, kt, yt, ct, it, st, ut, wl, D, wq, F0, F1, kss0c, &
        yss0c, css0c, iss0c, uss0c, kss1c, yss1c, css1c, iss1c, uss1c, &
        k0a, k0n, gyt, gylt, gyqt, gyct
    
    real*8, dimension(T+1) :: k, c, y, i, u, s, lnYL, kl, cl, yl, il, &
        ul, lnYLl, kq, cq, yq, iq, uq, lnYLlq, kc, cc, yc, ic, uc, lnYLlc
    
    real*8, dimension(T) :: gy, gyl, gyq, gyc
    
    integer, dimension(T+1) :: trend
    
    integer :: tt
    
    !-------------------------------------------------------------------
    
    ! Steady state (variables per unit of effective labor) under the 
    ! initial savings rate
    
    ! stock of capital
    kss0 = (s0*B/((1+n)*(1+gammap)-(1-delta)))**(1/(1-alpha))
    
    ! output
    yss0 = B*(kss0**alpha)
    
    ! consumption
    css0 = (1-s0)*yss0
    
    ! investment
    iss0 = s0*yss0
    
    ! utility
    if (sigma == 1) then
        uss0 = log(css0)
    else
        uss0 = ((css0**(1-sigma))-1)/(1-sigma)
    endif
    
    !-------------------------------------------------------------------
    
    ! Steady state under the new savings rate
    
    ! stock of capital
    kss1 = (s1*B/((1+n)*(1+gammap)-(1-delta)))**(1/(1-alpha))
    
    ! output
    yss1 = B*(kss1**alpha)
    
    ! consumption
    css1 = (1-s1)*yss1
    
    ! investment
    iss1 = s1*yss1
    
    ! utility
    if (sigma == 1)then
        uss1 = log(css1)
    else
        uss1 = ((css1**(1-sigma))-1)/(1-sigma)
    endif
    
    !-------------------------------------------------------------------
    
    ! Exact solution
    
    k(1) = kss0
    c(1) = css0
    y(1) = yss0
    i(1) = iss0
    u(1) = uss0
    s(1) = s0
    w = 0
    ka = kss0
    
    do tt = 1, T
        if (tt < 11)then
            kt = (s0/((1+n)*(1+gammap)))*B*(ka**alpha)+&
                ((1-delta)/((1+n)*(1+gammap)))*ka
            ka = kt
            k(tt+1) = kt
            yt = B*kt**alpha
            y(tt+1) = yt
            ct = (1-s0)*yt
            c(tt+1) = ct
            it = s0*yt
            i(tt+1) = it
            st = s0
            s(tt+1) = st
            if (sigma == 1)then
                ut = log(ct)
            else
                ut = ((ct**(1-sigma))-1)/(1-sigma)
            endif    
            u(tt+1) = ut
            w = w+(beta**tt)*ut
        else
            kt = (s1/((1+n)*(1+gammap)))*B*(ka**alpha)+&
                ((1-delta)/((1+n)*(1+gammap)))*ka
            ka = kt
            k(tt+1) = kt
            yt = B*kt**alpha
            y(tt+1) = yt
            ct = (1-s1)*yt
            c(tt+1) = ct
            it = s1*yt
            i(tt+1) = it
            st = s1
            s(tt+1) = s1
            if (sigma == 1)then
                ut = log(ct)
            else
                ut = ((ct**(1-sigma))-1)/(1-sigma)
            endif
            u(tt+1) = ut
            w = w+(beta**tt)*ut
        endif
    enddo
    
    ! log output per capita
    do tt = 1, T+1
        trend(tt) = tt
        lnYL(tt) = log(y(tt))+log(1+gammap)*tt
    enddo
     
    !-------------------------------------------------------------------
    
    ! Linear approximation
    
    kl(1) = kss0
    cl(1) = css0
    yl(1) = yss0
    il(1) = iss0
    ul(1) = uss0
    wl = 0
    ka = kss0
    D = alpha+(1-alpha)*(1-delta)/((1+n)*(1+gammap))
    
    do tt = 1, T
        if (tt < 11)then
            kt = kss0+D*(ka-kss0)
            ka = kt
            kl(tt+1) = kt
            yt = B*kt**alpha
            yl(tt+1) = yt
            ct = (1-s0)*yt
            cl(tt+1) = ct
            it = s0*yt
            il(tt+1) = it
            if (sigma==1)then
                ut = log(ct)
            else
                ut = ((ct**(1-sigma))-1)/(1-sigma)
            endif
            ul(tt+1) = ut
            wl = wl+(beta**tt)*ut
            else
            kt = kss1+D*(ka-kss1)
            ka = kt
            kl(tt+1) = kt
            yt = B*kt**alpha
            yl(tt+1) = yt
            ct = (1-s1)*yt
            cl(tt+1) = ct
            it = s1*yt
            il(tt+1) = it
            if (sigma==1)then
                ut = log(ct)
            else
                ut = ((ct**(1-sigma))-1)/(1-sigma)
            endif
            ul(tt+1) = ut
            wl = wl*(beta**tt)*ut
        endif
    enddo
    
    ! log output per capita
    do tt = 1, T+1
        lnYLl(tt) = log(yl(tt))+log(1+gammap)*tt
    enddo
    
    !-------------------------------------------------------------------
    
    ! Quadratic approximation
    kq(1) = kss0
    cq(1) = css0
    yq(1) = yss0
    iq(1) = iss0
    uq(1) = uss0
    wq = 0
    ka = kss0
    D = alpha+(1-alpha)*(1-alpha)/((1+n)*(1+gammap))
    F0 = -0.5*alpha*(1-alpha)*((1+n)*(1+gammap)-(1-delta))/((1+n)*(1+gammap)*kss0)
    F1 = -0.5*alpha*(1-alpha)*((1+n)*(1+gammap)-(1-delta))/((1+n)*(1+gammap)*kss1)
    
    do tt = 1, T
        if (tt < 11)then
            kt = kss0+D*(ka-kss0)+F0*(ka-kss0)**2
            ka = kt
            kq(tt+1) = kt
            yt = B*kt**alpha
            yq(tt+1) = yt
            ct = (1-s0)*yt
            cq(tt+1) = ct
            it = s0*yt
            iq(tt+1) = it
            if (sigma==1)then
                ut = log(ct)
            else
                ut = ((ct**(1-sigma))-1)/(1-sigma)
            endif
            uq(tt+1) = ut
            wq = wq+(beta**tt)*ut
        else
            kt = kss1+D*(ka-kss1)+F1*(ka-kss1)**2
            ka = kt
            kq(tt+1) = kt
            yt = B*kt**alpha
            yq(tt+1) = yt
            ct = (1-s1)*yt
            cq(tt+1) = ct
            it = s1*yt
            iq(tt+1) = it
            if (sigma==1)then
                ut = log(ct)
            else
                ut = ((ct**(1-sigma))-1)/(1-sigma)
            endif
            uq(tt+1) = ut
            wq = wq+(beta**tt)*ut
        endif
    enddo
        
    ! log output per capita
    do tt = 1, T+1
        lnYLlq(tt) = log(yq(tt))+log(1+gammap)*tt
    enddo
    
    !-------------------------------------------------------------------
    
    ! Exact continuous
    
    kss0c = (B*s0/(n+delta+gammap))**(1/(1-alpha))
    yss0c = B*(kss0c**alpha)
    css0c = (1-s0)*yss0c
    iss0c = s0*yss0c
    
    if (sigma==1)then
        uss0c = log(css0c)
    else
        uss0c = ((css0c**(1-sigma))-1)/(1-sigma)
    endif
    
    kss1c = (B*s1/(n+delta+gammap))**(1/(1-alpha))
    yss1c = B*(kss1c**alpha)
    css1c = (1-s1)*yss1c
    iss1c = s1*yss1c
    
    if (sigma==1)then
        uss1c = log(css1c)
    else
        uss1c = ((css1c**(1-sigma))-1)/(1-sigma)
    endif
    
    kc(1) = kss0
    cc(1) = css0
    yc(1) = yss0
    ic(1) = iss0
    uc(1) = uss0
    k0a = kss0c
    
    do tt = 1, T
        if (tt < 11)then
            kt = ((k0a**(1-alpha)-(s0*B/(n+delta+gammap)))*&
                exp((alpha-1)*(n+delta+gammap)*t)+s0*B/(n+delta+gammap))**&
                (1/(1-alpha))
            kc(tt+1) = kt
            yt = B*kt**alpha
            yc(tt+1) = yt
            ct = (1-s0)*yt
            cc(tt+1) = ct
            it = s0*yt
            ic(tt+1) = it
            k0n = kt
        else
            kt = ((k0n**(1-alpha)-(s1*B/(n+delta+gammap)))*&
            exp((alpha-1)*(n+delta+gammap)*(tt-10))+s1*B/(n+delta+gammap))**&
            (1/(1-alpha))
            kc(tt+1) = kt
            yt = B*kt**alpha
            yc(tt+1) = yt
            ct = (1-s1)*yt
            cc(tt+1) = ct
            it = s1*yt
            ic(tt+1) = it
        endif    
    enddo
    
    ! log output per capita
    do tt = 1, T+1
        lnYLlc(tt) = log(yc(tt))+log(1+gammap)*tt
    enddo
    
    do tt = 1, T
        gyt = ((y(tt+1)-y(tt))/y(tt))*100
        gy(tt) = gyt
        gylt = ((yl(tt+1)-yl(tt))/yl(tt))*100
        gyl(tt) = gylt
        gyqt = ((yq(tt+1)-yq(tt))/yq(tt))*100
        gyq(tt) = gyqt
        gyct = ((yc(tt+1)-yc(tt))/yc(tt))*100
        gyc(tt) = gyct
    enddo
    
    !===================================================================
    ! Save the output in a txt file
    10 format(i3.1, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, & 
    f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, & 
    f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x, &
    f10.4, 5x, f10.4, 5x, f10.4, 5x, f10.4, 5x)
    
    open(unit=1,file="change_savings_output.txt")
    do tt=1,100
        write(unit=1, fmt=10) tt, k(tt), kl(tt), kq(tt), kc(tt), &
        c(tt), cl(tt), cq(tt), cc(tt), y(tt), yl(tt), yq(tt), yc(tt), &
        gy(tt), gyl(tt), gyq(tt), gyc(tt), lnYL(tt), lnYLl(tt), &
        lnYLlq(tt), lnYLlc(tt)
    end do
    close(unit=1)
        
    !===================================================================
    
end
