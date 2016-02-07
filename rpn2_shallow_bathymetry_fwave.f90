subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

    use amr_module, only: mcapa

    use geoclaw_module, only: g => grav, dry_tolerance
    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad

    implicit none 

    ! Input
    integer, intent(in) :: ixy, maxm, mx, meqn, mwaves, mbc, maux
    real(kind=8), intent(in) :: ql(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: qr(meqn, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxl(maux, 1-mbc:maxm+mbc)
    real(kind=8), intent(in) :: auxr(maux, 1-mbc:maxm+mbc)

    ! Output
    real(kind=8), intent(in out) :: s(mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: fwave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: amdq(meqn,1-mbc:maxm+mbc)
    real(kind=8), intent(in out) :: apdq(meqn,1-mbc:maxm+mbc)
    
    ! Locals
    integer :: normal_dir, tangential_dir, i, mw
    real(kind=8) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r, u_l, u_r, v_l, v_r
    real(kind=8) :: b_l, b_r, phi_l, phi_r

    real(kind=8) :: h_hat, u_hat, v_hat, c_hat

    real(kind=8) :: delta(3), beta(3)
    real(kind=8) :: dxdc

    ! Normal vs. Tangential directions
    if (ixy == 1) then
        normal_dir = 2
        tangential_dir = 3
    else
        normal_dir = 3
        tangential_dir = 2
    endif


    ! Main loop
    do i = 2 - mbc, mx + mbc

        ! Skip this cell interface if both cells are dry
        if (qr(1, i - 1) < dry_tolerance .and. ql(1, i) < dry_tolerance) then
            cycle
        end if

        ! Extract states
        h_l = qr(1, i - 1)
        h_r = ql(1, i)
        hu_l = qr(normal_dir, i - 1)
        hu_r = ql(normal_dir, i)
        hv_l = qr(tangential_dir, i - 1)
        hv_r = ql(tangential_dir, i)
        b_l = auxr(1, i - 1)
        b_r = auxl(1, i)

        ! Determine velocities
        if (h_l < dry_tolerance) then
            h_l = 0.d0
            u_l = 0.d0
            v_l = 0.d0
        else
            u_l = hu_l / h_l
            v_l = hv_l / h_l
        end if
        if (h_r < dry_tolerance) then
            h_r = 0.d0
            u_r = 0.d0
            v_r = 0.d0
        else
            u_r = hu_r / h_r
            v_r = hv_r / h_r
        end if

        ! Momentum flux
        phi_l = h_l * u_l**2 + 0.5d0 * g * h_l**2
        phi_r = h_r * u_r**2 + 0.5d0 * g * h_r**2

        ! Compute Roe averages
        h_hat = 0.5d0 * (h_r + h_l)
        u_hat = (sqrt(h_r) * u_r + sqrt(h_l) * u_l) / (sqrt(h_r) + sqrt(h_l))
        v_hat = (sqrt(h_r) * v_r + sqrt(h_l) * v_l) / (sqrt(h_r) + sqrt(h_l))
        c_hat = sqrt(g * h_hat)
    
        ! Flux differences
        delta(1) = hu_r - hu_l
        delta(normal_dir) = phi_r - phi_l + 0.5d0 * g * (h_r + h_l) * (b_r - b_l)
        delta(tangential_dir) = h_r * u_r * v_r - h_l * u_l * v_l
   
        ! Wave speeds - Einfeldt speeds
        s(1, i) = min(u_hat - c_hat, u_l - sqrt(g * h_l))
        s(3, i) = max(u_hat + c_hat, u_r + sqrt(g * h_r))
        s(2, i) = 0.5d0 * (s(1, i) + s(3, i))

        ! Wave strengths
        beta(1) = (s(3, i) * delta(1) - delta(2)) / (s(3, i) - s(1, i))
        beta(3) = (delta(2) - s(1, i) * delta(1)) / (s(3, i) - s(1, i))
        beta(2) = delta(3) - v_l * beta(1) - v_r * beta(3)

        ! F-Waves
        fwave(:, 1, i) = [1.d0, s(1, i), v_l] * beta(1)
        fwave(:, 2, i) = [0.d0, 0.d0, 1.d0] * beta(2)
        fwave(:, 3, i) = [1.d0, s(3, i), v_r] * beta(3)

    end do
    ! End of main loop

    ! Capacity mapping for lat-long grids
    if (mcapa > 0) then
        do i = 2 - mbc, mx + mbc
            if(ixy == 1) then
                dxdc = earth_radius * deg2rad
            else
                dxdc = earth_radius * cos(auxl(3, i)) * deg2rad
            end if

            s(:, i) = dxdc * s(:, i)
            fwave(:, :, i) = dxdc * fwave(:, :, i)
        end do
    end if

    ! Accumulate fluctuations
    do i = 2 - mbc, mx + mbc
        do mw = 1, mwaves
            if (s(mw, i) < 0.d0) then
                amdq(:, i) = amdq(:, i) + fwave(:, mw, i)
            else if (s(mw, i) > 0.d0) then
                apdq(:, i) = apdq(:, i) + fwave(:, mw, i)
            else
                amdq(:, i) = amdq(:, i) + 0.5d0 * fwave(:, mw, i)
                apdq(:, i) = apdq(:, i) + 0.5d0 * fwave(:, mw, i)
            end if
        end do
    end do
    
    print *, s
    stop
    
end subroutine rpn2
