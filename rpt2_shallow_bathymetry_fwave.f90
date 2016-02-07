subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)

    use geoclaw_module, only: g => grav, dry_tolerance
    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad

    implicit none

    ! Input
    integer, intent(in) :: ixy, maxm, meqn, maux, mwaves, mbc, mx, imp

    real(kind=8), intent(in) :: ql(meqn, 1 - mbc:maxm + mbc)
    real(kind=8), intent(in) :: qr(meqn, 1 - mbc:maxm + mbc)
    real(kind=8), intent(in) :: aux1(maux, 1 - mbc:maxm + mbc)
    real(kind=8), intent(in) :: aux2(maux, 1 - mbc:maxm + mbc)
    real(kind=8), intent(in) :: aux3(maux, 1 - mbc:maxm + mbc)
    real(kind=8), intent(in) :: asdq(meqn, 1 - mbc:maxm + mbc)

    ! Output
    real(kind=8), intent(in out) :: bmasdq(meqn, 1 - mbc:maxm + mbc)
    real(kind=8), intent(in out) :: bpasdq(meqn, 1 - mbc:maxm + mbc)

    ! Locals
    integer :: normal_dir, transverse_dir, i, mw
    real(kind=8) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r, u_l, u_r, v_l, v_r
    real(kind=8) :: b_l, b_r, b_up, b_down, h_hat, u_hat, v_hat, eta

    real(kind=8) :: delta(3), s(3), beta(3), eig_vectors(3, 3)
    real(kind=8) :: dxdcm, dxdcp

    ! Normal vs. Tangential directions
    if (ixy == 1) then
        normal_dir = 2
        transverse_dir = 3
    else
        normal_dir = 3
        transverse_dir = 2
    endif

    ! Initialize output variables
    bmasdq = 0.d0
    bpasdq = 0.d0

    ! Main loop
    do i = 2 - mbc, mx + mbc

        ! Skip this cell interface if both cells are dry
        if (qr(1, i - 1) < dry_tolerance .and. ql(1, i) < dry_tolerance) then
            cycle
        end if

        ! Local states
        h_l = qr(1, i - 1)
        hu_l = qr(normal_dir, i - 1)
        hv_l = qr(transverse_dir, i - 1)
        b_l = aux2(1, i - 1)

        h_r = ql(1, i)
        hu_r = ql(normal_dir, i)
        hv_r = ql(transverse_dir, i)
        b_r = aux2(1, i)

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

        ! Check which side we are splitting for and compute sea-surface and
        ! bathy in above and below cells
        if (imp == 1) then
            ! Left-split
            eta = h_l + b_l
            b_up = aux1(1, i - 1)
            b_down = aux3(1, i - 1)
        else
            ! Right-split
            eta = h_r + b_r
            b_up = aux1(1, i)
            b_down = aux3(1, i)
        end if


        ! Skip this solve if eta is not greater than both bathy in the up
        ! and down directions
        if (eta < max(b_up, b_down)) then
            cycle
        end if

        ! Handle coordinate mapping
        if (coordinate_system == 2) then
            if (ixy == 2) then
                dxdcp = (earth_radius * deg2rad)
                dxdcm = dxdcp
            else
                if (imp == 1) then
                    dxdcp = earth_radius * cos(aux3(3, i-1)) * deg2rad
                    dxdcm = earth_radius * cos(aux1(3, i-1)) * deg2rad
                else
                    dxdcp = earth_radius * cos(aux3(3, i)) * deg2rad
                    dxdcm = earth_radius * cos(aux1(3, i)) * deg2rad
                end if
            end if
        else
            dxdcp = 1.d0
            dxdcm = 1.d0
        end if

        ! Compute relevant Riemann quantities
        s = 0.d0
        beta = 0.d0

        ! Compute Roe-averaged
        u_hat = (u_r * sqrt(h_r) + u_l * sqrt(h_l)) / (sqrt(h_r) + sqrt(h_l))
        v_hat = (v_r * sqrt(h_r) + v_l * sqrt(h_l)) / (sqrt(h_r) + sqrt(h_l))
        h_hat = 0.5d0 * (h_r + h_l)

        ! Compute speeds based on Roe-averages and Einfeldt
        s(1) = min(v_hat - sqrt(g * h_hat), v_l - sqrt(g * h_l))
        s(3) = max(v_hat + sqrt(g * h_hat), v_r + sqrt(g * h_r))
        s(2) = 0.5d0 * (s(1) + s(3))

        ! Compute fluctuations
        delta(1) = asdq(1, i)
        delta(2) = asdq(normal_dir, i)
        delta(3) = asdq(transverse_dir, i)

        beta(1) = s(3) * delta(1) / (s(3) - s(1)) - delta(3) / (s(3) - s(1))
        beta(2) = -s(2) * delta(1) + delta(2)
        beta(3) = delta(3) / (s(3) - s(1)) - s(1) * delta(1) / (s(3) - s(1))

        eig_vectors(:, 1) = [1.d0, s(2), s(1)]
        eig_vectors(:, 2) = [0.d0, 1.d0, 0.d0]
        eig_vectors(:, 3) = [1.d0, s(2), s(3)]

        ! Compute splitting
        do mw = 1, mwaves
            if (s(mw) < 0.d0) then
                bmasdq(1, i) = bmasdq(1, i)     &
                                 + dxdcm * s(mw) * beta(mw) * eig_vectors(1, mw)
                bmasdq(normal_dir, i) = bmasdq(normal_dir, i)     &
                                 + dxdcm * s(mw) * beta(mw) * eig_vectors(1, mw)
                bmasdq(transverse_dir, i) = bmasdq(transverse_dir, i)     &
                                 + dxdcm * s(mw) * beta(mw) * eig_vectors(1, mw)
            else if (s(mw) > 0.d0) then
                bpasdq(1, i) = bpasdq(1, i)    &
                                 + dxdcm * s(mw) * beta(mw) * eig_vectors(1, mw)
                bpasdq(normal_dir, i) = bpasdq(normal_dir, i)    &
                                 + dxdcm * s(mw) * beta(mw) * eig_vectors(1, mw)
                bpasdq(transverse_dir, i) = bpasdq(transverse_dir, i)    &
                                 + dxdcm * s(mw) * beta(mw) * eig_vectors(1, mw)
            end if
        end do

    end do
    ! Main cell interface loop

end subroutine rpt2