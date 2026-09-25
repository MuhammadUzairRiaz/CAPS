! CAPS k-space electrostatics (Fortran 2018, called from C++ through iso_c_binding).
!
!   caps_f_pme    smooth particle-mesh Ewald, reciprocal part (Essmann, Perera, Berkowitz, Darden, Lee and Pedersen,
!                 J. Chem. Phys. 103, 8577 (1995)): cardinal B-spline charge spreading of order p, 3-D FFT, the
!                 influence function with the B-spline moduli, energy, forces and the virial tensor
!   caps_f_ewald  the plain Ewald reciprocal sum over a box of wave vectors (reference and small systems)
!   caps_f_fft3d  the 3-D complex FFT itself (mixed radix 2, 3, 5 and any other prime), for tests
!
! Units: positions in Å, charges in e, kcoul the Coulomb constant (kcal/mol·Å/e²); energies in kcal/mol, forces in
! kcal/mol/Å. The cell is given by its three edge vectors a, b, c (rows of cell(3,3) in C order = columns here).
! The virial is Σ r·f in the pressure sense (xx yy zz xy xz yz), as the C++ evaluator keeps it.
module caps_fft
  use, intrinsic :: iso_c_binding
  implicit none
  private
  public :: fft1d, fft3d, good_size, smallest_factor, fft_lines

  real(c_double), parameter :: pi = 3.14159265358979323846_c_double

contains

  ! the smallest n' ≥ n whose prime factors are 2, 3, 5 and 7 (radix 7 runs through the generic butterfly)
  pure integer function good_size(n) result(m)
    integer, intent(in) :: n
    integer :: k
    m = max(n, 1)
    do
      k = m
      do while (mod(k, 2) == 0); k = k / 2; end do
      do while (mod(k, 3) == 0); k = k / 3; end do
      do while (mod(k, 5) == 0); k = k / 5; end do
      do while (mod(k, 7) == 0); k = k / 7; end do
      if (k == 1) return
      m = m + 1
    end do
  end function good_size

  ! 1-D DFT: X(m) = Σ_k x(k) exp(sign·2πi·m·k/n), unnormalised, out of place (in → out). Mixed-radix decimation in
  ! time with the structure of KISS FFT (M. Borgerding): factors 4, 2, 3, 5, … ; radix-2 and radix-4 butterflies
  ! written out, others generic. tw(0:n-1) = exp(sign·2πi·j/n).
  recursive subroutine kf_work(out, ofs, x, ifs, fstride, pf, mm, lev, tw, n, sign)
    complex(c_double_complex), intent(inout) :: out(0:)
    complex(c_double_complex), intent(in) :: x(0:), tw(0:)
    integer, intent(in) :: ofs, ifs, fstride, pf(:), mm(:), lev, n, sign
    integer :: p, m, j
    p = pf(lev)
    m = mm(lev)
    if (m == 1) then
      do j = 0, p - 1
        out(ofs + j) = x(ifs + j * fstride)
      end do
    else
      do j = 0, p - 1
        call kf_work(out, ofs + j * m, x, ifs + j * fstride, fstride * p, pf, mm, lev + 1, tw, n, sign)
      end do
    end if
    select case (p)
    case (2)
      call bfly2(out, ofs, fstride, tw, m)
    case (4)
      call bfly4(out, ofs, fstride, tw, m, sign)
    case (3)
      call bfly3(out, ofs, fstride, tw, m, n)
    case (5)
      call bfly5(out, ofs, fstride, tw, m, n)
    case default
      call bfly_generic(out, ofs, fstride, tw, m, p, n)
    end select
  end subroutine kf_work

  subroutine bfly2(f, o, fs, tw, m)
    complex(c_double_complex), intent(inout) :: f(0:)
    complex(c_double_complex), intent(in) :: tw(0:)
    integer, intent(in) :: o, fs, m
    complex(c_double_complex) :: t
    integer :: k
    do k = 0, m - 1
      t = f(o + k + m) * tw(k * fs)
      f(o + k + m) = f(o + k) - t
      f(o + k) = f(o + k) + t
    end do
  end subroutine bfly2

  subroutine bfly4(f, o, fs, tw, m, sign)
    complex(c_double_complex), intent(inout) :: f(0:)
    complex(c_double_complex), intent(in) :: tw(0:)
    integer, intent(in) :: o, fs, m, sign
    complex(c_double_complex) :: s0, s1, s2, s3, s4, s5, is4
    integer :: k
    do k = 0, m - 1
      s0 = f(o + k + m) * tw(k * fs)
      s1 = f(o + k + 2 * m) * tw(2 * k * fs)
      s2 = f(o + k + 3 * m) * tw(3 * k * fs)
      s5 = f(o + k) - s1
      f(o + k) = f(o + k) + s1
      s3 = s0 + s2
      s4 = s0 - s2
      f(o + k + 2 * m) = f(o + k) - s3
      f(o + k) = f(o + k) + s3
      is4 = cmplx(-aimag(s4), real(s4), c_double_complex)   ! i·s4
      if (sign > 0) then
        f(o + k + m) = s5 + is4
        f(o + k + 3 * m) = s5 - is4
      else
        f(o + k + m) = s5 - is4
        f(o + k + 3 * m) = s5 + is4
      end if
    end do
  end subroutine bfly4

  ! radix 3 (KISS FFT kf_bfly3): tw(n/3) = exp(sign·2πi/3)
  subroutine bfly3(f, o, fs, tw, m, n)
    complex(c_double_complex), intent(inout) :: f(0:)
    complex(c_double_complex), intent(in) :: tw(0:)
    integer, intent(in) :: o, fs, m, n
    complex(c_double_complex) :: s0, s1, s2, s3, w3
    real(c_double) :: epi3
    integer :: k
    w3 = tw(fs * m)
    epi3 = aimag(w3)
    do k = 0, m - 1
      s1 = f(o + k + m) * tw(k * fs)
      s2 = f(o + k + 2 * m) * tw(2 * k * fs)
      s3 = s1 + s2
      s0 = s1 - s2
      f(o + k + m) = f(o + k) - 0.5_c_double * s3
      s0 = s0 * epi3
      f(o + k) = f(o + k) + s3
      f(o + k + 2 * m) = f(o + k + m) + cmplx(aimag(s0), -real(s0), c_double_complex)
      f(o + k + m) = f(o + k + m) + cmplx(-aimag(s0), real(s0), c_double_complex)
    end do
  end subroutine bfly3

  ! radix 5 (KISS FFT kf_bfly5): ya = exp(sign·2πi/5), yb = exp(sign·4πi/5)
  subroutine bfly5(f, o, fs, tw, m, n)
    complex(c_double_complex), intent(inout) :: f(0:)
    complex(c_double_complex), intent(in) :: tw(0:)
    integer, intent(in) :: o, fs, m, n
    complex(c_double_complex) :: s0, s1, s2, s3, s4, s7, s8, s9, s10, s5, s6, s11, s12, ya, yb
    integer :: u
    ya = tw(fs * m)
    yb = tw(2 * fs * m)
    do u = 0, m - 1
      s0 = f(o + u)
      s1 = f(o + u + m) * tw(u * fs)
      s2 = f(o + u + 2 * m) * tw(2 * u * fs)
      s3 = f(o + u + 3 * m) * tw(3 * u * fs)
      s4 = f(o + u + 4 * m) * tw(4 * u * fs)
      s7 = s1 + s4
      s10 = s1 - s4
      s8 = s2 + s3
      s9 = s2 - s3
      f(o + u) = s0 + s7 + s8
      s5 = s0 + cmplx(real(s7) * real(ya) + real(s8) * real(yb), aimag(s7) * real(ya) + aimag(s8) * real(yb), c_double_complex)
      s6 = cmplx(aimag(s10) * aimag(ya) + aimag(s9) * aimag(yb), -(real(s10) * aimag(ya)) - real(s9) * aimag(yb), c_double_complex)
      f(o + u + m) = s5 - s6
      f(o + u + 4 * m) = s5 + s6
      s11 = s0 + cmplx(real(s7) * real(yb) + real(s8) * real(ya), aimag(s7) * real(yb) + aimag(s8) * real(ya), c_double_complex)
      s12 = cmplx(-(aimag(s10) * aimag(yb)) + aimag(s9) * aimag(ya), real(s10) * aimag(yb) - real(s9) * aimag(ya), c_double_complex)
      f(o + u + 2 * m) = s11 + s12
      f(o + u + 3 * m) = s11 - s12
    end do
  end subroutine bfly5

  subroutine bfly_generic(f, o, fs, tw, m, p, n)
    complex(c_double_complex), intent(inout) :: f(0:)
    complex(c_double_complex), intent(in) :: tw(0:)
    integer, intent(in) :: o, fs, m, p, n
    complex(c_double_complex) :: sc(0:p - 1)
    integer :: u, q1, q, k, idx
    do u = 0, m - 1
      do q1 = 0, p - 1
        sc(q1) = f(o + u + q1 * m)
      end do
      do q1 = 0, p - 1
        k = u + q1 * m
        f(o + k) = sc(0)
        idx = 0
        do q = 1, p - 1
          idx = idx + fs * k
          idx = mod(idx, n)
          f(o + k) = f(o + k) + sc(q) * tw(idx)
        end do
      end do
    end do
  end subroutine bfly_generic

  ! radices of n (4 first, then 2, 3, 5, …) and the remaining length after each
  pure subroutine factorize(n, pf, mm, nf)
    integer, intent(in) :: n
    integer, intent(out) :: pf(32), mm(32), nf
    integer :: r, p
    r = n
    nf = 0
    p = 4
    do while (r > 1)
      do while (mod(r, p) /= 0)
        select case (p)
        case (4)
          p = 2
        case (2)
          p = 3
        case default
          p = p + 2
        end select
        if (p * p > r) p = r
      end do
      r = r / p
      nf = nf + 1
      pf(nf) = p
      mm(nf) = r
    end do
  end subroutine factorize

  subroutine twiddles(n, sign, tab)
    integer, intent(in) :: n, sign
    complex(c_double_complex), intent(out) :: tab(0:n - 1)
    integer :: j
    do j = 0, n - 1
      tab(j) = cmplx(cos(2 * pi * j / n), sign * sin(2 * pi * j / n), c_double_complex)
    end do
  end subroutine twiddles

  subroutine fft1d(x, n, sign)
    integer, intent(in) :: n, sign
    complex(c_double_complex), intent(inout) :: x(0:n - 1)
    complex(c_double_complex) :: tw(0:n - 1), out(0:n - 1)
    integer :: pf(32), mm(32), nf
    if (n <= 1) return
    call twiddles(n, sign, tw)
    call factorize(n, pf, mm, nf)
    call kf_work(out, 0, x, 0, 1, pf(1:nf), mm(1:nf), 1, tw, n, sign)
    x = out
  end subroutine fft1d

  pure integer function smallest_factor(n) result(r)
    integer, intent(in) :: n
    r = 2
    do while (r * r <= n)
      if (mod(n, r) == 0) return
      r = r + 1
    end do
    r = n
  end function smallest_factor

  ! In-place 3-D transform of g(0:k1-1, 0:k2-1, 0:k3-1), one axis at a time.
  subroutine fft3d(g, k1, k2, k3, sign)
    integer, intent(in) :: k1, k2, k3, sign
    complex(c_double_complex), intent(inout) :: g(0:k1 - 1, 0:k2 - 1, 0:k3 - 1)
    complex(c_double_complex) :: t1(0:k1 - 1), t2(0:k2 - 1), t3(0:k3 - 1)
    complex(c_double_complex) :: line(0:max(k1, k2, k3) - 1), res(0:max(k1, k2, k3) - 1)
    integer :: p1(32), m1(32), n1, p2(32), m2(32), n2, p3(32), m3(32), n3
    integer :: i, j, k
    call twiddles(k1, sign, t1); call factorize(k1, p1, m1, n1)
    call twiddles(k2, sign, t2); call factorize(k2, p2, m2, n2)
    call twiddles(k3, sign, t3); call factorize(k3, p3, m3, n3)
    if (k1 > 1) then
      do k = 0, k3 - 1
        do j = 0, k2 - 1
          line(0:k1 - 1) = g(:, j, k)
          call kf_work(res, 0, line, 0, 1, p1(1:n1), m1(1:n1), 1, t1, k1, sign)
          g(:, j, k) = res(0:k1 - 1)
        end do
      end do
    end if
    if (k2 > 1) then
      do k = 0, k3 - 1
        do i = 0, k1 - 1
          line(0:k2 - 1) = g(i, :, k)
          call kf_work(res, 0, line, 0, 1, p2(1:n2), m2(1:n2), 1, t2, k2, sign)
          g(i, :, k) = res(0:k2 - 1)
        end do
      end do
    end if
    if (k3 > 1) then
      do j = 0, k2 - 1
        do i = 0, k1 - 1
          line(0:k3 - 1) = g(i, j, :)
          call kf_work(res, 0, line, 0, 1, p3(1:n3), m3(1:n3), 1, t3, k3, sign)
          g(i, j, :) = res(0:k3 - 1)
        end do
      end do
    end if
  end subroutine fft3d

  ! The lines along one axis (1, 2, 3) whose line index lies in [first, last): for axis 1 the lines are (j, k) with
  ! index j + k2·k, for axis 2 (i, k) with i + k1·k, for axis 3 (i, j) with i + k1·j. Several threads may transform
  ! disjoint ranges at once: nothing here is saved between calls.
  subroutine fft_lines(g, k1, k2, k3, axis, sign, first, last)
    integer, intent(in) :: k1, k2, k3, axis, sign, first, last
    complex(c_double_complex), intent(inout) :: g(0:k1 - 1, 0:k2 - 1, 0:k3 - 1)
    complex(c_double_complex), allocatable :: tw(:), line(:), res(:)
    integer :: pf(32), mm(32), nf, n, l, i, j, k
    select case (axis)
    case (1)
      n = k1
    case (2)
      n = k2
    case default
      n = k3
    end select
    if (n <= 1) return
    allocate (tw(0:n - 1), line(0:n - 1), res(0:n - 1))
    call twiddles(n, sign, tw)
    call factorize(n, pf, mm, nf)
    do l = first, last - 1
      select case (axis)
      case (1)
        j = mod(l, k2); k = l / k2
        line = g(:, j, k)
      case (2)
        i = mod(l, k1); k = l / k1
        line = g(i, :, k)
      case default
        i = mod(l, k1); j = l / k1
        line = g(i, j, :)
      end select
      call kf_work(res, 0, line, 0, 1, pf(1:nf), mm(1:nf), 1, tw, n, sign)
      select case (axis)
      case (1)
        g(:, j, k) = res
      case (2)
        g(i, :, k) = res
      case default
        g(i, j, :) = res
      end select
    end do
  end subroutine fft_lines

end module caps_fft


module caps_kspace
  use, intrinsic :: iso_c_binding
  use caps_fft
  implicit none
  private
  public :: caps_f_pme, caps_f_ewald, caps_f_fft3d, caps_f_good_size, caps_f_pme_spread, caps_f_fft_lines, caps_f_pme_influence, &
            caps_f_pme_gather, caps_f_bspline_moduli

  real(c_double), parameter :: pi = 3.14159265358979323846_c_double

contains

  integer(c_int) function caps_f_good_size(n) bind(C, name="caps_f_good_size")
    integer(c_int), value :: n
    caps_f_good_size = good_size(n)
  end function caps_f_good_size

  subroutine caps_f_fft3d(re, im, k1, k2, k3, sign) bind(C, name="caps_f_fft3d")
    integer(c_int), value :: k1, k2, k3, sign
    real(c_double), intent(inout) :: re(k1, k2, k3), im(k1, k2, k3)
    complex(c_double_complex), allocatable :: g(:, :, :)
    allocate (g(0:k1 - 1, 0:k2 - 1, 0:k3 - 1))
    g = cmplx(re, im, c_double_complex)
    call fft3d(g, k1, k2, k3, sign)
    re = real(g)
    im = aimag(g)
  end subroutine caps_f_fft3d

  ! h: cell edge vectors as columns (h(:,1) = a …); hi = h⁻¹, whose rows are the reciprocal vectors b*_1 … (no 2π)
  subroutine cell_inverse(cell, h, hi, vol)
    real(c_double), intent(in) :: cell(3, 3)
    real(c_double), intent(out) :: h(3, 3), hi(3, 3), vol
    h = cell                                ! C row a = Fortran column 1
    vol = h(1, 1) * (h(2, 2) * h(3, 3) - h(3, 2) * h(2, 3)) - h(1, 2) * (h(2, 1) * h(3, 3) - h(3, 1) * h(2, 3)) &
          + h(1, 3) * (h(2, 1) * h(3, 2) - h(3, 1) * h(2, 2))
    hi(1, 1) = (h(2, 2) * h(3, 3) - h(2, 3) * h(3, 2)) / vol
    hi(1, 2) = (h(1, 3) * h(3, 2) - h(1, 2) * h(3, 3)) / vol
    hi(1, 3) = (h(1, 2) * h(2, 3) - h(1, 3) * h(2, 2)) / vol
    hi(2, 1) = (h(2, 3) * h(3, 1) - h(2, 1) * h(3, 3)) / vol
    hi(2, 2) = (h(1, 1) * h(3, 3) - h(1, 3) * h(3, 1)) / vol
    hi(2, 3) = (h(1, 3) * h(2, 1) - h(1, 1) * h(2, 3)) / vol
    hi(3, 1) = (h(2, 1) * h(3, 2) - h(2, 2) * h(3, 1)) / vol
    hi(3, 2) = (h(1, 2) * h(3, 1) - h(1, 1) * h(3, 2)) / vol
    hi(3, 3) = (h(1, 1) * h(2, 2) - h(1, 2) * h(2, 1)) / vol
    vol = abs(vol)
  end subroutine cell_inverse

  ! Cardinal B-spline weights M_p(w + k), k = 0 … p−1, for the fractional offset w in [0, 1), and their derivatives
  ! (Essmann et al. 1995, eqs 4.1 and 4.4). theta(j) = M_p(w + p − j), the weight of grid point floor(u) − p + j.
  pure subroutine bspline(w, p, theta, dtheta)
    real(c_double), intent(in) :: w
    integer, intent(in) :: p
    real(c_double), intent(out) :: theta(p), dtheta(p)
    real(c_double) :: div
    integer :: k, j
    theta = 0
    theta(p) = 0
    theta(2) = w
    theta(1) = 1 - w
    do k = 3, p - 1
      div = 1.0_c_double / (k - 1)
      theta(k) = div * w * theta(k - 1)
      do j = 1, k - 2
        theta(k - j) = div * ((w + j) * theta(k - j - 1) + (k - j - w) * theta(k - j))
      end do
      theta(1) = div * (1 - w) * theta(1)
    end do
    ! derivatives from the order p − 1 values
    dtheta(1) = -theta(1)
    do j = 2, p
      dtheta(j) = theta(j - 1) - theta(j)
    end do
    ! the last order
    div = 1.0_c_double / (p - 1)
    theta(p) = div * w * theta(p - 1)
    do j = 1, p - 2
      theta(p - j) = div * ((w + j) * theta(p - j - 1) + (p - j - w) * theta(p - j))
    end do
    theta(1) = div * (1 - w) * theta(1)
  end subroutine bspline

  ! |b(m)|² of the Euler exponential spline for one axis (Essmann eq. 4.8), with the usual fix at the zeros for odd p
  subroutine bspline_moduli(k, p, bsp)
    integer, intent(in) :: k, p
    real(c_double), intent(out) :: bsp(0:k - 1)
    real(c_double) :: arr(p), darr(p), sr, si, arg
    integer :: m, j
    call bspline(0.0_c_double, p, arr, darr)   ! M_p(1), M_p(2), … at the integers
    do m = 0, k - 1
      sr = 0
      si = 0
      do j = 0, p - 2
        arg = 2 * pi * m * j / real(k, c_double)
        sr = sr + arr(j + 1) * cos(arg)
        si = si + arr(j + 1) * sin(arg)
      end do
      bsp(m) = sr * sr + si * si
    end do
    do m = 0, k - 1
      if (bsp(m) < 1e-7_c_double) bsp(m) = 0.5_c_double * (bsp(modulo(m - 1, k)) + bsp(modulo(m + 1, k)))
    end do
  end subroutine bspline_moduli

  ! Smooth PME, reciprocal part. x(3,n) positions, q(n) charges, cell(3,3), grid(3), order p, beta the Ewald
  ! coefficient (1/Å). Adds the forces to f(3,n); returns the energy and the virial (xx yy zz xy xz yz).
  subroutine caps_f_pme(n, x, q, cell, grid, p, beta, kcoul, energy, f, vir) bind(C, name="caps_f_pme")
    integer(c_int), value :: n, p
    real(c_double), intent(in) :: x(3, n), q(n), cell(3, 3)
    integer(c_int), intent(in) :: grid(3)
    real(c_double), value :: beta, kcoul
    real(c_double), intent(out) :: energy, vir(6)
    real(c_double), intent(inout) :: f(3, n)
    real(c_double) :: h(3, 3), hi(3, 3), vol, s(3), u, w
    real(c_double), allocatable :: th(:, :, :), dth(:, :, :), b1(:), b2(:), b3(:), phi(:, :, :)
    integer, allocatable :: idx(:, :, :)
    complex(c_double_complex), allocatable :: g(:, :, :)
    integer :: i, d, a, b, c, k1, k2, k3, m1, m2, m3, ia, ib, ic
    real(c_double) :: mv(3), m2len, fac, e_m, eterm, qa, grad(3), gs(3), st, ss
    k1 = grid(1); k2 = grid(2); k3 = grid(3)
    call cell_inverse(cell, h, hi, vol)
    allocate (th(p, 3, n), dth(p, 3, n), idx(p, 3, n))
    allocate (g(0:k1 - 1, 0:k2 - 1, 0:k3 - 1), phi(0:k1 - 1, 0:k2 - 1, 0:k3 - 1))
    ! weights and grid points of every atom
    do i = 1, n
      s = matmul(hi, x(:, i))
      do d = 1, 3
        u = (s(d) - floor(s(d))) * grid(d)
        w = u - floor(u)
        call bspline(w, p, th(:, d, i), dth(:, d, i))
        do a = 1, p                           ! theta(a) = M_p(w + p − a): grid point floor(u) − p + a
          idx(a, d, i) = modulo(int(floor(u)) - p + a, grid(d))
        end do
      end do
    end do
    ! spread the charges
    g = (0.0_c_double, 0.0_c_double)
    do i = 1, n
      do c = 1, p
        ic = idx(c, 3, i)
        do b = 1, p
          ib = idx(b, 2, i)
          st = q(i) * th(b, 2, i) * th(c, 3, i)
          do a = 1, p
            ia = idx(a, 1, i)
            g(ia, ib, ic) = g(ia, ib, ic) + st * th(a, 1, i)
          end do
        end do
      end do
    end do
    call fft3d(g, k1, k2, k3, +1)
    ! the influence function G(m) = kcoul B(m) exp(−π²m²/β²) / (π V m²); energy ½ Σ G |F(Q)|²; virial per m
    allocate (b1(0:k1 - 1), b2(0:k2 - 1), b3(0:k3 - 1))
    call bspline_moduli(k1, p, b1)
    call bspline_moduli(k2, p, b2)
    call bspline_moduli(k3, p, b3)
    energy = 0
    vir = 0
    fac = pi * pi / (beta * beta)
    do c = 0, k3 - 1
      m3 = c
      if (m3 >= (k3 + 1) / 2) m3 = m3 - k3
      do b = 0, k2 - 1
        m2 = b
        if (m2 >= (k2 + 1) / 2) m2 = m2 - k2
        do a = 0, k1 - 1
          m1 = a
          if (m1 >= (k1 + 1) / 2) m1 = m1 - k1
          if (a == 0 .and. b == 0 .and. c == 0) then
            g(a, b, c) = 0
            cycle
          end if
          mv = m1 * hi(1, :) + m2 * hi(2, :) + m3 * hi(3, :)
          m2len = dot_product(mv, mv)
          eterm = kcoul * exp(-fac * m2len) / (pi * vol * m2len * b1(a) * b2(b) * b3(c))
          e_m = 0.5_c_double * eterm * (real(g(a, b, c))**2 + aimag(g(a, b, c))**2)
          energy = energy + e_m
          ss = 2 * (1 + fac * m2len) / m2len
          vir(1) = vir(1) + e_m * (1 - ss * mv(1) * mv(1))
          vir(2) = vir(2) + e_m * (1 - ss * mv(2) * mv(2))
          vir(3) = vir(3) + e_m * (1 - ss * mv(3) * mv(3))
          vir(4) = vir(4) - e_m * ss * mv(1) * mv(2)
          vir(5) = vir(5) - e_m * ss * mv(1) * mv(3)
          vir(6) = vir(6) - e_m * ss * mv(2) * mv(3)
          g(a, b, c) = g(a, b, c) * eterm
        end do
      end do
    end do
    ! φ = the convolution θ_rec ⋆ Q; F_i = −q_i Σ ∇M_i(k) φ(k)
    call fft3d(g, k1, k2, k3, -1)
    phi = real(g)
    do i = 1, n
      gs = 0
      do c = 1, p
        ic = idx(c, 3, i)
        do b = 1, p
          ib = idx(b, 2, i)
          do a = 1, p
            ia = idx(a, 1, i)
            st = phi(ia, ib, ic)
            gs(1) = gs(1) + dth(a, 1, i) * th(b, 2, i) * th(c, 3, i) * st
            gs(2) = gs(2) + th(a, 1, i) * dth(b, 2, i) * th(c, 3, i) * st
            gs(3) = gs(3) + th(a, 1, i) * th(b, 2, i) * dth(c, 3, i) * st
          end do
        end do
      end do
      ! ∂/∂r = Σ_d K_d (row d of h⁻¹) ∂/∂u_d
      qa = q(i)
      grad = grid(1) * gs(1) * hi(1, :) + grid(2) * gs(2) * hi(2, :) + grid(3) * gs(3) * hi(3, :)
      f(:, i) = f(:, i) - qa * grad
    end do
  end subroutine caps_f_pme

  ! Plain Ewald reciprocal sum over |m_d| ≤ kmax(d): E = kcoul/(2πV) Σ_{m≠0} exp(−π²m²/β²)/m² |S(m)|².
  subroutine caps_f_ewald(n, x, q, cell, kmax, beta, kcoul, energy, f, vir) bind(C, name="caps_f_ewald")
    integer(c_int), value :: n
    real(c_double), intent(in) :: x(3, n), q(n), cell(3, 3)
    integer(c_int), intent(in) :: kmax(3)
    real(c_double), value :: beta, kcoul
    real(c_double), intent(out) :: energy, vir(6)
    real(c_double), intent(inout) :: f(3, n)
    real(c_double) :: h(3, 3), hi(3, 3), vol, mv(3), m2len, a, e_m, ss, sr, si, arg, fac, cs, sn
    integer :: m1, m2, m3, i
    call cell_inverse(cell, h, hi, vol)
    energy = 0
    vir = 0
    fac = pi * pi / (beta * beta)
    do m1 = -kmax(1), kmax(1)
      do m2 = -kmax(2), kmax(2)
        do m3 = -kmax(3), kmax(3)
          if (m1 == 0 .and. m2 == 0 .and. m3 == 0) cycle
          mv = m1 * hi(1, :) + m2 * hi(2, :) + m3 * hi(3, :)
          m2len = dot_product(mv, mv)
          a = kcoul * exp(-fac * m2len) / (2 * pi * vol * m2len)
          sr = 0
          si = 0
          do i = 1, n
            arg = 2 * pi * dot_product(mv, x(:, i))
            sr = sr + q(i) * cos(arg)
            si = si + q(i) * sin(arg)
          end do
          e_m = a * (sr * sr + si * si)
          energy = energy + e_m
          ss = 2 * (1 + fac * m2len) / m2len
          vir(1) = vir(1) + e_m * (1 - ss * mv(1) * mv(1))
          vir(2) = vir(2) + e_m * (1 - ss * mv(2) * mv(2))
          vir(3) = vir(3) + e_m * (1 - ss * mv(3) * mv(3))
          vir(4) = vir(4) - e_m * ss * mv(1) * mv(2)
          vir(5) = vir(5) - e_m * ss * mv(1) * mv(3)
          vir(6) = vir(6) - e_m * ss * mv(2) * mv(3)
          ! F_j = 4π a q_j m Im(S* e^{iθ_j}) … = 4π a q_j m (sr sin θ_j − si cos θ_j)
          do i = 1, n
            arg = 2 * pi * dot_product(mv, x(:, i))
            cs = cos(arg)
            sn = sin(arg)
            f(:, i) = f(:, i) + 4 * pi * a * q(i) * mv * (sr * sn - si * cs)
          end do
        end do
      end do
    end do
  end subroutine caps_f_ewald

  ! ---- staged PME, for the caller's threads: spread (serial) → FFT lines (+1) → influence (planes) → FFT lines (−1)
  !      → gather (atoms). th, dth (p·3·n) and idx (p·3·n) carry the B-spline weights from spread to gather.

  subroutine caps_f_bspline_moduli(k, p, bsp) bind(C, name="caps_f_bspline_moduli")
    integer(c_int), value :: k, p
    real(c_double), intent(out) :: bsp(0:k - 1)
    call bspline_moduli(k, p, bsp)
  end subroutine caps_f_bspline_moduli

  subroutine caps_f_pme_spread(n, x, q, cell, grid, p, g, th, dth, idx) bind(C, name="caps_f_pme_spread")
    integer(c_int), value :: n, p
    real(c_double), intent(in) :: x(3, n), q(n), cell(3, 3)
    integer(c_int), intent(in) :: grid(3)
    complex(c_double_complex), intent(out) :: g(0:grid(1) - 1, 0:grid(2) - 1, 0:grid(3) - 1)
    real(c_double), intent(out) :: th(p, 3, n), dth(p, 3, n)
    integer(c_int), intent(out) :: idx(p, 3, n)
    real(c_double) :: h(3, 3), hi(3, 3), vol, s(3), u, w, st
    integer :: i, d, a, b, c
    call cell_inverse(cell, h, hi, vol)
    g = (0.0_c_double, 0.0_c_double)
    do i = 1, n
      s = matmul(hi, x(:, i))
      do d = 1, 3
        u = (s(d) - floor(s(d))) * grid(d)
        w = u - floor(u)
        call bspline(w, p, th(:, d, i), dth(:, d, i))
        do a = 1, p
          idx(a, d, i) = modulo(int(floor(u)) - p + a, grid(d))
        end do
      end do
      do c = 1, p
        do b = 1, p
          st = q(i) * th(b, 2, i) * th(c, 3, i)
          do a = 1, p
            g(idx(a, 1, i), idx(b, 2, i), idx(c, 3, i)) = g(idx(a, 1, i), idx(b, 2, i), idx(c, 3, i)) + st * th(a, 1, i)
          end do
        end do
      end do
    end do
  end subroutine caps_f_pme_spread

  subroutine caps_f_fft_lines(g, k1, k2, k3, axis, sign, first, last) bind(C, name="caps_f_fft_lines")
    integer(c_int), value :: k1, k2, k3, axis, sign, first, last
    complex(c_double_complex), intent(inout) :: g(0:k1 - 1, 0:k2 - 1, 0:k3 - 1)
    call fft_lines(g, k1, k2, k3, axis, sign, first, last)
  end subroutine caps_f_fft_lines

  ! planes c in [first, last) along the third axis: G(m) = kcoul B(m) exp(−π²m²/β²)/(π V m²) applied in place;
  ! adds this range's energy and virial
  subroutine caps_f_pme_influence(g, grid, cell, b1, b2, b3, beta, kcoul, first, last, energy, vir) &
      bind(C, name="caps_f_pme_influence")
    integer(c_int), intent(in) :: grid(3)
    complex(c_double_complex), intent(inout) :: g(0:grid(1) - 1, 0:grid(2) - 1, 0:grid(3) - 1)
    real(c_double), intent(in) :: cell(3, 3), b1(0:grid(1) - 1), b2(0:grid(2) - 1), b3(0:grid(3) - 1)
    real(c_double), value :: beta, kcoul
    integer(c_int), value :: first, last
    real(c_double), intent(inout) :: energy, vir(6)
    real(c_double) :: h(3, 3), hi(3, 3), vol, mv(3), m2len, fac, eterm, e_m, ss
    integer :: a, b, c, m1, m2, m3, k1, k2, k3
    k1 = grid(1); k2 = grid(2); k3 = grid(3)
    call cell_inverse(cell, h, hi, vol)
    fac = pi * pi / (beta * beta)
    do c = first, last - 1
      m3 = c
      if (m3 >= (k3 + 1) / 2) m3 = m3 - k3
      do b = 0, k2 - 1
        m2 = b
        if (m2 >= (k2 + 1) / 2) m2 = m2 - k2
        do a = 0, k1 - 1
          m1 = a
          if (m1 >= (k1 + 1) / 2) m1 = m1 - k1
          if (a == 0 .and. b == 0 .and. c == 0) then
            g(a, b, c) = 0
            cycle
          end if
          mv = m1 * hi(1, :) + m2 * hi(2, :) + m3 * hi(3, :)
          m2len = dot_product(mv, mv)
          eterm = kcoul * exp(-fac * m2len) / (pi * vol * m2len * b1(a) * b2(b) * b3(c))
          e_m = 0.5_c_double * eterm * (real(g(a, b, c))**2 + aimag(g(a, b, c))**2)
          energy = energy + e_m
          ss = 2 * (1 + fac * m2len) / m2len
          vir(1) = vir(1) + e_m * (1 - ss * mv(1) * mv(1))
          vir(2) = vir(2) + e_m * (1 - ss * mv(2) * mv(2))
          vir(3) = vir(3) + e_m * (1 - ss * mv(3) * mv(3))
          vir(4) = vir(4) - e_m * ss * mv(1) * mv(2)
          vir(5) = vir(5) - e_m * ss * mv(1) * mv(3)
          vir(6) = vir(6) - e_m * ss * mv(2) * mv(3)
          g(a, b, c) = g(a, b, c) * eterm
        end do
      end do
    end do
  end subroutine caps_f_pme_influence

  ! atoms [first, last) (0-based): F_i = −q_i Σ ∇M_i(k) φ(k), φ = Re g after the inverse transform; adds to f
  subroutine caps_f_pme_gather(n, q, cell, grid, p, g, th, dth, idx, first, last, f) bind(C, name="caps_f_pme_gather")
    integer(c_int), value :: n, p, first, last
    real(c_double), intent(in) :: q(n), cell(3, 3)
    integer(c_int), intent(in) :: grid(3)
    complex(c_double_complex), intent(in) :: g(0:grid(1) - 1, 0:grid(2) - 1, 0:grid(3) - 1)
    real(c_double), intent(in) :: th(p, 3, n), dth(p, 3, n)
    integer(c_int), intent(in) :: idx(p, 3, n)
    real(c_double), intent(inout) :: f(3, n)
    real(c_double) :: h(3, 3), hi(3, 3), vol, gs(3), st, grad(3), tb, tc, db, dc
    integer :: i, a, b, c
    call cell_inverse(cell, h, hi, vol)
    do i = first + 1, last
      gs = 0
      do c = 1, p
        tc = th(c, 3, i)
        dc = dth(c, 3, i)
        do b = 1, p
          tb = th(b, 2, i)
          db = dth(b, 2, i)
          do a = 1, p
            st = real(g(idx(a, 1, i), idx(b, 2, i), idx(c, 3, i)))
            gs(1) = gs(1) + dth(a, 1, i) * tb * tc * st
            gs(2) = gs(2) + th(a, 1, i) * db * tc * st
            gs(3) = gs(3) + th(a, 1, i) * tb * dc * st
          end do
        end do
      end do
      grad = grid(1) * gs(1) * hi(1, :) + grid(2) * gs(2) * hi(2, :) + grid(3) * gs(3) * hi(3, :)
      f(:, i) = f(:, i) - q(i) * grad
    end do
  end subroutine caps_f_pme_gather

end module caps_kspace
