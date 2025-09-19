! Fortran program to simulate technology evolution
! N0 constant, E0 constant
! For each value of q, nt technologies are generated and sorted
! Then, adding one by one, checks the equilibrium until tau_i > Tau
! If (1.0_dp - q*Chi) <= 0.0, then skip this technology
! by Matteo & Debora - 2025-07-30

program igs
  implicit none
  integer, parameter :: dp = selected_real_kind(15)
  integer, parameter :: indntmin=3, indntmax=7, nnt = 100
  integer, parameter :: ntmax = 10**indntmax
  integer            :: nt_array(nnt), ii
  real(dp)           :: t(ntmax), th(ntmax), g(ntmax)
  real(dp)           :: E, N, Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Lambda, Omega
  integer            :: i, j, m, idum, ns, nreal, outunit, nt
  real(dp)           :: q, d, E0, N0, tm, gm, tau, B, A
  real(dp)           :: Es, check
  real(dp)           :: ran2
  character(len=50)  :: filename

  ! Random number generator seed
  idum = -7120187
  
  ! Set parameters
  d = 0.1_dp; N0 = 100.0_dp; E0 = 10.0_dp; tm = 0.5_dp; gm = 0.5_dp
  ! Fixed q
  q = 0.85_dp

  do i = 1, nnt
      nt_array(i) = nint(10.0_dp ** ( indntmin + (indntmax - indntmin) * real(i - 1, dp) / (nnt - 1) ))
  end do

  outunit = 40
  write(filename, '(A,I0,A)') 'output_tloop_q_', int(q*1000), '.dat'
  open(unit=outunit, file=filename, status='unknown')
  write(outunit, '(A,I0)') '# Seed: ', idum
  write(outunit, '(A,2F8.2,1X,F8.4,1X,F8.4)') '# Parameters: E0, N0, d, q = ', E0, N0, d, q
  write(outunit, '(A)') '# Format: nt, E/E0, N/N0, tau, t(ns), ns'


  ii = 0
1 nreal = 100 ! Number of realizations per nt

  ii = ii + 1
  nt = nt_array(ii)
  print *, 'Number of technologies:', nt
5 E = E0
  N = N0
  ! Draw technologies
  do i = 1, nt
     th(i) = -tm * log(ran2(idum))
     g(i)  = -gm * log(ran2(idum))
     t(i)  = g(i) / th(i)
  end do
  call sort2(nt, t, g)
  do i = 1, nt
     th(i) = g(i) / t(i)
  end do
  
  ! Initialize state variables
  Gamma = 0.0_dp; Theta = 0.0_dp; Tau = 0.0_dp
  Delta = 0.0_dp; Phi = 0.0_dp; Psi = 0.0_dp
  Chi = 0.0_dp; Sigma = 0.0_dp; Omega = 0.0_dp; Lambda = 0.0_dp

  ns = 0
  
  ! Loop over technologies
10 ns = ns + 1

  Gamma = Gamma + g(ns)
  Theta = Theta + th(ns)
  Delta = Delta + g(ns) * th(ns)
  Phi   = Phi   + g(ns) * g(ns)
  Psi   = Psi   + th(ns) * th(ns)
  Chi   = Psi / (1.0 + Theta)
  Lambda = Delta / (1.0 + Theta)
  Omega = Phi - Gamma * Delta / (1.0 + Theta)
  Sigma = Chi * Gamma - Delta

  ! Solution for E
  A = d * Delta * N0 / ((1.0 + Theta) * (1.0 - q * Chi) * E0)
  B = d * (q * Delta * Sigma / ((1.0 + Theta) * (1.0 - q * Chi)) - Omega)
  A = 1.0 + A + B
  E = E0 * (A - sqrt(A * A - 4.0 * B)) / (2.0 * B)
  N = (N0 + q * Sigma * (E0 - E)) / (1.0 - q * Chi)
  tau = (N / (E0 - E) + Gamma) / (1.0 + Theta)

  if (tau < 0.0 .or. E > E0) goto 10
  if (t(ns + 1) < tau) goto 10
  if (1.0_dp - q*Chi > 0) write(outunit,*) nt, E / E0, N / N0, tau, t(ns), ns
  nreal = nreal - 1
  if (nreal > 1) goto 5
  if (nt < ntmax) goto 1

end program igs

subroutine sort2(n, arr, brr)
  integer, parameter :: dp = selected_real_kind(15)
  integer, intent(in) :: n
  real(dp), intent(inout) :: arr(n), brr(n)
  integer, parameter :: M=7, NSTACK=50
  integer :: i, ir, j, jstack, k, l, istack(NSTACK)
  real(dp) :: a, b, temp

  jstack = 0
  l = 1
  ir = n
1 if (ir - l < M) then
      do j = l + 1, ir
        a = arr(j)
        b = brr(j)
        do i = j - 1, 1, -1
          if (arr(i) <= a) exit
          arr(i + 1) = arr(i)
          brr(i + 1) = brr(i)
        end do
        arr(i + 1) = a
        brr(i + 1) = b
      end do
      if (jstack == 0) return
      ir = istack(jstack)
      l = istack(jstack - 1)
      jstack = jstack - 2
    else
      k = (l + ir) / 2
      temp = arr(k)
      arr(k) = arr(l + 1)
      arr(l + 1) = temp
      temp = brr(k)
      brr(k) = brr(l + 1)
      brr(l + 1) = temp
      if (arr(l + 1) > arr(ir)) then
        temp = arr(l + 1)
        arr(l + 1) = arr(ir)
        arr(ir) = temp
        temp = brr(l + 1)
        brr(l + 1) = brr(ir)
        brr(ir) = temp
      end if
      if (arr(l) > arr(ir)) then
        temp = arr(l)
        arr(l) = arr(ir)
        arr(ir) = temp
        temp = brr(l)
        brr(l) = brr(ir)
        brr(ir) = temp
      end if
      if (arr(l + 1) > arr(l)) then
        temp = arr(l + 1)
        arr(l + 1) = arr(l)
        arr(l) = temp
        temp = brr(l + 1)
        brr(l + 1) = brr(l)
        brr(l) = temp
      end if
      i = l + 1
      j = ir
      a = arr(l)
      b = brr(l)
3     i = i + 1
      if (arr(i) < a) goto 3
4     j = j - 1
      if (arr(j) > a) goto 4
      if (j < i) goto 5
      temp = arr(i)
      arr(i) = arr(j)
      arr(j) = temp
      temp = brr(i)
      brr(i) = brr(j)
      brr(j) = temp
      goto 3
5     arr(l) = arr(j)
      arr(j) = a
      brr(l) = brr(j)
      brr(j) = b
      jstack = jstack + 2
      if (jstack > NSTACK) stop 'NSTACK too small in sort2'
      if (ir - i + 1 >= j - l) then
        istack(jstack) = ir
        istack(jstack - 1) = i
        ir = j - 1
      else
        istack(jstack) = j - 1
        istack(jstack - 1) = l
        l = i
      end if
    end if
    goto 1
end subroutine sort2

real(kind=selected_real_kind(15)) function ran2(idum)
  integer, parameter :: dp = selected_real_kind(15)
  integer, intent(inout) :: idum
  integer, parameter :: IM1=2147483563, IM2=2147483399
  real(dp), parameter :: AM=1.0/IM1
  integer, parameter :: IMM1=IM1-1
  integer, parameter :: IA1=40014, IA2=40692, IQ1=53668, IQ2=52774
  integer, parameter :: IR1=12211, IR2=3791
  integer, parameter :: NTAB=32, NDIV=1+IMM1/NTAB
  real(dp), parameter :: EPS=1.2e-7, RNMX=1.0-EPS
  integer :: idum2, j, k, iv(NTAB), iy
  save :: iv, iy, idum2
  data idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum <= 0) then
    idum = max(-idum, 1)
    idum2 = idum
    do j = NTAB + 8, 1, -1
      k = idum / IQ1
      idum = IA1 * (idum - k * IQ1) - k * IR1
      if (idum < 0) idum = idum + IM1
      if (j <= NTAB) iv(j) = idum
    end do
    iy = iv(1)
  end if
  k = idum / IQ1
  idum = IA1 * (idum - k * IQ1) - k * IR1
  if (idum < 0) idum = idum + IM1
  k = idum2 / IQ2
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2
  if (idum2 < 0) idum2 = idum2 + IM2
  j = 1 + iy / NDIV
  iy = iv(j) - idum2
  iv(j) = idum
  if (iy < 1) iy = iy + IMM1
  ran2 = min(AM * real(iy), RNMX)
  return
end function ran2