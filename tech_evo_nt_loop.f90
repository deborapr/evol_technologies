! Fortran program to simulate technology evolution - Sweep nt from 10^2 to 10^6 for fixed q values
! Saves only final values per nt in a single file for each q
! Also saves aggregates per nt for each q
! By Debora Princepe - 2025-07-16
! run: gfortran -O3 -march=native -funroll-loops -ffast-math tech_evo_nt_loop.f90 -o ntloop


program tech_evo_nt_loop
    implicit none

    integer, parameter    :: dp = selected_real_kind(15)
    integer, parameter    :: n_nt = 100
    integer, parameter    :: n_q = 1 
    ! real(dp), dimension(n_q), parameter :: q_array = [0.25_dp, 0.5_dp, 0.75_dp]
    real(dp), dimension(n_q), parameter :: q_array = [0.85_dp]
    integer, dimension(n_nt) :: nt_array

    integer, allocatable  :: seed(:)
    real(dp), allocatable :: gamma_i(:), theta_i(:), tau_i(:)
    character(len=100)    :: filename_final, filename_agg, timestamp
    character(len=32)     :: arg
    integer, allocatable  :: sorted_idx(:)
    integer               :: values(8), num_args
    integer :: seed_size, i, nt, nt_idx, q_idx, unit_final, unit_agg, ns, iter, max_iter= 10000
    logical :: converged_E
    real(dp) :: mean = 0.5_dp, q, d, E0, N0, den, A, B
    real(dp) :: Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, Tau
    real(dp) :: N, E, N_new, E_new
    real(dp) :: tol = 1e-5_dp
    logical :: use_saved_seed = .true.
    character(len=20) :: seed_filename = 'random_seed.dat'

    ! Data structures
    integer, allocatable  :: valid_tech(:)
    integer, allocatable  :: pending_tech(:)
    integer :: num_valid, num_pending
    real(dp), allocatable :: N_sol(:), E_sol(:), tau_sol(:)

    ! Loop parameters
    real(dp) :: nt_min = 1000.0_dp, nt_max = 10000000.0_dp

    ! Precomputed constants
    real(dp) :: inv_E0, N0_over_E0

    ! Cache for root finding
    real(dp) :: last_E_solution
    logical :: has_cached_E

    ! Timing variables
    real(dp) :: cpu_start, cpu_end, elapsed_time

    ! Read parameters (d, E0, N0 only - q and nt will be looped)
    num_args = command_argument_count()
    if (num_args /= 3) then
        print *, "Use: ./techevo_nt_loop d E0 N0"
        print *, "Example: ./techevo_nt_loop 0.8 10.0 1000.0"
        print *, "Note: nt will be varied from 10^2 to 10^6, q will be {0.25, 0.5, 0.75}"
        stop
    end if
    call get_command_argument(1, arg); read(arg, *) d
    call get_command_argument(2, arg); read(arg, *) E0
    call get_command_argument(3, arg); read(arg, *) N0

    call date_and_time(VALUES=values)
    write(timestamp, '(I4,2I2.2,A,2I2.2,2I2.2)') &
        values(1), values(2), values(3), '_', &
        values(5), values(6), values(7), values(8)/100

    call handle_random_seed(seed, seed_filename, use_saved_seed)

    do i = 1, n_nt
        nt_array(i) = int(exp(log(nt_min) + (log(nt_max) - log(nt_min)) * (i-1) / (n_nt-1)))
    end do

    print *, "Sweeping nt from 10^2 to 10^6 for q values: ", q_array

    ! Main loops: q values, then nt values
    do q_idx = 1, n_q
        q = q_array(q_idx)

        ! Create output files for this q
        write(filename_final, '(A,F5.2,A,A,A)') 'final_values_q_', q, '_', trim(timestamp), '.dat'
        write(filename_agg, '(A,F5.2,A,A,A)') 'final_sums_q_', q, '_', trim(timestamp), '.dat'
        open(newunit=unit_final, file=filename_final, status='replace')
        write(unit_final, '(A,F8.5)') '# Final values per nt for q = ', q
        write(unit_final, '(A,3F12.5)') '# Parameters: d, E0, N0 = ', d, E0, N0
        write(unit_final, '(A)') '# Format: nt, ns_final, N_final, E_final, tau_final'

        open(newunit=unit_agg, file=filename_agg, status='replace')
        write(unit_agg, '(A,F8.5)') '# Sums per nt for q = ', q
        write(unit_agg, '(A,3F12.5)') '# Parameters: d, E0, N0 = ', d, E0, N0
        write(unit_agg, '(A)') '# Format: nt, Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega'

        do nt_idx = 1, n_nt
            nt = nt_array(nt_idx)

            ! Allocate arrays for current nt
            if (allocated(gamma_i)) deallocate(gamma_i, theta_i, tau_i, sorted_idx)
            if (allocated(valid_tech)) deallocate(valid_tech, pending_tech)
            if (allocated(N_sol)) deallocate(N_sol, E_sol, tau_sol)

            allocate(gamma_i(nt), theta_i(nt), tau_i(nt), sorted_idx(nt))
            allocate(valid_tech(nt), pending_tech(nt))
            allocate(N_sol(nt), E_sol(nt), tau_sol(nt))

            ! Precompute constants
            inv_E0 = 1.0_dp / E0
            N0_over_E0 = N0 * inv_E0

            ! Generate exponential random variables
            call generate_exponential_array(gamma_i, mean, nt)
            call generate_exponential_array(theta_i, mean, nt)
            tau_i = gamma_i / theta_i

            ! Sort by tau_i
            sorted_idx = [(i, i = 1, nt)]
            call quicksort_optimized(tau_i, sorted_idx, 1, nt)
            call reorder_by_index(gamma_i, sorted_idx, nt)
            call reorder_by_index(theta_i, sorted_idx, nt)
            call reorder_by_index(tau_i, sorted_idx, nt)

            ! Initialize for this q and nt
            num_valid = 0
            num_pending = 0
            ns = 0
            Gamma = 0.0_dp; Theta = 0.0_dp; Tau = 0.0_dp
            Delta = 0.0_dp; Phi = 0.0_dp; Psi = 0.0_dp
            Chi = 0.0_dp; Sigma = 0.0_dp; Omega = 0.0_dp
            N = N0; E = E0
            last_E_solution = 0.0_dp
            has_cached_E = .false.

            call cpu_time(cpu_start)

            ! Main simulation loop
            do i = 1, nt
                if (try_add_technology(i, gamma_i(i), theta_i(i), tau_i(i), &
                                      Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, &
                                      N, E, N_new, E_new, Tau, q, d, N0, E0, &
                                      inv_E0, N0_over_E0, tol, max_iter, &
                                      last_E_solution, has_cached_E)) then

                    num_valid = num_valid + 1
                    valid_tech(num_valid) = i

                    ns = ns + 1
                    N_sol(ns) = N_new
                    E_sol(ns) = E_new
                    tau_sol(ns) = tau_i(i)

                    N = N_new; E = E_new

                    call smart_reevaluate_pending(pending_tech, num_pending, gamma_i, theta_i, tau_i, &
                        valid_tech, num_valid, &
                        Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, &
                        N, E, Tau, q, d, N0, E0, &
                        inv_E0, N0_over_E0, tol, max_iter, &
                        N_sol, E_sol, tau_sol, ns, &
                        last_E_solution, has_cached_E)
                else
                    if (num_pending < nt) then
                        num_pending = num_pending + 1
                        pending_tech(num_pending) = i
                    end if
                end if
            end do

            call cpu_time(cpu_end)
            elapsed_time = cpu_end - cpu_start

            ! Write final values to files for this q and nt
            write(unit_final, '(I8,I8,3E16.8)') nt, ns, N, E, Tau
            write(unit_agg, '(I8,8E16.8)') nt, Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega

            flush(unit_final)
            flush(unit_agg)

            if (mod(nt,1000) == 0 .or. nt <= 1000) then
                write(*, '(A,I0,A,I0,A,F8.3,A)') &
                    "nt=", nt, ", ns=", ns, " (", elapsed_time, "s)"
            end if
        end do

        close(unit_final)
        close(unit_agg)

        write(*, '(A,F5.2,A)') "Files created for q=", q, ":"
        write(*, '(A)') "  "//trim(filename_final)
        write(*, '(A)') "  "//trim(filename_agg)
        print *, ""
    end do

    print *, "All simulations completed!"

    ! Cleanup
    deallocate(seed, gamma_i, theta_i, tau_i, sorted_idx)
    deallocate(valid_tech, pending_tech, N_sol, E_sol, tau_sol)

contains


!=============================== subroutines ===============================!

    subroutine handle_random_seed(seed, filename, use_saved)
        integer, allocatable, intent(out) :: seed(:)
        character(len=*), intent(in) :: filename
        logical, intent(in) :: use_saved
        integer :: seed_size

        if (use_saved) then
            call load_seed_from_file(seed, filename)
            if (allocated(seed)) then
                call random_seed(put=seed)
                return
            end if
        end if

        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        call random_seed(get=seed)
        if (.not. use_saved) call save_seed_to_file(seed, filename)
    end subroutine handle_random_seed

    subroutine generate_exponential_array(arr, lambda, n)
        integer, intent(in) :: n
        real(dp), intent(in) :: lambda
        real(dp), intent(out) :: arr(n)
        real(dp) :: u_arr(n)
        integer :: i
        
        call random_number(u_arr)
        do i = 1, n
            arr(i) = -lambda * log(u_arr(i))
        end do
    end subroutine generate_exponential_array

    subroutine reorder_by_index(arr, idx, n)
        integer, intent(in) :: n, idx(n)
        real(dp), intent(inout) :: arr(n)
        real(dp) :: temp_arr(n)
        integer :: i
        
        temp_arr = arr
        do i = 1, n
            arr(i) = temp_arr(idx(i))
        end do
    end subroutine reorder_by_index

    logical function try_add_technology(tech_idx, gamma_val, theta_val, tau_val, &
        Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, &
        N, E, N_new, E_new, Tau, q, d, N0, E0, &
        inv_E0, N0_over_E0, tol, max_iter, &
        cached_E, has_cache) result(success)
        integer, intent(in) :: tech_idx, max_iter
        real(dp), intent(in) :: gamma_val, theta_val, tau_val
        real(dp), intent(inout) :: Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega
        real(dp), intent(in) :: N, E, q, d, N0, E0, inv_E0, N0_over_E0, tol
        real(dp), intent(out) :: N_new, E_new, Tau
        real(dp), intent(inout) :: cached_E
        logical, intent(inout) :: has_cache
        logical :: converged_E
        real(dp) :: den, A, B, temp_gamma, temp_theta, temp_delta, temp_phi, temp_psi

        success = .false.

        temp_gamma = Gamma + gamma_val
        temp_theta = Theta + theta_val
        temp_delta = Delta + theta_val * gamma_val
        temp_phi = Phi + gamma_val * gamma_val
        temp_psi = Psi + theta_val * theta_val

        Chi = temp_psi / (1.0_dp + temp_theta)
        Sigma = temp_psi * temp_gamma / (1.0_dp + temp_theta) - temp_delta
        Omega = temp_phi - temp_gamma * temp_delta / (1.0_dp + temp_theta)

        if ((1.0_dp - q*Chi) <= 0.0_dp) return

        den = (1.0_dp + temp_theta) * (1.0_dp - q*Chi)
        A = d * temp_delta * N0_over_E0 / den
        B = d * (q * temp_delta * Sigma / den - Omega)

        if (has_cache) then
            E_new = cached_E
        else
            E_new = E
        end if

        call find_zero_optimized(E_new, tol, max_iter, E_new, converged_E, A, B, inv_E0)

        if (.not. converged_E) return

        cached_E = E_new
        has_cache = .true.

        N_new = (N0 + q*Sigma*(E0 - E_new)) / (1.0_dp - q*Chi)
        Tau = (N_new/(E0-E_new) + temp_gamma) / (1.0_dp + temp_theta)

        if (tau_val > Tau .or. N_new < 0.0_dp .or. E_new < 0.0_dp) return

        Gamma = temp_gamma
        Theta = temp_theta
        Delta = temp_delta
        Phi = temp_phi
        Psi = temp_psi
        success = .true.
    end function try_add_technology

    subroutine smart_reevaluate_pending(pending_list, num_pending, gamma_arr, theta_arr, tau_arr, &
        valid_list, num_valid, &
        Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, &
        N, E, Tau, q, d, N0, E0, &
        inv_E0, N0_over_E0, tol, max_iter, &
        N_sol, E_sol, tau_sol, ns, &
        cached_E, has_cache)
        integer, intent(inout) :: pending_list(:), num_pending
        real(dp), intent(in) :: gamma_arr(:), theta_arr(:), tau_arr(:)
        integer, intent(inout) :: valid_list(:), num_valid
        real(dp), intent(inout) :: Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega
        real(dp), intent(inout) :: N, E, Tau
        real(dp), intent(in) :: q, d, N0, E0, inv_E0, N0_over_E0, tol
        integer, intent(in) :: max_iter
        real(dp), intent(inout) :: N_sol(:), E_sol(:), tau_sol(:)
        integer, intent(inout) :: ns
        real(dp), intent(inout) :: cached_E
        logical, intent(inout) :: has_cache

        integer :: i, tech_idx, new_pending_count
        real(dp) :: N_new, E_new
        logical :: added_any

        added_any = .true.

        do while (added_any .and. num_pending > 0)
            added_any = .false.
            new_pending_count = 0

            do i = 1, num_pending
                tech_idx = pending_list(i)

                if (try_add_technology(tech_idx, gamma_arr(tech_idx), theta_arr(tech_idx), tau_arr(tech_idx), &
                    Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, &
                    N, E, N_new, E_new, Tau, q, d, N0, E0, &
                    inv_E0, N0_over_E0, tol, max_iter, &
                    cached_E, has_cache)) then

                    num_valid = num_valid + 1
                    valid_list(num_valid) = tech_idx

                    ns = ns + 1
                    N_sol(ns) = N_new
                    E_sol(ns) = E_new
                    tau_sol(ns) = tau_arr(tech_idx)

                    N = N_new; E = E_new
                    added_any = .true.
                else
                    new_pending_count = new_pending_count + 1
                    pending_list(new_pending_count) = tech_idx
                end if
            end do

            num_pending = new_pending_count
        end do
    end subroutine smart_reevaluate_pending

    recursive subroutine quicksort_optimized(arr, idx, first, last)
        integer, intent(inout) :: idx(:)
        real(dp), intent(in) :: arr(:)
        integer, intent(in) :: first, last
        integer :: i, j, pivot_idx, temp_idx
        real(dp) :: pivot

        if (last - first < 10) then
            call insertion_sort(arr, idx, first, last)
            return
        end if

        pivot_idx = median_of_three_pivot(arr, idx, first, last)
        pivot = arr(idx(pivot_idx))

        i = first; j = last
        do
            do while (i <= last .and. arr(idx(i)) < pivot)
                i = i + 1
            end do
            do while (j >= first .and. pivot < arr(idx(j)))
                j = j - 1
            end do

            if (i >= j) exit

            temp_idx = idx(i); idx(i) = idx(j); idx(j) = temp_idx
            i = i + 1; j = j - 1
        end do

        if (first < i-1) call quicksort_optimized(arr, idx, first, i-1)
        if (j+1 < last) call quicksort_optimized(arr, idx, j+1, last)
    end subroutine quicksort_optimized

    subroutine insertion_sort(arr, idx, first, last)
        real(dp), intent(in) :: arr(:)
        integer, intent(inout) :: idx(:)
        integer, intent(in) :: first, last
        integer :: i, j, key

        do i = first + 1, last
            key = idx(i)
            j = i - 1
            do while (j >= first .and. arr(idx(j)) > arr(key))
                idx(j + 1) = idx(j)
                j = j - 1
            end do
            idx(j + 1) = key
        end do
    end subroutine insertion_sort

    integer function median_of_three_pivot(arr, idx, first, last) result(pivot_idx)
        real(dp), intent(in) :: arr(:)
        integer, intent(in) :: idx(:), first, last
        integer :: mid

        mid = (first + last) / 2

        if (arr(idx(first)) > arr(idx(mid))) then
            if (arr(idx(mid)) > arr(idx(last))) then
                pivot_idx = mid
            else if (arr(idx(first)) > arr(idx(last))) then
                pivot_idx = last
            else
                pivot_idx = first
            end if
        else
            if (arr(idx(first)) > arr(idx(last))) then
                pivot_idx = first
            else if (arr(idx(mid)) > arr(idx(last))) then
                pivot_idx = last
            else
                pivot_idx = mid
            end if
        end if
    end function median_of_three_pivot

    subroutine find_zero_optimized(x0, tol, max_iter, x, converged, A, B, inv_E0)
        real(dp), intent(in) :: x0, tol, A, B, inv_E0
        integer, intent(in) :: max_iter
        real(dp), intent(out) :: x
        logical, intent(out) :: converged

        real(dp) :: fx, h, df, step, x_old
        integer :: iter
        real(dp), parameter :: eps = epsilon(1.0_dp)

        x = x0
        converged = .false.

        do iter = 1, max_iter
            x_old = x
            fx = E_residual_inline(x, A, B, inv_E0)

            if (abs(fx) < tol) then
                converged = .true.
                exit
            end if

            h = sqrt(eps) * max(abs(x), 1.0_dp)
            df = (E_residual_inline(x + h, A, B, inv_E0) - fx) / h

            if (abs(df) > sqrt(eps)) then
                step = -fx/df
                step = sign(min(abs(step), 10.0_dp * max(abs(x), 1.0_dp)), step)
                x = x + step
            else
                x = x + sign(h, -fx)
            end if

            if (abs(x - x_old) < tol * max(abs(x), 1.0_dp)) then
                converged = .true.
                exit
            end if
        end do
    end subroutine find_zero_optimized

    real(dp) function E_residual_inline(E_trial, A, B, inv_E0) result(res)
        real(dp), intent(in) :: E_trial, A, B, inv_E0
        real(dp) :: e_ratio

        e_ratio = E_trial * inv_E0
        res = e_ratio*(1.0_dp + A + B) - B*e_ratio*e_ratio - 1.0_dp
    end function E_residual_inline

    subroutine save_seed_to_file(seed, filename)
        integer, intent(in) :: seed(:)
        character(len=*), intent(in) :: filename
        integer :: unit, i

        open(newunit=unit, file=filename, status='replace')
        write(unit, '(I0)') size(seed)
        do i = 1, size(seed)
            write(unit, '(I0)') seed(i)
        end do
        close(unit)
    end subroutine save_seed_to_file

    subroutine load_seed_from_file(seed, filename)
        integer, allocatable, intent(out) :: seed(:)
        character(len=*), intent(in) :: filename
        integer :: unit, seed_size, i, ios

        open(newunit=unit, file=filename, status='old', iostat=ios)
        if (ios /= 0) return

        read(unit, *) seed_size
        allocate(seed(seed_size))
        do i = 1, seed_size
            read(unit, *) seed(i)
        end do
        close(unit)
    end subroutine load_seed_from_file

end program tech_evo_nt_loop