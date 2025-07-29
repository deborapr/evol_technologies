! Fortran program to simulate technology evolution
! N0 constant, E0 constant
! This program generates nt random technologies and sorts them
! Then, adding one by one, checks the equilibrium until tau_i > Tau
! If (1.0_dp - q*Chi) <= 0.0, then skip this technology
! by Debora Princepe - 2025-07-02

program tech_evo_optimized
    implicit none

    integer, parameter    :: nt = 1000000
    integer, parameter    :: dp = selected_real_kind(15)
    integer, allocatable  :: seed(:)
    real(dp), allocatable :: gamma_i(:), theta_i(:), tau_i(:)
    character(len=100)    :: filename, timestamp
    character(len=32)     :: arg
    integer, allocatable  :: sorted_idx(:)
    integer               :: values(8), num_args
    integer :: seed_size, i, unit, ns, iter, max_iter = 10000
    logical :: converged_E
    real(dp) :: mean = 0.5_dp
    real(dp) :: q, d, E0, N0, den, A, B
    real(dp) :: Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, Tau
    real(dp) :: N, E, N_new, E_new
    real(dp) :: tol = 1e-5_dp
    logical :: use_saved_seed = .true.
    character(len=20) :: seed_filename = 'random_seed.dat'
    
    ! Optimized data structures
    integer, allocatable  :: valid_tech(:)  ! Compact list of valid technologies
    integer, allocatable  :: pending_tech(:) ! Technologies to re-evaluate
    integer :: num_valid, num_pending
    real(dp), allocatable :: N_sol(:), E_sol(:), tau_sol(:)
    
    ! Precomputed constants
    real(dp) :: inv_E0, N0_over_E0
    
    ! Cache for expensive calculations
    real(dp) :: last_E_solution = 0.0_dp
    logical :: has_cached_E = .false.

    ! Read parameters
    num_args = command_argument_count()
    if (num_args /= 4) then
        print *, "Use: ./techevo q d E0 N0"
        print *, "Example: ./techevo 0.8 0.8 10.0 1000.0"
        stop
    end if
    call get_command_argument(1, arg); read(arg, *) q
    call get_command_argument(2, arg); read(arg, *) d
    call get_command_argument(3, arg); read(arg, *) E0
    call get_command_argument(4, arg); read(arg, *) N0

    ! Precompute constants
    inv_E0 = 1.0_dp / E0
    N0_over_E0 = N0 * inv_E0

    ! Allocate arrays
    allocate(gamma_i(nt), theta_i(nt), tau_i(nt), sorted_idx(nt))
    allocate(valid_tech(nt), pending_tech(nt))
    allocate(N_sol(nt), E_sol(nt), tau_sol(nt))

    ! Initialize
    num_valid = 0
    num_pending = 0
    sorted_idx = [(i, i = 1, nt)]

    ! Seed handling (optimized)
    call handle_random_seed(seed, seed_filename, use_saved_seed)
    
    ! Generate random variables (vectorized where possible)
    call generate_exponential_variables(gamma_i, theta_i, tau_i, nt, mean)

    ! Sort by tau_i using optimized quicksort
    call quicksort_optimized(tau_i, sorted_idx, 1, nt)

    ! Rearrange arrays efficiently (vectorized)
    call rearrange_arrays(gamma_i, theta_i, tau_i, sorted_idx, nt)

    ! Initialize state variables
    ns = 0
    Gamma = 0.0_dp; Theta = 0.0_dp; Tau = 0.0_dp
    Delta = 0.0_dp; Phi = 0.0_dp; Psi = 0.0_dp
    Chi = 0.0_dp; Sigma = 0.0_dp; Omega = 0.0_dp
    N = N0; E = E0

    ! Main optimization: process technologies in batches and use smart re-evaluation
    do i = 1, nt
        ! Try to add current technology
        if (try_add_technology(i, gamma_i(i), theta_i(i), tau_i(i), &
                              Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, &
                              N, E, N_new, E_new, Tau, q, d, N0, E0, &
                              inv_E0, N0_over_E0, tol, max_iter, &
                              last_E_solution, has_cached_E)) then
            
            num_valid = num_valid + 1
            valid_tech(num_valid) = i
            
            ! Update solution arrays
            ns = ns + 1
            N_sol(ns) = N_new
            E_sol(ns) = E_new
            tau_sol(ns) = tau_i(i)
            
            ! Update state
            N = N_new; E = E_new

            ! Smart re-evaluation: only check promising pending technologies
            call smart_reevaluate_pending(pending_tech, num_pending, gamma_i, theta_i, tau_i, &
                                        valid_tech, num_valid, &
                                        Gamma, Theta, Delta, Phi, Psi, Chi, Sigma, Omega, &
                                        N, E, Tau, q, d, N0, E0, &
                                        inv_E0, N0_over_E0, tol, max_iter, &
                                        N_sol, E_sol, tau_sol, ns, &
                                        last_E_solution, has_cached_E)

            ! ! Progress reporting (less frequent)
            ! if (mod(ns, 1000) == 0 .or. ns < 100) then
            !     write(*, '(A,I8,A,F12.4,A,F12.4,A,F12.6)') &
            !         "Tech ", ns, ": N=", N, " E=", E, " Tau=", Tau
            ! end if
        else
            ! Add to pending list for later re-evaluation
            if (num_pending < nt) then
                num_pending = num_pending + 1
                pending_tech(num_pending) = i
            end if
        end if
    end do

    ! Final summary
    write(*,*)
    write(*, '(A,I8,A,F12.4,A,F12.4)') "Final: ", ns, " technologies, N=", N, " E=", E

    ! Save results with optimized I/O
    call date_and_time(VALUES=values)
    write(timestamp, '(I4,2I2.2,A,2I2.2,2I2.2)') &
        values(1), values(2), values(3), '_', &
        values(5), values(6), values(7), values(8)/100
    
    call save_results_optimized(ns, N_sol, E_sol, tau_sol, gamma_i, theta_i, tau_i, &
                               valid_tech, num_valid, seed, E0, N0, q, d, timestamp)

    ! Cleanup
    deallocate(seed, gamma_i, theta_i, tau_i, sorted_idx)
    deallocate(valid_tech, pending_tech, N_sol, E_sol, tau_sol)

contains

    !=================== OPTIMIZED SUBROUTINES ===================

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

    subroutine generate_exponential_variables(gamma_arr, theta_arr, tau_arr, n, lambda)
        real(dp), intent(out) :: gamma_arr(:), theta_arr(:), tau_arr(:)
        integer, intent(in) :: n
        real(dp), intent(in) :: lambda
        real(dp), allocatable :: random_vals(:)
        integer :: i
        
        allocate(random_vals(2*n))
        call random_number(random_vals)
        
        ! Vectorized exponential generation
        do i = 1, n
            gamma_arr(i) = -lambda * log(random_vals(2*i-1))
            theta_arr(i) = -lambda * log(random_vals(2*i))
            tau_arr(i) = gamma_arr(i) / theta_arr(i)
        end do
        
        deallocate(random_vals)
    end subroutine generate_exponential_variables

    subroutine rearrange_arrays(gamma_arr, theta_arr, tau_arr, idx, n)
        real(dp), intent(inout) :: gamma_arr(:), theta_arr(:), tau_arr(:)
        integer, intent(in) :: idx(:)
        integer, intent(in) :: n
        real(dp), allocatable :: temp_g(:), temp_t(:), temp_tau(:)
        
        allocate(temp_g(n), temp_t(n), temp_tau(n))
        
        ! Vectorized rearrangement
        temp_g = gamma_arr(idx)
        temp_t = theta_arr(idx)
        temp_tau = tau_arr(idx)
        
        gamma_arr = temp_g
        theta_arr = temp_t
        tau_arr = temp_tau
        
        deallocate(temp_g, temp_t, temp_tau)
    end subroutine rearrange_arrays

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
        
        ! Temporarily add technology
        temp_gamma = Gamma + gamma_val
        temp_theta = Theta + theta_val
        temp_delta = Delta + theta_val * gamma_val
        temp_phi = Phi + gamma_val * gamma_val
        temp_psi = Psi + theta_val * theta_val
        
        Chi = temp_psi / (1.0_dp + temp_theta)
        Sigma = temp_psi * temp_gamma / (1.0_dp + temp_theta) - temp_delta
        Omega = temp_phi - temp_gamma * temp_delta / (1.0_dp + temp_theta)
        
        ! First check if denominator would be positive
        if ((1.0_dp - q*Chi) <= 0.0_dp) return
        
        ! Compute derived quantities
        den = (1.0_dp + temp_theta) * (1.0_dp - q*Chi)
        A = d * temp_delta * N0_over_E0 / den
        B = d * (q * temp_delta * Sigma / den - Omega)
        
        ! Smart initial guess for root finding
        if (has_cache) then
            E_new = cached_E
        else
            E_new = E
        end if
        
        call find_zero_optimized(E_new, tol, max_iter, E_new, converged_E, A, B, inv_E0)
        
        if (.not. converged_E) return
        
        ! Cache the solution
        cached_E = E_new
        has_cache = .true.
        
        N_new = (N0 + q*Sigma*(E0 - E_new)) / (1.0_dp - q*Chi)
        Tau = (N_new/(E0-E_new) + temp_gamma) / (1.0_dp + temp_theta)

        if (tau_val > Tau .or. N_new < 0.0_dp .or. E_new < 0.0_dp) return

        ! Update all aggregate variables
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
        
        ! Iterate until no more technologies can be added
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
                    
                    ! Technology was successfully added
                    num_valid = num_valid + 1
                    valid_list(num_valid) = tech_idx
                    
                    ns = ns + 1
                    N_sol(ns) = N_new
                    E_sol(ns) = E_new
                    tau_sol(ns) = tau_arr(tech_idx)
                    
                    N = N_new; E = E_new
                    added_any = .true.
                else
                    ! Keep in pending list
                    new_pending_count = new_pending_count + 1
                    pending_list(new_pending_count) = tech_idx
                end if
            end do
            
            num_pending = new_pending_count
        end do
    end subroutine smart_reevaluate_pending

    ! Optimized quicksort with better pivot selection
    recursive subroutine quicksort_optimized(arr, idx, first, last)
        integer, intent(inout) :: idx(:)
        real(dp), intent(in) :: arr(:)
        integer, intent(in) :: first, last
        integer :: i, j, pivot_idx, temp_idx
        real(dp) :: pivot
        
        ! Early return if invalid bounds or single element
        if (first >= last) return
        
        if (last - first < 10) then
            ! Use insertion sort for small arrays
            call insertion_sort(arr, idx, first, last)
            return
        end if
        
        ! Better pivot selection (median-of-three)
        pivot_idx = median_of_three_pivot(arr, idx, first, last)
        pivot = arr(idx(pivot_idx))
        
        ! Initialize partition indices within bounds
        i = first
        j = last
        
        do
            ! Ensure indices stay within bounds
            do while (i <= j .and. i <= last .and. arr(idx(i)) < pivot)
                i = i + 1
            end do
            do while (j >= i .and. j >= first .and. arr(idx(j)) > pivot)
                j = j - 1
            end do
            
            if (i >= j) exit
            
            ! Swap only if both indices are valid
            if (i <= last .and. j >= first) then
                temp_idx = idx(i)
                idx(i) = idx(j)
                idx(j) = temp_idx
                i = i + 1
                j = j - 1
            end if
        end do
        
        ! Recursive calls only if valid ranges exist
        if (first < j) call quicksort_optimized(arr, idx, first, j)
        if (i < last) call quicksort_optimized(arr, idx, i, last)
    end subroutine quicksort_optimized

    integer function median_of_three_pivot(arr, idx, first, last) result(pivot_idx)
        real(dp), intent(in) :: arr(:)
        integer, intent(in) :: idx(:), first, last
        integer :: mid
        
        ! Ensure mid point is valid
        mid = first + (last - first) / 2
        
        ! Safety check for valid array bounds
        if (first < 1 .or. last > size(idx)) then
            pivot_idx = first  ! Default to first if bounds invalid
            return
        end if
        
        ! Compare elements safely
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

    subroutine insertion_sort(arr, idx, first, last)
        real(dp), intent(in) :: arr(:)
        integer, intent(inout) :: idx(:)
        integer, intent(in) :: first, last
        integer :: i, j, key
        
        do i = first + 1, last
            key = idx(i)
            j = i - 1
            ! Add bounds check to prevent j from going below first
            do while (j >= first .and. arr(idx(j)) > arr(key))
                idx(j + 1) = idx(j)
                j = j - 1
                ! Add explicit exit condition when j would go below first
                if (j < first) exit
            end do
            ! Ensure j+1 is valid before assignment
            if (j < first - 1) then
                idx(first) = key
            else
                idx(j + 1) = key
            end if
        end do
    end subroutine insertion_sort

    ! Optimized root finding
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
            
            ! Inline residual calculation for speed
            fx = E_residual_inline(x, A, B, inv_E0)
            
            if (abs(fx) < tol) then
                converged = .true.
                exit
            end if
            
            ! Numerical derivative
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

    ! Inlined residual function
    real(dp) function E_residual_inline(E_trial, A, B, inv_E0) result(res)
        real(dp), intent(in) :: E_trial, A, B, inv_E0
        real(dp) :: e_ratio
        
        e_ratio = E_trial * inv_E0
        res = e_ratio*(1.0_dp + A + B) - B*e_ratio*e_ratio - 1.0_dp
    end function E_residual_inline

    ! Optimized file I/O
    subroutine save_results_optimized(ns, N_sol, E_sol, tau_sol, gamma_arr, theta_arr, tau_arr, &
                                    valid_idx, num_valid, seed, E0, N0, q, d, timestamp)
        integer, intent(in) :: ns, num_valid
        real(dp), intent(in) :: N_sol(:), E_sol(:), tau_sol(:)
        real(dp), intent(in) :: gamma_arr(:), theta_arr(:), tau_arr(:)
        integer, intent(in) :: valid_idx(:), seed(:)
        real(dp), intent(in) :: E0, N0, q, d
        character(len=*), intent(in) :: timestamp
        
        character(len=100) :: filename
        integer :: unit, i
        
        ! Save simulation results
        write(filename, '(A,A,A)') 'sim_results_opt_', trim(timestamp), '.dat'
        open(newunit=unit, file=filename, status='replace')
        write(unit, '(A,4F12.5)') '# Parameters: E0, N0, q, d = ', E0, N0, q, d
        write(unit, '(A,I0)') '# Total valid technologies: ', ns
        write(unit, '(A)') '# Format: tech_number, N, E, tau'
        do i = 1, ns
            write(unit, '(I8,3E16.8)') i, N_sol(i), E_sol(i), tau_sol(i)
        end do
        close(unit)
        
        ! Save technology data (only valid ones to save space)
        write(filename, '(A,A,A)') 'valid_tech_opt_', trim(timestamp), '.dat'
        open(newunit=unit, file=filename, status='replace')
        write(unit, '(A)') '# Valid technologies: gamma, theta, tau'
        do i = 1, num_valid
            write(unit, '(3E16.8)') gamma_arr(valid_idx(i)), theta_arr(valid_idx(i)), tau_arr(valid_idx(i))
        end do
        close(unit)
        
        print *, "Results saved with prefix 'opt_", trim(timestamp), "'"
    end subroutine save_results_optimized

    ! Simplified I/O functions
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

end program tech_evo_optimized