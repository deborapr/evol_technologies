# evol_technologies

Codes for simulations of the study _"Exploring the conditions for sustainability with open-ended innovation"_ 

In this work, we introduce a model of genuine open-ended technological innovation, in which each innovation event corresponds to a random draw from the space of possible technologies, characterized by parameters that determine their environmental impact and interactions with the population. Innovation is endogenous: a new technology can successfully invade only if it satisfies viability conditions that depend on the current state of both the environment and the population. Further details can be found in the main text.

**Algorithm for numerical simulations**

Time is measured in terms of innovations, i.e., at time $t$ the number of technologies introduced is $t$. In order to convert to real time, one needs to draw exponential random time intervals from the exponential distribution with mean $(\lambda N)^{-1}$, and sum them up.

The dynamics starts from the state $t=0$, with $E(t=0)=E_0$ and $N(t=0)=N_0(E_0)$. Then, the process proceeds sequentially: a technology is drawn at random with parameters given from exponential probability distributions, and its stability is checked via the condition $\gamma_t/\theta_t \leq \tau$. If condition is not met, the time is incremented, $t=t+1$, and a new technology is drawn. If the condition is satisfied, the new equilibrium values of $E$ and $N$ are computed, along with the updated threshold $\tau$. One must then verify whether the status (active/inactive) of any previously introduced technology has changed, computing all relevant quantities self-consistently. After convergence, time is incremented again and the process continues with the next technology.

Computationally, we implement the following algorithm: 

    1. Build an ensemble of _T_ technologies \gamma_t,\theta_t, with parameters $q$ and $d$ fixed and common to all of them;
    2. Calculate $\tau_t = \gamma_t/\theta_t$, and sort all technologies in increasing order of $\tau_t$;
    3. Looping over the sorted list from $t=1$ to $T$, attempt to activate each technology:
        a. (*)Compute the accumulated variables $\Theta$, $\Gamma$, $\Delta$, $\Psi$, $\Phi$, as well as the derived quantities $\chi$  and $\Sigma$, and determine the corresponding equilibrium values of $E$ and $N$;
        b. For the technology under evaluation, check the condition in Eq.~\ref{eq:tau_theta_gamma}. If the condition is not satisfied, exclude the technology from the current set, recompute the accumulated variables without it, and list it as pending;
        c. If the condition is satisfied, re-evaluate all previously pending technologies in the same order, using the updated system state (repeat from step * for each one).
    4.  Continue this process until the entire ensemble has been evaluated. The algorithm yields the final state with $n_s$ active technologies, the equilibrium values $E$ and $N$, and the highest accepted value of $\tau_t$.
