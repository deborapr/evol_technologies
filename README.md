# evol_technologies

Codes for simulations of the study _"Exploring the conditions for sustainability with open-ended innovation"_ 

In this work, we introduce a model of genuine open-ended technological innovation, in which each innovation event corresponds to a random draw from the space of possible technologies, characterized by parameters that determine their environmental impact and interactions with the population. Innovation is endogenous: a new technology can successfully invade only if it satisfies viability conditions that depend on the current state of both the environment and the population. Further details can be found in the main text.

**Algorithm for numerical simulations**

Time is measured in terms of innovations; that is, at time $t$, the number of technologies introduced is $t$. To convert to real time, one must draw random time intervals from an exponential distribution with mean $(\lambda N)^{-1}$ and sum them cumulatively.

For each realization, the dynamics starts from the state $t=0$, with $E(t=0)=E_0$ and $N(t=0)=N_0(E_0)$, and fixed model parameters $q$ and $d$. 
Starting from a state with no active technologies ($n_s = 0$), we attempt to activate technologies sequentially.
At each step, we compute the updated equilibrium values of the environment and population, and check if the next technology satisfies the activation condition, i.e., whether it can successfully invade. The process stops as soon as this condition is no longer met.

Computationally, the procedure is as follows:

1. Generate an ensemble of $T$ technologies $(\gamma_t, \theta_t)$, sampled independently from exponential distributions with means $\bar{\gamma}$ and $\bar{\theta}$, respectively.
2. Compute $\tau_t = \gamma_t / \theta_t$, and sort the technologies in increasing order of $\tau_t$.
3. Loop through the sorted list, attempting to activate one technology at a time:
    1. **(*)** Increment the counter: $n_s \leftarrow n_s + 1$, and consider the technology  $\tau_{n_s}$. Compute the accumulated variables $\Theta$, $\Gamma$, $\Delta$, $\Psi$, and $\Phi$, the derived quantities $\chi$, $\Omega$, and $\Sigma$, and determine the corresponding equilibrium values of $E$ and $N$;
    2. Compute the updated threshold $\tau$, and check whether the next technology satisfies the activation condition given by Eq.~\ref{eq:tau_theta_gamma}: $\tau_{n_s + 1} \leq \tau$;
    3. If the condition is violated, stop the process and proceed to the next realization. Store the current values of the environment, population, and the last successfully activated $\tau_t$. Otherwise, return to step **(*)**.

The algorithm yields, for each realization, the final number of active technologies $n_s$, the equilibrium values of $E$ and $N$, and the highest accepted value of $\tau_t$. We performed 100 realizations for each set of parameters. 

