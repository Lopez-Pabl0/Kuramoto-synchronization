import matplotlib.pyplot as plt
import numpy as np

import functions as f

SEED=12345  # fixed seed -> reproducible results
rng=np.random.default_rng(SEED)

N=50000
dt=0.1
T=100
steps=f.n_steps(T,dt)


width=1.0   # Lorentz width Δ (not used for "Arcsin")
freq="Arcsin"

if freq=="Lorentz":
    # K values scaled with the width (Kc=2*width), so the sweep covers the same range of K/Kc for any width
    k_vector = width*np.array([0, 0.5, 1, 1.5, 1.8, 2, 2.2, 2.5, 3, 4, 5])
    t_discard=40.0
    theta0,omega=f.initialize(N,width,0,frequency=freq,rng=rng)
    r_mean, r_desv,r_hist,theta=f.evolution_variable_k(theta0,omega,k_vector,dt,steps,t_discard)

    #Representation of Time Evolution
    times = np.arange(steps) * dt
    fig, ax = plt.subplots()
    for i in range(0,len(k_vector),2):
        ax.plot(times, r_hist[i, :], label=f'K = {k_vector[i]:g}')

    ax.set_xlabel("Time")
    ax.set_xlim(left=0)
    ax.set_ylabel("r(t)")
    ax.set_title(rf"Evolution of r for different K with $\Delta$={width:g}")
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()
    plt.show()

    #Representation of Relation r-K
    Kc = 2 * width
    K_cont = np.linspace(k_vector.min(), k_vector.max(), 500)
    r_theo_cont = np.zeros_like(K_cont)  # Calculate r_theo for K continuous
    mask = K_cont > Kc #Create a vector of True or False
    r_theo_cont[mask] = np.sqrt(1 - Kc / K_cont[mask]) #Now apply the formula only where K>Kc

    fig, ax = plt.subplots()
    ax.plot(K_cont, r_theo_cont, 'k--', label=r"$r_{\mathrm{theo}}=\sqrt{1-\frac{K_c}{K}}$")
    ax.errorbar(k_vector, r_mean, yerr=r_desv, fmt='o', capsize=4, label=r"$r_{\mathrm{sim}}$")
    ax.axvline(Kc, color='gray', ls=':', label=fr"$K_c={Kc}$")

    ax.set_xlabel("K")
    ax.set_xlim(left=-0.01)
    ax.set_ylabel("r")
    ax.set_title(r"Simulation of $r(K)$ and $r_{\mathrm{theo}}(K)$")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    plt.show()

elif freq=="Arcsin":
     k_vector=[0,0.5,1,1.25,1.5,1.75,2,2.25,2.5,3,4]
     t_discard=70.0
     # Forward branch: start from the incoherent state (r~0) and sweep K upwards.
     # Each K starts from the final phases of the previous one (adiabatic sweep).
     theta0,omega=f.initialize(N,None,0,frequency=freq,rng=rng)
     r_mean_up, r_desv_up,r_hist_up,theta_end=f.evolution_variable_k(theta0,omega,k_vector,dt,steps,t_discard,continuation=True)
     # Backward branch: SAME frequencies omega, starting from the final (synchronized) phases of the
     # forward branch at the largest K, so the cycle is not broken.
     k_down = k_vector[::-1]
     r_mean_down, r_desv_down,r_hist_down,theta_end=f.evolution_variable_k(theta_end,omega,k_down,dt,steps,t_discard,continuation=True)

     #Representation of Time Evolution
     times = np.arange(steps) * dt
     fig, ax = plt.subplots()
     for i in range(0,len(k_vector),2):
         ax.plot(times, r_hist_up[i, :], label=f'K = {k_vector[i]}')

     ax.set_xlabel("Time")
     ax.set_xlim(left=0)
     ax.set_ylabel("r(t)")
     ax.set_title("Evolution of r for different K--Forward branch")
     ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
     fig.tight_layout()
     plt.show()

     fig, ax = plt.subplots()
     for i in range(0,len(k_vector),2):
         ax.plot(times, r_hist_down[i, :], label=f'K = {k_down[i]}')

     ax.set_xlabel("Time")
     ax.set_xlim(left=0)
     ax.set_ylabel("r(t)")
     ax.set_title("Evolution of r for different K--Backward branch")
     ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
     fig.tight_layout()
     plt.show()

     fig, ax = plt.subplots()
     ax.errorbar(k_vector, r_mean_up, yerr=r_desv_up, color='blue', alpha=0.6, fmt='o', capsize=4, label=r"$r_{\mathrm{sim}}$ Forward branch")
     ax.errorbar(k_down, r_mean_down, yerr=r_desv_down, color='red', alpha=0.6, fmt='o', capsize=4, label=r"$r_{\mathrm{sim}}$ Backward branch")
     ax.axvline(2.0, color='gray', ls=':', label=fr"$K_c={2:.1f}$")

     ax.set_xlabel("K")
     ax.set_xlim(left=-0.01)
     ax.set_ylabel("r")
     ax.set_title('Hysteresis Study of r as K Varies for the arcsin(-1,1) Distribution')
     ax.legend()
     ax.grid(True)
     fig.tight_layout()
     plt.show()


else:
     raise ValueError(f"{freq} must be an eligible frequency ")
