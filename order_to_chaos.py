import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j0

import functions as f

SEED=12345  # fixed seed -> reproducible results
rng=np.random.default_rng(SEED)

k=0.0 #do not change it
N=50000
dt=0.1
T=20
steps=f.n_steps(T,dt)



freq="Lorentz"
width=[0.1,0.5,1.0,2.0]

times = np.arange(steps) * dt
t = np.linspace(0, T, 1000)  # continuous time for the theoretical curve

if freq=="Lorentz":
    r_hist=np.zeros((len(width), steps))
    psi_hist=np.zeros((len(width), steps))
    for i,w in enumerate(width):
        theta0,omega=f.initialize(N,w,1,frequency=freq,rng=rng) # use 1 for an order initial phase r(t=0)=1
        r_hist[i, :], psi_hist[i, :], theta=f.evolution_parameter(theta0,omega,k,dt,steps)
    #Representation
    fig, ax = plt.subplots()
    for i, w in enumerate(width):
            color = f"C{i}"
            ax.plot(times, r_hist[i, :], color=color, label=rf'Simulation for $\Delta$ = {w:g}')
            ax.plot(t, np.exp(-w*t), '--', color=color, alpha=0.6,
                    label=rf'Theoretical Curve $e^{{-\Delta t}}$ for $\Delta = {w:g}$')
    ax.set_xlabel("Time")
    ax.set_xlim(left=0)
    ax.set_ylabel("r(t)")
    ax.set_ylim(0,1.0)
    ax.set_title(f"Time evolution of the parameter r for K={k}")
    ax.legend()
    fig.tight_layout()
    plt.show()

elif freq=="Arcsin":
    theta0,omega=f.initialize(N,None,1,frequency=freq,rng=rng)
    r_hist, psi_hist, theta=f.evolution_parameter(theta0,omega,k,dt,steps)
    #Representation
    r_theo = np.abs(j0(t))
    fig, ax = plt.subplots()
    ax.plot(times, r_hist, color='#5dade2', alpha=0.8,label='Simulation')
    ax.plot(t, r_theo, '--',color='blue', alpha=0.6, label=r'Theoretical Curve $|J_0(t)|$')

    ax.set_xlabel("Time")
    ax.set_xlim(left=0)
    ax.set_ylabel("r(t)")
    ax.set_ylim(0,1.0)
    ax.set_title(f"Time evolution of the parameter r for K={k}")
    ax.legend()
    fig.tight_layout()
    plt.show()
else:
     raise ValueError(f"{freq} must be an eligible frequency ")
