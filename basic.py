import numpy as np

import functions as f

SEED=12345  # fixed seed -> reproducible results
rng=np.random.default_rng(SEED)

k=3

T=50
dt=0.1
steps=f.n_steps(T,dt)

N=50000
freq="Lorentz"  # "Lorentz" or "Arcsin"
width=1.0       # Lorentz width Δ (not used for "Arcsin")

theta0,omega=f.initialize(N,width,1,frequency=freq,rng=rng)
r,psi=f.parameter(theta0)
#Show the initial situation
f.plot_phases_static(theta0, r, psi, title="Initial Situation")
f.plot_frequency(omega,freq,width)

#Let the system evolve
r_hist, psi_hist,theta=f.evolution_parameter(theta0, omega, k, dt, steps)
#Show the final situation
f.plot_r_t(r_hist, steps, dt,k)
