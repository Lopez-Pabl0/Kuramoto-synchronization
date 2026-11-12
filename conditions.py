import matplotlib.pyplot as plt
import numpy as np

import functions as f

SEED=12345  # fixed seed -> reproducible results
rng=np.random.default_rng(SEED)

width=1.0
k=5.0   # To test the parameters dt & N it's better to have a K>K_c

dt=0.1
T=50
steps=f.n_steps(T,dt)


N=[100,1000,5000,50000]

r_hist=np.zeros((len(N), steps))
psi_hist=np.zeros((len(N), steps))

for i,n in enumerate(N):
    theta0, omega=f.initialize(n,width,initial_phase=0, frequency="Lorentz",rng=rng)
    r_hist[i, :], psi_hist[i, :], theta=f.evolution_parameter(theta0,omega,k,dt,steps)

#Representation
times = np.arange(steps) * dt
fig, ax = plt.subplots()
for i, n in enumerate(N):
    ax.plot(times, r_hist[i, :], label=f'N = {n}')

ax.set_xlabel("Time")
ax.set_xlim(left=0)
ax.set_ylabel("r(t)")
ax.set_title(f"Time evolution of the parameter r for K={k}")
ax.legend()

fig.tight_layout()
plt.show()

#Now we can try to change the value of dt (this is much slower be patient)
N=50000

T=50

dt_values=[0.5,0.1,0.01,0.001]
#now the _hist has a variable length on columns so its better to use lists, also we are not going to make operations just representation
r_hist = []
psi_hist = []

# Same initial phases and frequencies for every dt, so the only difference between curves is the time step
theta0, omega=f.initialize(N,width,initial_phase=0, frequency="Lorentz",rng=rng)
for dt_i in dt_values:
    steps=f.n_steps(T,dt_i)
    r_aux, psi_aux, theta=f.evolution_parameter(theta0,omega,k,dt_i,steps)
    r_hist.append(r_aux)
    psi_hist.append(psi_aux)

#Representation
fig, ax = plt.subplots()
for i, dt_i in enumerate(dt_values):
    times = np.arange(len(r_hist[i])) * dt_i
    ax.plot(times, r_hist[i], label=f'dt = {dt_i}')

ax.set_xlabel("Time")
ax.set_xlim(left=0)
ax.set_ylabel("r(t)")
ax.set_title(f"Time evolution of the parameter r for K={k} and N={N}")
ax.legend()

fig.tight_layout()
plt.show()
