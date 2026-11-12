import matplotlib.pyplot as plt
import numpy as np

import functions as f

SEED=12345  # fixed seed -> reproducible results
rng=np.random.default_rng(SEED)

k=10.0 #do not change it
N=50000
dt=0.1
T=20
steps=f.n_steps(T,dt)



r_hist=[]
psi_hist=[]
frequencies = [
    ("Lorentz", 1.0),
    ("Lorentz", 2.0),
    ("Arcsin", None)
]

for freq, w in frequencies:
    theta0,omega=f.initialize(N,w,0,frequency=freq,rng=rng) # use 0 for a disordered initial phase r(t=0)~0
    r_aux, psi_aux, theta=f.evolution_parameter(theta0,omega,k,dt,steps)
    r_hist.append(r_aux)
    psi_hist.append(psi_aux)

times = np.arange(steps) * dt
lorentz_idx = [i for i, (freq, w) in enumerate(frequencies) if freq=="Lorentz"]
arcsin_idx = [i for i, (freq, w) in enumerate(frequencies) if freq=="Arcsin"]
# ============================================================
# Lorentz saturation
# ============================================================
colors = ["blue", "red", "purple", "orange"]
fig, ax = plt.subplots()
for n, i in enumerate(lorentz_idx):
    w = frequencies[i][1]
    color = colors[n % len(colors)]
    r_sat = np.sqrt(max(0.0, 1 - (2 * w) / k))  # Theoretical saturation r(K) with Kc=2Δ (r=0 if K<=Kc)
    ax.plot(times, r_hist[i], color=color, label=rf'Simulation for $\Delta$ = {w:g}')
    ax.axhline(r_sat, linestyle="--", color=color, alpha=0.6, label=rf'Saturation for $\Delta$ = {w:g}')

ax.set_xlabel("Time")
ax.set_xlim(left=0)
ax.set_ylabel("r(t)")
ax.set_ylim(0, 1.0)
ax.set_title(f"Time evolution of the parameter r for K={k}")
ax.legend()
fig.tight_layout()
plt.show()
# ============================================================
# Comparison: Lorentz (every width of the list) vs Arcsin
# ============================================================

fig, ax = plt.subplots()

for n, i in enumerate(lorentz_idx):
    w = frequencies[i][1]
    ax.plot(times, r_hist[i], color=colors[n % len(colors)], label=rf'Lorentz $\Delta$ = {w:g}')
for i in arcsin_idx:
    ax.plot(times, r_hist[i], color="green", label='Arcsin(-1,1)')

ax.set_xlabel("Time")
ax.set_xlim(left=0)
ax.set_ylabel("r(t)")
ax.set_ylim(0, 1.05)
ax.set_title(f"Comparison of frequency distributions for K={k}")
ax.legend()
fig.tight_layout()
plt.show()
