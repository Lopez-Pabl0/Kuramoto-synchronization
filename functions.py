import matplotlib.pyplot as plt  #graphics
import numpy as np  #vectors, matrix, random numbers


def n_steps(T, dt):
    """Number of stored time points so that the evolution covers [0, T] (round avoids float errors such as 0.3/0.1=2.999...)"""
    return round(T/dt) + 1

def initialize (N,width,initial_phase,frequency="",rng=None):
    # rng: np.random.Generator; pass one created with a fixed seed to make the run reproducible
    if rng is None:
        rng = np.random.default_rng()
    #phases (theta)
    if initial_phase==0:
        theta0 = np.pi*rng.uniform(-1, 1, N) # initial phases between -π and π, disorder r~0
    else:
        theta0=np.ones(N)*0.25*np.pi #every initial phase in π/4, order r=1

    #frequencies (omega)
    eps = 1e-12
    x=rng.uniform(0, 1, N) #random number between 0 and 1
    x = np.clip(x, eps, 1-eps) #This avoids possible problems with tg of 0 or 1//Cut the interval to np.clip(array, min_val, max_val)
    if frequency=="Lorentz":
        omega=width*np.tan(np.pi*(x-0.5))  #transformation for Lorentz
    elif frequency=="Arcsin":
        omega=np.sin(np.pi*(x-0.5))   # transformation for Arcsin (width is not used)
    else:
        raise ValueError(f"{frequency} must be an eligible frequency ")
    return theta0, omega

def parameter(theta):
    """Returns (r, psi) from the theta phases set."""
    z = np.mean(np.exp(1j * theta)) #1j=the complex unit; after doing the mean We have a complex scalar
    return np.abs(z), np.angle(z)

def kuramoto_parameter_eq(theta, omega, k):
    """
    dtheta/dt = omega + K*r* sin( Psi - theta_i)
    Orden O(N)  """
    r, psi = parameter(theta)

    return omega + k * r * np.sin(psi - theta)

def evolution_parameter(theta0, omega, k, dt, steps):
    theta = theta0.copy()
    r_hist = np.zeros(steps)
    psi_hist = np.zeros(steps)
    r_hist[0],psi_hist[0]=parameter(theta0)
    for j in range(1,steps) :
        # RK4
        k1 = kuramoto_parameter_eq(theta, omega, k)
        k2 = kuramoto_parameter_eq(theta + 0.5 * dt * k1, omega, k)
        k3 = kuramoto_parameter_eq(theta + 0.5 * dt * k2, omega, k)
        k4 = kuramoto_parameter_eq(theta + dt * k3, omega, k)
        theta=theta + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        theta=np.mod(theta, 2 * np.pi) # np.mod(x, m) returns the remainder of dividing x by m [0,2π)
        r, psi = parameter(theta)
        r_hist[j] = r
        psi_hist[j] = psi
    return r_hist, psi_hist, theta

def evolution_variable_k(theta0, omega, k_vector, dt, steps,t_discard,continuation=False):
    """
    continuation=False: every K starts from theta0.
    continuation=True:  every K starts from the final phases of the previous K
                        (adiabatic sweep, needed to follow a branch in the hysteresis study).
    Returns also the final phases so another sweep can start from them.
    """
    burn_idx = round(t_discard / dt)  # steps to discard (t<t_discard)
    if burn_idx >= steps:
        raise ValueError(f"t_discard={t_discard} must be smaller than the total time T={(steps-1)*dt:g}, "
                         "otherwise there is no stationary zone to average r")
    r_hist=np.zeros((len(k_vector), steps))
    r_mean = np.zeros(len(k_vector))
    r_desv = np.zeros(len(k_vector))
    theta = theta0.copy()
    for i, k in enumerate(k_vector):
        if not continuation:
            theta = theta0.copy()
        r_aux=np.zeros(steps)
        r_aux[0], psi = parameter(theta)
        aux = []
        for j in range(1, steps):
            # RK4
            k1 = kuramoto_parameter_eq(theta, omega, k)
            k2 = kuramoto_parameter_eq(theta + 0.5*dt*k1, omega, k)
            k3 = kuramoto_parameter_eq(theta + 0.5*dt*k2, omega, k)
            k4 = kuramoto_parameter_eq(theta + dt*k3,      omega, k)
            theta = theta + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
            theta = np.mod(theta, 2*np.pi)
            r, psi = parameter(theta)
            r_aux[j]=r
            if j >= burn_idx:        # equivalent to t >= t_discard
                aux.append(r)
        r_hist[i, :]=r_aux
        aux = np.array(aux, dtype=float)
        if aux.size > 0:
            r_mean[i] = np.mean(aux)
            r_desv[i] = np.std(aux)
        else:
            r_mean[i] = np.nan
            r_desv[i] = np.nan
    return r_mean, r_desv, r_hist, theta




def plot_phases_static(theta, r, psi, title=""):
    x, y = np.cos(theta), np.sin(theta)
    fig, ax = plt.subplots(figsize=(5,5))
    circle = plt.Circle((0,0), 1, color='gray', fill=False, linestyle='--')
    ax.add_artist(circle)
    ax.scatter(x, y, s=18)
    ax.arrow(0,0, r*np.cos(psi), r*np.sin(psi), head_width=0.01, head_length=0.01)
    ax.set_aspect('equal'); ax.set_xlim(-1.2,1.2); ax.set_ylim(-1.2,1.2)
    ax.axhline(0, color='gray', linewidth=0.8, linestyle='-') #axis-x
    ax.axvline(0, color='gray', linewidth=0.8, linestyle='-') #axis-y
    ax.set_title(f"{title}  r={r:.3f}")
    fig.tight_layout()
    plt.show()

def plot_frequency(omega, frequency, width=None):
    if frequency=="Lorentz":
        w_min, w_max = -10, 10
        w = np.linspace(w_min, w_max, 800)
        g_w = width/(np.pi*(w**2+width**2)) # g(ω) theoretical
    elif frequency=="Arcsin":
        w_min, w_max = -1, 1
        w = np.linspace(w_min, w_max, 802)[1:-1]  # without ω=±1, where g(ω) diverges
        g_w = 1/(np.pi*np.sqrt(1-w**2))  # g(ω) theoretical
    else:
        raise ValueError(f"{frequency} must be an eligible frequency ")

    # Graph
    plt.figure(figsize=(8,5))
    # Normalize with ALL the frequencies (density=True would only use those inside the range, inflating the histogram for heavy tails)
    bins = 200
    bin_width = (w_max - w_min) / bins
    plt.hist(omega, bins=bins, range=(w_min, w_max), weights=np.full(len(omega), 1/(len(omega)*bin_width)), alpha=0.7, label='Histogram')
    plt.plot(w, g_w, 'k-', linewidth=2, label='g(ω) theoretical')
    plt.xlim(w_min, w_max)
    plt.xlabel('ω')
    plt.ylabel('Probability Density')
    plt.title('Frequencies distribution ω vs g(ω) Theoretical')
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_r_t(r_hist, steps, dt,k):
    times = np.arange(steps) * dt
    fig, ax = plt.subplots()
    ax.plot(times, r_hist, label=f'K = {k:.2f}')
    ax.set_xlabel("Time")
    ax.set_xlim(left=0)
    ax.set_ylabel("r(t)")
    ax.set_ylim(0, 1.0)
    ax.set_title("Time evolution of the parameter r")
    ax.legend()
    fig.tight_layout()
    plt.show()
