
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
from pathlib import Path

# Earth and orbit parameters
mu = 3.986e14
Re = 6.378e6
h = 500e3
r = Re + h
n = np.sqrt(mu / r**3)
inc = np.radians(97.4)

# CubeSat parameters
I = np.diag([0.042, 0.042, 0.0067])
w0 = np.radians([10.0, 10.0, 10.0])
mmax = 0.2

# Choose the controller gain
k = 5e4

T = 3 * 2 * np.pi / n
results_dir = Path("results")
results_dir.mkdir(exist_ok=True)


def B_inertial(t):
    u = n * t
    rh = np.array([
        np.cos(u),
        np.sin(u) * np.cos(inc),
        np.sin(u) * np.sin(inc)
    ])
    z = np.array([0.0, 0.0, 1.0])
    return (7.8e15 / r**3) * (
        z - 3 * np.dot(z, rh) * rh
    )


def qdot(q, w):
    x, y, z, s = q
    wx, wy, wz = w

    return 0.5 * np.array([
        s * wx + y * wz - z * wy,
        s * wy + z * wx - x * wz,
        s * wz + x * wy - y * wx,
        -x * wx - y * wy - z * wz
    ])


def calculate_control(q, w, t):
    q = q / np.linalg.norm(q)
    Bb = Rotation.from_quat(q).as_matrix().T @ B_inertial(t)

    # Simplified B-dot estimate
    Bdot = -np.cross(w, Bb)

    m_cmd = -k * Bdot
    m_norm = np.linalg.norm(m_cmd)

    # Limit total dipole magnitude
    m = m_cmd * min(1.0, mmax / max(m_norm, 1e-15))

    return Bb, Bdot, m


def dynamics(t, X):
    q = X[:4] / np.linalg.norm(X[:4])
    w = X[4:]

    Bb, Bdot, m = calculate_control(q, w, t)
    tau = np.cross(m, Bb)

    wdot = np.linalg.solve(
        I, tau - np.cross(w, I @ w)
    )

    return np.concatenate([qdot(q, w), wdot])


# Initial state
X0 = np.concatenate([[0.0, 0.0, 0.0, 1.0], w0])

sol = solve_ivp(
    dynamics,
    [0, T],
    X0,
    max_step=2.0,
    rtol=1e-8,
    atol=1e-10
)

if not sol.success:
    raise RuntimeError(sol.message)

# Calculate magnetic dipole magnitude at each time
dipole = []

for i, t in enumerate(sol.t):
    q = sol.y[:4, i]
    w = sol.y[4:, i]

    _, _, m = calculate_control(q, w, t)
    dipole.append(np.linalg.norm(m))

dipole = np.array(dipole)
time_min = sol.t / 60

# Determine saturation
saturated = dipole >= mmax * (1 - 1e-9)
sat_percentage = 100 * np.mean(saturated)

print(f"Gain: {k:.0f}")
print(f"Maximum dipole: {dipole.max():.4f} A m^2")
print(f"Samples at saturation: {sat_percentage:.2f}%")

# Plot magnetic dipole magnitude
plt.figure(figsize=(10, 5))
plt.plot(time_min, dipole, label="Dipole magnitude")
plt.axhline(
    mmax, color="red", linestyle="--",
    label="Maximum = 0.2 A m^2"
)
plt.xlabel("Time (min)")
plt.ylabel("Magnetic Dipole Magnitude (A m^2)")
plt.title("Magnetorquer Dipole and Saturation")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(results_dir / "dipole_magnitude.png", dpi=200)
plt.show()

print("Plot saved to results/dipole_magnitude.png")
