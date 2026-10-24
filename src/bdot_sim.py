
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
from pathlib import Path

# ==========================================
# 1. EARTH AND ORBIT PARAMETERS
# ==========================================
mu = 3.986e14
Re = 6.378e6
h = 500e3

r = Re + h
n = np.sqrt(mu / r**3)
inc = np.radians(97.4)

# ==========================================
# 2. CUBESAT PARAMETERS
# ==========================================
I = np.diag([0.042, 0.042, 0.0067])
w0 = np.radians([20, 0, 0])

# B-dot controller
k = 5e4
mmax = 0.2  # A m^2

# ==========================================
# 3. EARTH MAGNETIC FIELD MODEL
# ==========================================
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


# ==========================================
# 4. QUATERNION ATTITUDE EQUATION
# ==========================================
def qdot(q, w):
    x, y, z, s = q
    wx, wy, wz = w

    return 0.5 * np.array([
        s * wx + y * wz - z * wy,
        s * wy + z * wx - x * wz,
        s * wz + x * wy - y * wx,
        -x * wx - y * wy - z * wz
    ])


# ==========================================
# 5. CUBESAT ROTATIONAL DYNAMICS
# ==========================================
def f(t, X):
    q = X[:4] / np.linalg.norm(X[:4])
    w = X[4:]

    # Magnetic field in body coordinates
    Bb = Rotation.from_quat(q).as_matrix().T @ B_inertial(t)

    # Simplified B-dot estimate
    # Orbital variation of the field is neglected
    Bdot = -np.cross(w, Bb)

    # Magnetic dipole command
    m_cmd = -k * Bdot

    # Saturate the total dipole magnitude
    m_norm = np.linalg.norm(m_cmd)
    m = m_cmd * min(1.0, mmax / max(m_norm, 1e-15))

    # Magnetic torque
    tau = np.cross(m, Bb)

    # Euler rotational equations
    wdot = np.linalg.solve(
        I, tau - np.cross(w, I @ w)
    )

    return np.concatenate([qdot(q, w), wdot])


# ==========================================
# 6. RUN THE SIMULATION
# ==========================================
X0 = np.concatenate([
    [0.0, 0.0, 0.0, 1.0],
    w0
])

T = 3 * 2 * np.pi / n

sol = solve_ivp(
    f,
    [0, T],
    X0,
    max_step=2.0,
    rtol=1e-8,
    atol=1e-10
)

if not sol.success:
    raise RuntimeError(sol.message)

# Angular rate in degrees per second
rate = np.degrees(
    np.linalg.norm(sol.y[4:7], axis=0)
)

time_min = sol.t / 60

# ==========================================
# 7. CALCULATE DETUMBLING TIME
# ==========================================
threshold = 0.5

indices = np.where(rate < threshold)[0]

if len(indices) > 0:
    idx = indices[0]

    # Interpolate between adjacent samples
    if idx > 0:
        t1, t2 = time_min[idx - 1], time_min[idx]
        r1, r2 = rate[idx - 1], rate[idx]

        if r1 != r2:
            detumble_time = t1 + (
                (threshold - r1) * (t2 - t1) / (r2 - r1)
            )
        else:
            detumble_time = t2
    else:
        detumble_time = time_min[0]

    print(f"Time below 0.5 deg/s: {detumble_time:.2f} minutes")
else:
    detumble_time = None
    print("0.5 deg/s threshold not reached.")

# ==========================================
# 8. SAVE THE ANGULAR-RATE PLOT
# ==========================================
results_dir = Path("results")
results_dir.mkdir(exist_ok=True)

plt.figure(figsize=(9, 5))
plt.plot(time_min, rate, label=f"Gain k = {k:.0e}")
plt.axhline(
    threshold, color="red", linestyle="--",
    label="Detumbling threshold: 0.5 deg/s"
)

plt.xlabel("Time (min)")
plt.ylabel("Angular Rate (deg/s)")
plt.title("CubeSat B-dot Detumbling")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(results_dir / "rate.png", dpi=200)
plt.show()

# ==========================================
# 9. PRINT SUMMARY
# ==========================================
print("\n--- SIMULATION SUMMARY ---")
print(f"Initial rate: {rate[0]:.3f} deg/s")
print(f"Final rate: {rate[-1]:.6f} deg/s")
print(f"Simulation duration: {T / 60:.2f} minutes")
print(f"Controller gain: {k:.0f}")
print(f"Maximum magnetic dipole: {mmax} A m^2")
print(f"Results saved to: {results_dir / 'rate.png'}")
