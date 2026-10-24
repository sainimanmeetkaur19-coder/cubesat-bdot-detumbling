
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation
from pathlib import Path
import csv

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
threshold = 0.5

# Three gains to compare
gains = [1e4, 5e4, 2e5]

# Simulation duration: three orbits
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


def dynamics(t, X, k):
    q = X[:4] / np.linalg.norm(X[:4])
    w = X[4:]

    Bb = Rotation.from_quat(q).as_matrix().T @ B_inertial(t)

    # Simplified B-dot estimate
    Bdot = -np.cross(w, Bb)

    m_cmd = -k * Bdot
    m_norm = np.linalg.norm(m_cmd)
    m = m_cmd * min(1.0, mmax / max(m_norm, 1e-15))

    tau = np.cross(m, Bb)
    wdot = np.linalg.solve(
        I, tau - np.cross(w, I @ w)
    )

    return np.concatenate([qdot(q, w), wdot])


X0 = np.concatenate([[0.0, 0.0, 0.0, 1.0], w0])

plt.figure(figsize=(10, 6))
table = []

for k in gains:
    sol = solve_ivp(
        lambda t, X: dynamics(t, X, k),
        [0, T],
        X0,
        max_step=2.0,
        rtol=1e-8,
        atol=1e-10
    )

    if not sol.success:
        raise RuntimeError(sol.message)

    time_min = sol.t / 60
    rate = np.degrees(np.linalg.norm(sol.y[4:7], axis=0))

    # First sampled time below threshold
    indices = np.where(rate < threshold)[0]

    if len(indices):
        idx = indices[0]

        if idx > 0 and rate[idx - 1] != rate[idx]:
            t1, t2 = time_min[idx - 1], time_min[idx]
            r1, r2 = rate[idx - 1], rate[idx]
            detumble_time = t1 + (
                (threshold - r1) * (t2 - t1) / (r2 - r1)
            )
        else:
            detumble_time = time_min[idx]
        result = round(detumble_time, 2)
    else:
        result = None

    table.append([k, result, round(rate[-1], 6)])
    plt.plot(time_min, rate, label=f"k = {k:.0f}")

    print(
        f"Gain {k:.0f}: "
        f"detumbling time = "
        f"{result if result is not None else 'not reached'} min"
    )

plt.axhline(
    threshold, color="black", linestyle="--",
    label="Threshold = 0.5 deg/s"
)
plt.xlabel("Time (min)")
plt.ylabel("Angular Rate (deg/s)")
plt.title("B-dot Controller Gain Comparison")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(results_dir / "gain_comparison.png", dpi=200)
plt.show()

# Save comparison table
with open(
    results_dir / "detumbling_times.csv",
    "w",
    newline=""
) as f:
    writer = csv.writer(f)
    writer.writerow([
        "Gain",
        "Detumbling Time (min)",
        "Final Rate (deg/s)"
    ])
    writer.writerows(table)

print("\nComparison table saved to results/detumbling_times.csv")
print("Plot saved to results/gain_comparison.png")
