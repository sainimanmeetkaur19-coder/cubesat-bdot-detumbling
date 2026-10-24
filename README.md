# CubeSat B-Dot Detumbling and Hohmann Transfer Orbit

A CubeSat attitude-control and orbital-maneuver simulation project using **Python** and **GMAT (General Mission Analysis Tool)**. The project studies magnetic detumbling using a B-dot controller and a Hohmann transfer from Low Earth Orbit (LEO) to a geostationary-altitude orbit.

## Project Objectives

- Simulate CubeSat angular-rate reduction using a B-dot detumbling controller.
- Compare controller gains and their effects on detumbling time.
- Analyse magnetorquer dipole-moment commands and saturation.
- Calculate the velocity increments required for a Hohmann transfer.
- Validate the transfer trajectory using GMAT and compare it with analytical calculations.

## Technologies Used

- **Python:** NumPy, SciPy, Matplotlib
- **Orbital simulation:** GMAT
- **Version control:** Git and GitHub

## 1. B-Dot Detumbling Simulation

The B-dot controller uses changes in the measured magnetic field to generate magnetic dipole commands. The interaction between the magnetic dipole and the Earth's magnetic field produces torque that reduces the spacecraft's angular velocity.

### Simulation Parameters

| Parameter | Value |
|---|---:|
| Initial orbit altitude | 500 km |
| Orbital inclination used in the model | 97.4° |
| Moment of inertia | diag(0.042, 0.042, 0.0067) kg·m² |
| Maximum magnetic dipole moment | 0.2 A·m² |
| Nominal controller gain | 50,000 |
| Detumbling threshold | 0.5°/s |

### Results

For the nominal initial angular velocity of [10°, 10°, 10°]/s, the simulation reached the detumbling threshold of 0.5°/s in approximately **57.91 minutes**.

An alternative initial tumble case, [20°, 0°, 0°]/s, reached the same threshold in approximately **85.47 minutes**.

### Controller Gain Comparison

| Controller gain | Detumbling time |
|---:|---:|
| 10,000 | 263.4 min |
| 50,000 | 57.91 min |
| 200,000 | 51.02 min |

The results show that increasing the controller gain reduced the simulated detumbling time for the tested cases. However, the improvement was smaller at the highest gain, where magnetorquer saturation and the simplified model can affect performance.

### Generated Results

- `results/rate.png` — Angular-rate reduction over time.
- `results/gain_comparison.png` — Detumbling performance for different controller gains.
- `results/dipole_magnitude.png` — Commanded magnetic dipole magnitude and saturation.
- `results/detumbling_times.csv` — Numerical comparison of detumbling times.

The dipole analysis showed a maximum commanded magnitude of 0.2 A·m², with approximately 13.5% of the evaluated samples at saturation.

## 2. Hohmann Transfer Orbit

A Hohmann transfer was analysed between an initial 500 km circular orbit and a final orbit with an apogee altitude of 35,786 km.

### Analytical Parameters

| Parameter | Value |
|---|---:|
| Initial orbit altitude | 500 km |
| Initial orbital radius | 6,878 km |
| Target apogee altitude | 35,786 km |
| Target orbital radius | 42,164 km |
| Transfer semi-major axis | 24,521 km |
| First-burn magnitude (analytical) | 2.3698 km/s |
| Second-burn magnitude (analytical) | 1.4463 km/s |
| Total ideal velocity increment | 3.8161 km/s |
| Ideal transfer time | Approximately 5.31 hours |

### GMAT Verification

The transfer was modelled in GMAT using an impulsive maneuver and propagation to the transfer orbit's apogee.

| Parameter | Result |
|---|---:|
| Target apogee altitude | 35,786 km |
| Final GMAT first-burn value | 2.3672 km/s |
| Reported GMAT apogee altitude | 35,786.52 km |
| Absolute altitude difference | 0.52 km |
| Relative altitude difference | Approximately 0.00145% |

The GMAT-reported apogee is close to the target altitude. The GMAT burn value differs from the ideal analytical value, so the comparison should be interpreted in the context of the actual initial orbital elements, maneuver frame, and simulation settings. The burn value and achieved apogee are simulation results, not flight-validated measurements.

The GMAT report is saved at `results/hohmann_transfer_report.txt`, and the mission script is saved at `gmat/hohmann_gto.script`.

## Repository Structure

```text
cubesat-bdot-detumbling/
├── README.md
├── LICENSE
├── .gitignore
├── src/
│   ├── bdot_sim.py
│   ├── gain_comparison.py
│   ├── dipole_analysis.py
│   └── hohmann_transfer.py
├── gmat/
│   └── hohmann_gto.script
├── results/
│   ├── rate.png
│   ├── gain_comparison.png
│   ├── dipole_magnitude.png
│   ├── detumbling_times.csv
│   └── hohmann_transfer_report.txt
└── docs/
    └── CubeSat_BDot_Detumbling_and_Hohmann_Transfer_Report.pdf
```

*The tree above describes the intended project structure; include only files that actually exist in the repository.*

## How to Run

### Python simulations

1. Install Python 3 and the required libraries:

   ```bash
   python -m pip install numpy scipy matplotlib
   ```

2. From the repository root, run the scripts:

   ```bash
   python src/bdot_sim.py
   python src/gain_comparison.py
   python src/dipole_analysis.py
   python src/hohmann_transfer.py
   ```

3. Review the generated plots and numerical outputs in `results/`.

### GMAT simulation

1. Open GMAT.
2. Load `gmat/hohmann_gto.script`.
3. Confirm the spacecraft orbital elements, maneuver settings, and propagation stop conditions.
4. Run the mission and compare the resulting apogee with the target altitude.

## Model Limitations

The detumbling simulation uses a simplified Earth magnetic-dipole model and neglects the orbital contribution to the body-frame magnetic-field derivative used by the B-dot controller. It assumes ideal measurements and actuators, with no sensor noise, actuator dynamics, or external disturbance torques. The inertia tensor is an estimate rather than a verified CAD-derived value.

The orbital transfer calculation assumes an ideal two-impulse Hohmann transfer. GMAT results depend on the initial orbital elements, coordinate and maneuver frames, propagation settings, and burn configuration. Therefore, the results are preliminary simulations and should not be interpreted as validated flight predictions.

## Conclusion

This project combines spacecraft attitude-control simulation with orbital-transfer analysis. The Python models demonstrate B-dot detumbling behaviour and controller-gain sensitivity, while GMAT provides an independent numerical environment for evaluating the transfer trajectory. Together, they illustrate key concepts in CubeSat attitude dynamics, magnetic actuation, orbital mechanics, and mission analysis.
