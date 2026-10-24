
import numpy as np

# Earth constants (km-based units)
mu = 398600.4418       # km^3/s^2
Re = 6378.0            # km

# Orbit altitudes
h_initial = 500.0      # km
h_apogee = 35786.0     # km

# Orbital radii from Earth's centre
r1 = Re + h_initial
r2 = Re + h_apogee

# Semi-major axis of transfer ellipse
a_transfer = (r1 + r2) / 2

# Initial circular orbital velocity
v_circular = np.sqrt(mu / r1)

# Transfer velocity at perigee
v_transfer_perigee = np.sqrt(
    mu * (2 / r1 - 1 / a_transfer)
)

# First burn: enter transfer orbit
delta_v1 = v_transfer_perigee - v_circular

# Transfer velocity at apogee
v_transfer_apogee = np.sqrt(
    mu * (2 / r2 - 1 / a_transfer)
)

# Circular orbital velocity at target radius
v_circular_final = np.sqrt(mu / r2)

# Second burn: circularize at apogee (optional)
delta_v2 = v_circular_final - v_transfer_apogee

# Transfer time: half the orbital period of the ellipse
transfer_time = np.pi * np.sqrt(a_transfer**3 / mu)

print("=== HOHMANN TRANSFER RESULTS ===")
print(f"Initial orbit radius: {r1:.2f} km")
print(f"Target apogee radius: {r2:.2f} km")
print(f"Transfer semi-major axis: {a_transfer:.2f} km")
print(f"Initial circular velocity: {v_circular:.4f} km/s")
print(f"Transfer perigee velocity: {v_transfer_perigee:.4f} km/s")
print(f"First burn delta-v: {delta_v1:.4f} km/s")
print(f"Transfer apogee velocity: {v_transfer_apogee:.4f} km/s")
print(f"Final circular velocity: {v_circular_final:.4f} km/s")
print(f"Second burn delta-v: {delta_v2:.4f} km/s")
print(f"Total delta-v for circularization: {delta_v1 + delta_v2:.4f} km/s")
print(f"Transfer time: {transfer_time / 3600:.2f} hours")
