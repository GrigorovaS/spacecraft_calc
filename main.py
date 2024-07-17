from orbit_calculations import runge_kutta_orbit, plot_earth, N2, h2, position0, velocity0, positions1, T
from hohmann_transfer import calculate_hohmann_transfer, convert_seconds_to_hms, mass, force_of_engine, force_of_engine2
from sgp4.api import Satrec, jday, SatrecArray, accelerated
from find_aqua import tle_line1_old, tle_line2_old, tle_line1_new, tle_line2_new, jd, fr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from tle_transformation import tle_to_keplerian, keplerian_to_cartesian, extract_tle_time

# Problem 1
# Perform orbit calculations
final_position, final_velocity, final_position_mag, final_velocity_mag = runge_kutta_orbit(position0, velocity0, N2, h2, T)

# Problem 2
# Perform Hohmann transfer calculations
burn_duration, total_delta_v, total_delta_v_analytical = calculate_hohmann_transfer(final_position_mag, final_velocity_mag, mass, force_of_engine)

# Perform time format conversion
days, hours, minutes, seconds = convert_seconds_to_hms(burn_duration)

print(f"Total Delta V for Hohmann Transfer: {np.round(total_delta_v)} m/s")
print(f"Total Analytical Delta V for Hohmann Transfer: {np.round(total_delta_v_analytical)} m/s")
print(f"Total Burn Duration with F = 10 mN:{np.round(days, 0)}d {np.round(hours, 0)}h {np.round(minutes, 0)}m {np.round(seconds, 0)}s")

# Perform Hohmann transfer calculations for F = 500 N
burn_duration2, total_delta_v, total_delta_v_analytical = calculate_hohmann_transfer(final_position_mag, final_velocity_mag, mass, force_of_engine2)

# Perform time format conversion
days2, hours2, minutes2, seconds2 = convert_seconds_to_hms(burn_duration2)

print(f"Total Burn Duration with F = 500 N:{np.round(days2, 0)}d {np.round(hours2, 0)}h {np.round(minutes2, 0)}m {np.round(seconds2, 0)}s")

# Problem 3
# Initialize record from the TLE lines
aqua_old = Satrec.twoline2rv(tle_line1_old, tle_line2_old)
aqua_new = Satrec.twoline2rv(tle_line1_new, tle_line2_new)

# Find Aqua's position and velocity from old and new tle
e_old, position_old, velocity_old = aqua_old.sgp4(jd, fr)
e_new, position_new, velocity_new = aqua_new.sgp4(jd, fr)

print(f"Position using old TLE: {position_old}km")
print(f"Position using new TLE: {position_new}km")

# Problem 4
# Parse TLE and calculate Keplerian elements
keplerian_elements = tle_to_keplerian(tle_line1_new, tle_line2_new)

print("Keplerian Elements:")
print("Semi-major axis (m):", keplerian_elements[0])
print("Eccentricity:", keplerian_elements[1])
print("Inclination (rad):", keplerian_elements[2])
print("RAAN (rad):", keplerian_elements[3])
print("Argument of Perigee (rad):", keplerian_elements[4])
print("Mean Anomaly (rad):", keplerian_elements[5])

# Convert to Cartesian coordinates
position, velocity = keplerian_to_cartesian(*keplerian_elements)

print("Position Vector (m):", position)
print("Velocity Vector (m/s):", velocity)


# Perform orbit calculations for one day from the latest TLE epoch time
T2 = 24 * 3600
tle_position, tle_velocity, tle_position_mag, tle_velocity_mag = runge_kutta_orbit(position, velocity, N2, h2, T2)

print(f"Position after one day propagation with the Runge-Kutta propagator using new TLE: {tle_position[-1] / 1000} km")
print(f"Velocity after one day propagation with the Runge-Kutta propagator using new TLE: {tle_velocity[-1] / 1000} km/s")

# Convert TLE time to yy, mm, dd, hh, mins, ss
yy, mm, dd, hh, mins, ss = extract_tle_time(tle_line1_new)

# Calculate day after based on TLE line 1
day_after = (extract_tle_time(tle_line1_new)[2] + 1)

# Convert to julian day
jd0, fr0 = jday(yy, mm, dd, hh, mins, ss)
jd24, fr24 = jday(yy, mm, day_after, hh, mins, ss) 

# Propagate trajectory for one day from the latest TLE epoch time using sgp4
e24, position24, velocity24 = aqua_new.sgp4(jd24, fr24)

print(f"Position after one day propagation with SGP4 using new TLE: {position24}km")
print(f"Velocity after one day propagation with SGP4 using new TLE: {velocity24}km/s")

# Propagete the trajectory over the course of one day from the latest TLE epoch
# time using sgp4
# Create arrays to be used 
total_steps = 10000  # Total steps required
step_size = 0.0001  # Increment size for fr

# Initialize arrays
jd_array = np.zeros(total_steps, dtype=int)
fr_array = np.zeros(total_steps)

# Populate arrays
current_fr = fr0
current_jd = jd0

for i in range(N2):
    # Check if fr exceeds a full day
    if current_fr >= 1.0:
        current_fr -= 1.0  # Reset fr, keeping the overflow
        current_jd += 1  # Move to the next Julian Day
    
    # Store values
    jd_array[i] = current_jd
    fr_array[i] = current_fr
    
    # Increment fr
    current_fr += step_size

e_over24, position_over24, velocity_over24 = aqua_new.sgp4_array(jd_array, fr_array)
position_over24_m = position_over24 * 1000

# Visualization for Problem 1
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
plot_earth(ax)

# Plot all orbits  
print("Plotting orbits")

### The plot of the magnitude of the deviation between the two positions as 
### a function of time is omitted due to an error
### in the initial position in one of the orbits which I cannot seem to find; 
### therefore, the plot in the task is not an accurate
### representation for the performance of the propagators. I have included a 
### plot of the two orbits instead
colors = plt.cm.viridis(np.linspace(0, 1, N2)) # Supposed to be a gradient colour but it doesn't really show on the plot
for i in range(N2 - 1):
    # Original orbit over 13 hours
    ax.plot(positions1[i:i+2, 0], positions1[i:i+2, 1], positions1[i:i+2, 2], color=colors[i], linewidth=2)
    # Aqua orbit propagated with SGP4 over 24 hours
    ax.plot(position_over24_m[i:i+2, 0], position_over24_m[i:i+2, 1], position_over24_m[i:i+2, 2], color='green', linewidth=2)
    # Aqua orbit propagated with RK4 over 24 hours
    ax.plot(tle_position[i:i+2, 0], tle_position[i:i+2, 1], tle_position[i:i+2, 2], color='magenta', linewidth=2)

# Mark the start and end positions
# Original
ax.scatter(*positions1[0], color='red', label='Original RK4 Start', s=100)
ax.scatter(*positions1[-1], color='blue', label='Original RK4 End', s=100)
# SGP4
ax.scatter(*position_over24_m[0], color='yellow', label='Aqua SGP4 Start', s=100)
ax.scatter(*position_over24_m[-1], color='cyan', label='Aqua SGP4 End', s=100)
# RK4
ax.scatter(*tle_position[0], color='orange', label='Aqua RK4 Start', s=100)
ax.scatter(*tle_position[-1], color='purple', label='Aqua RK4 End', s=100)

ax.set_xlabel('X Position (m)')
ax.set_ylabel('Y Position (m)')
ax.set_zlabel('Z Position (m)')
ax.set_title('Spacecraft Orbits around Earth Over 13 and 24 Hours')
plt.legend()
plt.show()




