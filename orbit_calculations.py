import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants and Initial Conditions
G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M = 5.97219e24   # mass of the Earth, kg
position0 = np.array([-5529.203, -2217.254, 3399.353]) * 1000  # Convert from km to meters
velocity0 = np.array([3.049295, 2.478617, 6.576514]) * 1000  # Convert from km/s to m/s
T = 13 * 3600.0  # Total simulation time in seconds for 13 hours

# Acceleration 
def acceleration(r):
    r_mag = np.linalg.norm(r)
    return -G * M * r / r_mag**3

# Runge-Kutta 4th Order Method for orbit calculation
def runge_kutta_orbit(pos0, vel0, N, h, T):
    positions = np.zeros((N, 3))
    velocities = np.zeros((N, 3))
    positions[0], velocities[0] = pos0, vel0

    for i in range(1, N):
        r = positions[i-1]
        v = velocities[i-1]

        k1_r = h * v
        k1_v = h * acceleration(r)

        k2_r = h * (v + 0.5 * k1_v)
        k2_v = h * acceleration(r + 0.5 * k1_r)

        k3_r = h * (v + 0.5 * k2_v)
        k3_v = h * acceleration(r + 0.5 * k2_r)

        k4_r = h * (v + k3_v)
        k4_v = h * acceleration(r + k3_r)

        positions[i] = r + (k1_r + 2*k2_r + 2*k3_r + k4_r) / 6
        velocities[i] = v + (k1_v + 2*k2_v + 2*k3_v + k4_v) / 6

    return positions, velocities, np.linalg.norm(positions[-2]), np.linalg.norm(velocities[-2])

# Plot Earth as a transparent sphere with grid outline
def plot_earth(ax, radius=6.371e6):
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x = radius * np.cos(u) * np.sin(v)
    y = radius * np.sin(u) * np.sin(v)
    z = radius * np.cos(v)
    ax.plot_wireframe(x, y, z, color="green", linewidth=0.5)

# Main simulation with two different resolutions for error estimation
N1 = 5000  # Original number of steps
h1 = T / N1  # Original step size
positions1, velocities1, final_position1, final_velocity1 = runge_kutta_orbit(position0, velocity0, N1, h1, T)

N2 = 10000  # Higher resolution with more steps
h2 = T / N2  # Smaller step size
positions2, velocities2, final_position2, final_velocity2 = runge_kutta_orbit(position0, velocity0, N2, h2, T)

# Estimate the error by comparing final positions
# Analytical error is expected to be o(h^4) for the full operation
error_estimate = np.linalg.norm(positions1[-1] - positions2[-2]) #positions2[multiplier_of_N1*(N1 - multiplier_of_N1)]
print(f"Estimated Error in Final Position using RK4 propagator between step h1 = {h1}s and h2 = {h2}s: {np.round(error_estimate,2)} meters")


## Visualization
#fig = plt.figure(figsize=(10, 8))
#ax = fig.add_subplot(111, projection='3d')
#plot_earth(ax)
#
## Plot the original orbit
#colors = plt.cm.viridis(np.linspace(0, 1, N1))
#for i in range(N1 - 1):
#    ax.plot(positions1[i:i+2, 0], positions1[i:i+2, 1], positions1[i:i+2, 2], color=colors[i], linewidth=2)
#
## Mark the start and end positions
#ax.scatter(*positions1[0], color='red', label='Start', s=100)
#ax.scatter(*positions1[-1], color='blue', label='End', s=100)
#
#ax.set_xlabel('X Position (m)')
#ax.set_ylabel('Y Position (m)')
#ax.set_zlabel('Z Position (m)')
#ax.set_title('Spacecraft Orbit around Earth Over 13 Hours')
#plt.legend()
#plt.show()

if __name__ == "__main__":
    final_positions, final_velocities, position_mag, velocity_mag = runge_kutta_orbit(position0, velocity0, N2, h2)
    print("Final Position:", final_positions[-1])
    print("Final Velocity:", final_velocities[-1])


